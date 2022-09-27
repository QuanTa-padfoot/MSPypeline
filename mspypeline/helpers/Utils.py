from typing import Optional, Dict, Tuple, Iterator, Union, Iterable, Sized, Callable
import pandas as pd
from collections import defaultdict as ddict
from itertools import combinations
from collections import deque
import numpy as np
from matplotlib.lines import Line2D
import os
from mspypeline.helpers import get_logger

logger = get_logger(str(os.path.basename(__file__).split(".")[0]))


def get_number_rows_cols_for_fig(obj: Union[int, Sized]) -> Tuple[int, int]:
    if isinstance(obj, Sized):
        obj = len(obj)
    n_rows, n_cols = 0, 0
    while n_rows * n_cols < obj:
        if n_rows <= n_cols:
            n_rows += 1
        else:
            n_cols += 1
    return n_rows, n_cols


def fill_dict(d: dict, s: str, s_split=None):
    if s_split is None:
        s_split = s.split("_")
    if len(s_split) > 1:
        fill_dict(d[s_split[0]], s, s_split[1:])
    else:
        d[s_split[0]] = s


def default_to_regular(d: ddict) -> dict:
    if isinstance(d, ddict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d


def get_analysis_design(names: Iterable[str]) -> dict:
    factory = lambda: ddict(factory)
    analysis_design = factory()
    for name in names:
        fill_dict(analysis_design, name)

    return default_to_regular(analysis_design)


def plot_annotate_line(ax, row1, row2, x, data,  fs: int = None, maxasterix: int = 3):
    """
    adjusted function
    from: https://stackoverflow.com/questions/11517986/indicating-the-statistically-significant-difference-in-bar-graph
    Annotate plot with p-values with line indicators.

    :param ax: axis of plot to put the annotaion line
    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                #text = ">" + text
                break

        if len(text) == 0:
            text = 'n. s.'

    linex = [x, x]
    liney = [row1, row2]
    midy = (row1 + row2) / 2

    ax.plot(linex, liney, c='black', linewidth=1)

    kwargs = dict(ha='center', va='center')
    if fs is not None:
        kwargs['fontsize'] = fs

    ax.text(x, midy, text, rotation=-90, **kwargs)


def venn_names(named_sets: Dict[str, set]) -> Iterator[Tuple[set, set, set]]:
    names = set(named_sets)
    for i in range(1, len(named_sets) + 1):
        for to_intersect in combinations(sorted(named_sets), i):
            others = names.difference(to_intersect)
            intersected = set.intersection(*(named_sets[k] for k in to_intersect))
            unioned = set.union(*(named_sets[k] for k in others)) if others else set()
            yield to_intersect, others, intersected - unioned


def install_r_dependencies(r_package_names, r_bioconducter_package_names):
    from rpy2.robjects.packages import importr
    import rpy2.robjects.packages as rpackages

    r_packages_uninstalled = [x for x in r_package_names if not rpackages.isinstalled(x)]
    r_bioconducter_packages_uninstalled = [x for x in r_bioconducter_package_names if not rpackages.isinstalled(x)]

    logger.debug("Installing remaining R packages: %s, remaining bioconductor packages: %s",
                 r_packages_uninstalled, r_bioconducter_packages_uninstalled)
    if r_packages_uninstalled:
        utils = importr('utils')
        utils.chooseCRANmirror(ind=1)
        for p in r_packages_uninstalled:
            utils.install_packages(p)

    if r_bioconducter_packages_uninstalled:
        biocm = importr("BiocManager")
        for p in r_bioconducter_packages_uninstalled:
            biocm.install(p)


def get_non_na_percentage(x: int) -> float:
    return 1 / (1 + np.exp(0.5 * x - 3.5)) * 0.5 + 0.5


def get_number_of_non_na_values(x: int, offset: int = 0) -> int:
    percentage = get_non_na_percentage(x)
    return max(int(np.round(percentage * x)) - offset, 3 - offset)


def get_intersection_and_unique(v1: pd.DataFrame, v2: pd.DataFrame, na_function: Callable = get_number_of_non_na_values):
    """
    Given two dataframes with identical index, determines which rows (proteins) are present in both dataframes, or
    unique to either dataframe. The number of missing values are counted row-wise and all rows above a threshold are
    marked as positive. The threshold is determined by the na_function, which determines the allowed number of missing
    values based on the number of columns. Then all rows which are positive in both dataframes are marked as
    intersection. Rows that are positive in one dataframe but are missing completely in the other are marked as unique.
    Parameters
    ----------
    v1
        first dataframe
    v2
        second dataframe
    na_function
        function which takes and returns an int

    Returns
    -------
    Three masks: The intersection, unique in v1 and unique in v2

    """
    # get number of allowed non na values for both dataframes
    non_na_group_1 = na_function(v1.shape[1])
    non_na_group_2 = na_function(v2.shape[1])
    # find rows which fulfill the requirement
    mask_1 = (v1 > 0).sum(axis=1) >= non_na_group_1
    mask_2 = (v2 > 0).sum(axis=1) >= non_na_group_2
    # combine both rows
    mask = np.logical_and(mask_1, mask_2)
    # determine missing
    missing_1 = (v1 > 0).sum(axis=1) == 0
    missing_2 = (v2 > 0).sum(axis=1) == 0
    # determine exclusive
    exclusive_1 = np.logical_and(mask_1, missing_2)
    exclusive_2 = np.logical_and(mask_2, missing_1)
    return mask, exclusive_1, exclusive_2


def dict_depth(d: dict) -> int:
    level = 0
    queue = deque([(id(d), d, level)])
    memo = set()
    while queue:
        id_, o, level = queue.popleft()
        if id_ in memo:
            continue
        memo.add(id_)
        if isinstance(o, dict):
            queue += ((id(v), v, level + 1) for v in o.values())
    return level


def get_legend_elements(labels: list, color_map: Optional[dict] = None, marker_size:Optional[int] = None):
    """
        Returns custom legend elements based on a list of labels and an optional color map.
        These elements can be passed to a legend via the 'handles' parameter

        Parameters
        ----------
        labels:
            list of strings
        color_map:
            dict of strings, with keys being the name of a label and values the corresponding color

    """
    if color_map is None:
        color_map = {name: f"C{i}" for i, name in enumerate(labels)}
    if marker_size is None:
        marker_size=10
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=name.replace("_", " "),
                              markerfacecolor=color_map.get(name, "blue"), markersize=marker_size)
                       for name in labels]
    return legend_elements


def get_plot_name_suffix(df_to_use: Optional[str] = None, level: Optional[int] = None) -> str:
    """
    Generate a suffix for the plot name

    Parameters
    ----------
    df_to_use
        dataframe that was used
    level
        level on which data was aggregated

    Returns
    -------
    a string which can be used as a suffix for file paths

    """
    s = "" if df_to_use is None else f"_{df_to_use}"
    s += "" if level is None else f"_level_{level}"
    return s


class DataDict(dict):
    """
    Overwrites the standard dictionary to provide an additional DataSource.
    When a missing key is looked up the DataSource is searched for a method named:
    e.g. looking up key=parameters, looking for method named "preprocess_parameters",
    which is expected to return data, which will then be stored under the key.
    This allows data from disk to be loaded on demand instead of loading all possible data at the beginning.
    """
    def __init__(self, data_source, *args, **kwargs):
        """
        Parameters
        ----------
        data_source
            class which will be searched for methods
        args
            passed to dict.__init__
        kwargs
            passed to dict.__init__
        """
        super().__init__(*args, **kwargs)
        self.data_source = data_source

    def __missing__(self, key):
        try:
            self.data_source.logger.debug("Reading %s from disk", key)
            data = getattr(self.data_source, f"preprocess_{key}")()
            self[key] = data
            return data
        except FileNotFoundError as e:
            raise KeyError("Missing file:", key, e)
        except AttributeError as e:
            raise KeyError("Missing function to load:", key, e)


def format_docstrings(**mapping):
    def docstring_decorator(fn):
        fn.__doc__ = fn.__doc__.format(**mapping)
        return fn
    return docstring_decorator


def add_end_docstrings(*docstr):
    def docstring_decorator(fn):
        fn.__doc__ = fn.__doc__ + "".join(docstr)
        return fn
    return docstring_decorator


def make_contrasts(exp1, exp2):
    contrast_matrix = pd.DataFrame(data=(1, -1), columns=[f"{exp2}-{exp1}"], index=[exp2, exp1])
    contrast_matrix.name = "Contrasts"
    return contrast_matrix
