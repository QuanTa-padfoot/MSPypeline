import os
from copy import copy
from itertools import combinations
from typing import Tuple, Optional, Union, Callable, Dict, Iterable
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.lines as mlines
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection
from adjustText import adjust_text
from matplotlib.colorbar import ColorbarBase
from matplotlib_venn import venn2, venn3
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde
from scipy import stats
from sklearn.decomposition import PCA
import functools
import warnings

from mspypeline.helpers import get_number_rows_cols_for_fig, plot_annotate_line, get_legend_elements, \
    get_plot_name_suffix, venn_names, format_docstrings

FIG_FORMAT = ".pdf"
rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42

TECHREP_SUFFIX = ", \n technical replicates aggregated"


def linear(x, m, b):
    return m * x + b


def collect_plots_to_pdf(path: str, *args, dpi: int = 200):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if not path.endswith(".pdf"):
        path += ".pdf"
    with PdfPages(path) as pdf:
        for plot in args:
            figure = None
            if callable(plot):
                figure, axes = plot()
            elif isinstance(plot, Iterable):
                for x in plot:
                    if isinstance(x, plt.Figure):
                        figure = x
                        break
            elif isinstance(plot, plt.Figure):
                figure = plot
            if figure is not None:
                pdf.savefig(figure, dpi=dpi)
                plt.close(figure)


_get_path_and_name_kwargs_doc = """
        Uses following kwargs to generate a path and file name.
        
        * save_path (str): directory where file should be saved
        * df_to_use (str): used to generate a name suffix
        * level (int): used to generate a name suffix
        * split_files (bool): should the file be saved in a subdirectory

        Additionally, all parts of the passed name enclosed by curly brackets will be replaced from passed kwargs.
"""


@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def get_path_and_name_from_kwargs(name: str, **kwargs) -> Tuple[str, str]:
    """
    creates a the path and file name for a plot

    Parameters
    ----------
    name
    kwargs
        {kwargs}

    Returns
    -------
    tuple(str, str)
        The path where the file will be saved and a name for the file without a file type.
    """
    save_path = kwargs.get("save_path", None)
    # get values with default
    # suffix
    df_to_use = kwargs.get("df_to_use", None)
    level = kwargs.get("level", None)

    # get path modifications from the name
    prefix, name = os.path.split(name)
    # modify the save path if specified
    split_files = kwargs.get("split_files", True)
    if save_path is not None and split_files:
        save_path = os.path.join(save_path, prefix)

    suffix = get_plot_name_suffix(df_to_use=df_to_use, level=level)
    name = f"{name}{suffix}"
    name = name.format_map(kwargs)
    return save_path, name


def save_plot_func(
        fig: plt.Figure, path: str, plot_name: str, func: Callable, fig_format: str = FIG_FORMAT,
        dpi: int = 200, **kwargs
) -> None:
    """
    Saves figure in path. Directories will be created if they not exist.

    Parameters
    ----------
    fig
        figure to be saved
    path
        path to the plot
    plot_name
        name of the saved figure
    func
        function used to save the plots
    fig_format
        figure format of the plot. default is PDF.
    dpi
        DPI of saved figure
    kwargs
        accepts kwargs

    """
    if path is not None:
        try:
            os.makedirs(path, exist_ok=True)
            res_path = os.path.join(path, plot_name)
            fig.savefig(res_path + fig_format, dpi=dpi, bbox_inches="tight")
        except PermissionError:
            warnings.warn(f"Permission error in function {str(func).split(' ')[1]}. Did you forget to close the file?")


def save_plot(plot_name: str):
    """
    Decorator to save figures, which are returned by the decorated function. Assumes that a tuple of figure, axes is
    returned. Plot is saved by using get_path_and_name_from_kwargs and save_plot_func.

    Parameters
    ----------
    plot_name
        string to be saved to. Sting can contain "/" to indicate a folder structure where the plot
        should be saved. Also can contain preformatted parts like: "plot_{name}". The {name} will then be replaced
        by a passed kwarg "name".

    """
    def decorator_save_plot(func):
        func.__doc__ = func.__doc__.format(name=plot_name)

        @functools.wraps(func)
        def wrapper_save_plot(*args, **kwargs):
            # run original function
            ret = func(*args, **kwargs)
            if ret is not None:
                path, pn = get_path_and_name_from_kwargs(name=plot_name, **kwargs)
                save_plot_func(ret[0], path, pn, func, **kwargs)
            return ret
        return wrapper_save_plot
    return decorator_save_plot


def save_csv_fn(save_path: str, csv_name: str, df: Union[pd.Series, pd.DataFrame]):
    if save_path is not None:
        os.makedirs(save_path, exist_ok=True)
        df.to_csv(os.path.join(save_path, csv_name) + ".csv", header=True)


def save_csvs(name_map: Dict[str, str]):
    """
    Saves all dataframes as csv which are specified as dict keys. The values are the file names.

    Parameters
    ----------
    name_map
        mapping of kwarg name to file name

    """
    def decorator_save_csvs(func):
        @functools.wraps(func)
        def wrapper_save_csvs(*args, **kwargs):
            for kwarg_name, file_name in name_map.items():
                df = kwargs.get(kwarg_name, None)
                if df is not None:
                    save_path, csv_name = get_path_and_name_from_kwargs(file_name, **kwargs)
                    save_csv_fn(save_path, csv_name, df)
            return func(*args, **kwargs)
        return wrapper_save_csvs
    return decorator_save_csvs


def save_venn_to_txt(name_map: Dict[str, str]):
    def decorator_save_venn(func):
        @functools.wraps(func)
        def wrapper_save_venn(*args, **kwargs):
            for kwarg_name, file_name in name_map.items():
                named_sets = kwargs.get(kwarg_name, None)
                if named_sets is not None:
                    if len(named_sets) > 6:
                        warnings.warn(f"Skipping save_venn_to_txt because more than 6 experiments were passed at once. "
                                      f"({len(named_sets)}")
                        continue
                    save_path, txt_name = get_path_and_name_from_kwargs(file_name, **kwargs)
                    if save_path is None:
                        continue
                    os.makedirs(save_path, exist_ok=True)
                    for intersected, unioned, result in venn_names(named_sets):
                        # create name based on the intersections and unions that were done
                        intersected_name = "&".join(sorted(intersected))
                        unioned_name = "-" + "-".join(sorted(unioned)) if unioned else ""
                        res_path = os.path.join(save_path, f"{txt_name}_{intersected_name}{unioned_name}.txt")
                        # write all names line by line into the file
                        with open(res_path, "w") as out:
                            for re in result:
                                out.write(re + "\n")
            return func(*args, **kwargs)
        return wrapper_save_venn
    return decorator_save_venn


class QuantileNormalize(colors.Normalize):
    """
    Can be used to color the values of a dataset according to the quantiles.
    """
    def __init__(self, quantiles: pd.Series):
        """
        Parameters
        ----------
        quantiles
            calculated quantiles of the dataset
        """
        self.quantiles = quantiles
        assert self.quantiles.index.min() >= 0
        assert self.quantiles.index.max() <= 1
        self.outputs = np.concatenate((self.quantiles.index.values, [1]))
        super().__init__(0, 1)

    def __call__(self, value, clip=None):
        index = np.searchsorted(self.quantiles, value)
        return np.ma.array(self.outputs[index], fill_value=0)


def split_plot(n_rows, n_cols, figsize=(7, 7), plot_name="", data=None, plot_fn=None, save_path=None):
    n_figures = int(np.ceil(len(data.columns) / (n_rows * n_cols)))

    if save_path is not None:
        with PdfPages(os.path.join(save_path, plot_name + ".pdf")) as pdf:
            for n_figure in range(n_figures):
                fig, axarr = plt.subplots(n_rows, n_cols, figsize=figsize)
                for i, (pos, ax) in enumerate(np.ndenumerate(axarr)):
                    idx = n_figure * (n_rows * n_cols) + i
                    try:
                        experiment = data.columns[idx]
                    except IndexError:
                        break
                    plot_fn(ax, data[experiment], experiment)
                fig.tight_layout(rect=[0, 0.03, 1, 0.95])

                pdf.savefig(fig)
                plt.close(fig)


@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
@save_csvs({"unique_g1": "csv_significant/volcano_plot_data_{g1}_vs_{g2}_unique_{g1}",
            "unique_g2": "csv_significant/volcano_plot_data_{g1}_vs_{g2}_unique_{g2}"})
def save_volcano_results(
        volcano_data: pd.DataFrame, interesting_proteins, unique_g1: pd.Series = None,
        unique_g2: pd.Series = None, g1: str = "group1", g2: str = "group2", adj_pval: bool = False,
        intensity_label: str = "Intensity", show_suptitle: bool = True, pval_threshold: float = 0.05,
        fchange_threshold: float = 2, scatter_size: float = 20, n_labelled_proteins: int = 10, close_plots: str = "all",
        exp_has_techrep: bool = False, **kwargs
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes, plt.Axes]]:
    """
    Saves multiple csv files and images containing the information of the volcano plot

    Parameters
    ----------
    volcano_data
        DataFrame containing data for the volcano plot with columns logFC and column specified under col. The index
        should be protein names or gene names
    interesting_proteins
        mapping of pathways that shoul be annotated in the volcano plot
    unique_g1
        Series containing intensities of proteins or genes unique to group one
    unique_g2
        Series containing intensities of proteins or genes unique to group two
    g1
        name of first sample that should be analysed (downregulated)
    g2
        name of second sample that should be analysed (upregulated)
    adj_pval
        should adjusted p values be used
    intensity_label
        from which intensities were the fold changes calculated
    show_suptitle
        should the figure title be shown
    pval_threshold
        maximum p value to be considered significant
    fchange_threshold
        minimum fold change threshold (before log2 transformation) to be labelled significant
    scatter_size
        size of the points in the scatter plots
    n_labelled_proteins
        number of points that will be annotated in th plot
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    if close_plots is not None:
        plt.close(close_plots)

    col_mapping = {"adjpval": "adjusted p value", "pval": "unadjusted p value"}
    if adj_pval:
        col = "adjpval"
    else:
        col = "pval"

    # given pathway mapping(s)/proteins of interest create a list later used to be annotated in the volcano plot
    pathway_label = ""
    POI_vals = []
    if interesting_proteins.values():
        all_pathway_proteins = set.union(*(set(x) for x in interesting_proteins.values()))
        POI_vals = volcano_data[volcano_data.index.isin(list(all_pathway_proteins))]
        pathway_label = "_".join(interesting_proteins.keys())

    def get_volcano_significances(fchange, pval, pval_threshold, fchange_threshold):
        if pval > pval_threshold or abs(fchange) < np.log2(fchange_threshold):
            return "ns"
        elif fchange >= 0:
            return "up"
        elif fchange < 0:
            return "down"
        else:
            raise ValueError(f"heisenbug: fold change: {fchange}, p value: {pval}")

    g1_name = g1.replace("_", " ")
    g2_name = g2.replace("_", " ")

    # add the measured regulation to the data based on the given thresholds
    volcano_data["regulation"] = [get_volcano_significances(log_fold_change, p_val, pval_threshold, fchange_threshold)
                                  for log_fold_change, p_val in zip(volcano_data["logFC"], volcano_data[col])]

    # save the volcano data csv in full and only the significant part
    save_path, csv_name = get_path_and_name_from_kwargs(
        "csv_full/volcano_plot_data_{g1}_vs_{g2}_full_{p}_{pathway_label}", g1=g1, g2=g2,
        p=col_mapping[col].replace(' ', '_'), pathway_label=pathway_label, **kwargs)
    save_csv_fn(save_path, csv_name, volcano_data)
    save_path, csv_name = get_path_and_name_from_kwargs(
        "csv_significant/volcano_plot_data_{g1}_vs_{g2}_significant_{p}_{pathway_label}", g1=g1, g2=g2,
        p=col_mapping[col].replace(' ', '_'), pathway_label=pathway_label, **kwargs)
    save_csv_fn(save_path, csv_name, volcano_data[volcano_data[col] < 0.05])

    significance_to_color = {"down": "blue", "ns": "gray", "up": "red"}
    significance_to_label = {"down": f"higher in  {g1_name}", "ns": "non-significant", "up": f"higher in  {g2_name}"}

    # plot
    fig = plt.figure(figsize=(7, 7))

    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 8, 1])
    ax_unique_down: plt.Axes = plt.subplot(gs[0])
    ax: plt.Axes = plt.subplot(gs[1])
    ax_unique_up: plt.Axes = plt.subplot(gs[2])

    # hide the spines between ax and ax2
    ax_unique_down.spines['right'].set_visible(False)
    ax_unique_up.spines['left'].set_visible(False)
    ax_unique_down.yaxis.tick_left()
    ax_unique_up.yaxis.tick_right()
    ax_unique_up.yaxis.set_label_position("right")
    # hide the xticks
    ax_unique_down.tick_params(which='both', bottom=False, labelbottom=False)
    ax_unique_up.tick_params(which='both', bottom=False, labelbottom=False)

    # non sign gray, left side significant blue, right side red
    for regulation in significance_to_color:
        mask = [x == regulation for x in volcano_data["regulation"]]
        ax.scatter(volcano_data["logFC"][mask], -np.log10(volcano_data[col])[mask], s=scatter_size, edgecolors="none",
                   color=significance_to_color[regulation], alpha=0.6,
                   label=f"{sum(mask)} {significance_to_label[regulation]}")
    # get axis bounds for vertical and horizontal lines
    ymin, ymax = ax.get_ybound()
    xmin, xmax = ax.get_xbound()
    m = max(abs(xmin), xmax)
    # center the plot around 0
    ax.set_xlim(left=-1 * m, right=m)
    # update the x bounds
    xmin, xmax = ax.get_xbound()
    axline_kwargs = dict(linestyle="--", color="black", alpha=0.5, linewidth=1)
    # add line at significance threshold
    if any(volcano_data[col] < 0.05):
        x_offset = (np.log2(fchange_threshold) / xmax) / 2
        ax.axhline(-np.log10(0.05), **axline_kwargs, xmin=0, xmax=0.5 - x_offset)
        ax.axhline(-np.log10(0.05), **axline_kwargs, xmin=0.5 + x_offset, xmax=1)

    # add lines for minimum fold change threshold
    y_percentage = (-np.log10(0.05) + abs(ymin)) / (ymax + abs(ymin))
    if fchange_threshold > 0:
        ax.axvline(-np.log2(fchange_threshold), **axline_kwargs, ymin=y_percentage, ymax=1)
        ax.axvline(np.log2(fchange_threshold), **axline_kwargs, ymin=y_percentage, ymax=1)
    # plot unique values with mean intensity at over maximum
    ax_unique_down.scatter([0] * len(unique_g1), unique_g1, s=scatter_size, color="dodgerblue", edgecolors="none",
                           alpha=0.6, label=f"{len(unique_g1)} unique in  {g1_name}")
    ax_unique_up.scatter([0] * len(unique_g2), unique_g2, s=scatter_size, color="coral", edgecolors="none",
                         alpha=0.6, label=f"{len(unique_g2)} unique in  {g2_name}")
    # adjust bounds for unique axis
    ymin_down, ymax_down = ax_unique_down.get_ybound()
    ymin_up, ymax_up = ax_unique_up.get_ybound()
    ax_unique_down.set_ylim(bottom=min(ymin_down, ymin_up), top=max(ymax_down, ymax_up))
    ax_unique_up.set_ylim(bottom=min(ymin_down, ymin_up), top=max(ymax_down, ymax_up))

    # figure stuff
    if show_suptitle:
        fig.suptitle(f"{g1_name} vs {g2_name}" + (TECHREP_SUFFIX if exp_has_techrep else ""))
    ax.set_xlabel(f"{intensity_label} Fold Change")
    ax.set_ylabel(r"-$Log_{10}$" + f" {col_mapping[col]}")
    ax_unique_down.set_ylabel(intensity_label)
    ax_unique_up.set_ylabel(intensity_label)
    fig.legend(bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=False)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    # save intermediate results
    path, plot_name = get_path_and_name_from_kwargs(name="plots/volcano_{g1}_{g2}_no_annotation_{p}_{pathway_label}",
                                                    g1=g1, g2=g2, p=col_mapping[col].replace(' ', '_'),
                                                    pathway_label=pathway_label, **kwargs)
    save_plot_func(fig, path, plot_name, save_volcano_results, **kwargs)

    # add text labels to the most significantly regulated genes
    significant_upregulated = volcano_data[
        (volcano_data["logFC"] > np.log2(fchange_threshold)) & (volcano_data[col] < 0.05)
    ].sort_values(by=[col], ascending=True).head(n_labelled_proteins)
    significant_downregulated = volcano_data[
        (volcano_data["logFC"] < -np.log2(fchange_threshold)) & (volcano_data[col] < 0.05)
    ].sort_values(by=[col], ascending=True).head(n_labelled_proteins)
    significant = pd.concat([significant_upregulated, significant_downregulated])

    texts = []
    # if list for proteins of interest (pathway list chosen) is given, annotate those proteins,
    # otherwise annotate most significant proteins
    if len(POI_vals) > 0:
        for log_fold_change, p_val, gene_name in zip(POI_vals["logFC"], POI_vals[col], POI_vals.index):
            texts.append(ax.text(log_fold_change, -np.log10(p_val), gene_name, ha="center", va="center", fontsize=8,
                                 color="black"))
    else:
        for log_fold_change, p_val, gene_name in zip(significant["logFC"], significant[col], significant.index):
            texts.append(ax.text(log_fold_change, -np.log10(p_val), gene_name, ha="center", va="center", fontsize=8,
                                 color="black"))
    adjust_text(texts, arrowprops=dict(width=0.15, headwidth=0, color='gray', alpha=0.6), ax=ax)

    # save the final result
    path, plot_name = get_path_and_name_from_kwargs(name="plots/volcano_{g1}_{g2}_annotation_{p}_{pathway_label}",
                                                    g1=g1, g2=g2, p=col_mapping[col].replace(' ', '_'),
                                                    pathway_label=pathway_label, **kwargs)
    save_plot_func(fig, path, plot_name, save_volcano_results, **kwargs)
    # TODO scatter plot of significant genes
    return fig, (ax, ax_unique_down, ax_unique_up)


@save_plot("pca_overview")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_pca_results(
        pca_data: pd.DataFrame, pca_fit: PCA = None, normalize: bool = True, intensity_label: str = "Intensity",
        color_map: Optional[dict] = None, show_suptitle: bool = True, marker_size: int = 150,
        legend_marker_size: int = 12, close_plots: str = "all", exp_has_techrep: bool = False, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Saves image containing the pca results with prefix: {{name}}

    Parameters
    ----------
    pca_data
        DataFrame containing transformed/dimensionally-reduced data with which PCA was performed
    pca_fit
        PCA object that was fitted to normalized input data
    normalize
        should the transformed data be normalized with the singular values before plotting
    intensity_label
        figure title
    color_map
        mapping from column name to color if custom colors are wanted
    show_suptitle:
        should the figure title be shown
    marker_size
        size of the points in the scatter plots
    legend_marker_size
        size of the legend marker
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    if close_plots is not None:
        plt.close(close_plots)
    n_components = pca_data.shape[0]
    singular_values = np.ones(n_components)
    base_color_map = {value: f"C{i}" for i, value in enumerate(pca_data.columns.get_level_values(0).unique())}
    color_map = {} if color_map is None else color_map
    base_color_map.update(color_map)
    if normalize and pca_fit is None:
        warnings.warn("Normalizing not possible when pca_fit is None")
    elif normalize and pca_fit is not None:
        singular_values = pca_fit.singular_values_

    if pca_fit is not None:
        pca_var = np.round(pca_fit.explained_variance_ratio_ * 100, decimals=1)
        pca_var = [f"{var} %" for var in pca_var]
    else:
        pca_var = ["nan" for _ in range(len(pca_fit.n_components_))]

    if n_components == 2:
        fig, axarr = plt.subplots(1, 1, figsize=(14, 14))
        ax = axarr
        ax.scatter(
            pca_data.loc["PC_1"] / singular_values[0],
            pca_data.loc["PC_2"] / singular_values[1],
            s=marker_size, edgecolors="none",
            c=[base_color_map.get(name, "blue") for name in pca_data.columns.get_level_values(0)])
        ax.set_xlabel("PC 1 - {0}".format(pca_var[0]), fontsize=22, labelpad=20)
        ax.set_ylabel("PC 2 - {0}".format(pca_var[1]), fontsize=22, labelpad=20)
        ax.tick_params(axis="both", labelsize=18)
    else:
        fig, axarr = plt.subplots(n_components, n_components, figsize=(14, 14))
        for row in range(n_components):
            row_pc = row + 1
            for col in range(n_components):
                col_pc = col + 1
                if row > col:
                    ax = axarr[col, row]
                    ax.scatter(
                        pca_data.loc[f"PC_{row_pc}"] / singular_values[row],
                        pca_data.loc[f"PC_{col_pc}"] / singular_values[col],
                        s=marker_size, edgecolors=None,
                        c=[base_color_map.get(name, "blue") for name in pca_data.columns.get_level_values(0)])
                    ax.set_xlabel(f"PC {row_pc} - {pca_var[0]}", fontsize=22, labelpad=20)
                    ax.set_ylabel(f"PC {col_pc} - {pca_var[1]}", fontsize=22, labelpad=20)
                    ax.tick_params(axis="both", labelsize=18)

    if show_suptitle:
        fig.suptitle(intensity_label + (TECHREP_SUFFIX if exp_has_techrep else ""), fontsize=30)
    legend_elements = get_legend_elements(labels=pca_data.columns.get_level_values(0).unique(),
                                          color_map=base_color_map, marker_size=legend_marker_size)
    fig.legend(handles=legend_elements, bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=False, fontsize=20)
    fig.tight_layout(rect=[0.03, 0.03, 1, 0.95])
    return fig, axarr


@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
@save_csvs({"protein_intensities": "csv_intensity/{pathway}_protein_intensities",
            "significances": "csv_pval/{pathway}_pvalues"})
def save_pathway_analysis_results(
        protein_intensities: pd.DataFrame, significances: pd.DataFrame = None, pathway: str = "",
        show_suptitle: bool = True, threshold: float = 0.05, intensity_label: str = "Intensity",
        color_map: Optional[dict] = None, close_plots: str = "all", exp_has_techrep: bool = False, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Saves plots into the pathway_analysis dir.
    
    Parameters
    ----------
    protein_intensities
        data of intensities
    significances
        significances between different conditions
    pathway
        name of the pathway
    show_suptitle
        should the pathway name be shown as figure title
    threshold
        maximum p value indicating significance
    intensity_label
        from which intensity was the data. will be shown on x axis
    color_map
        a mapping from the column names to a color
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    Returns
    -------

    """
    if close_plots is not None:
        plt.close(close_plots)
    level_keys = list(protein_intensities.columns.get_level_values(0).unique())
    n_rows, n_cols = get_number_rows_cols_for_fig(protein_intensities.index)
    fig, axarr = plt.subplots(n_rows, n_cols, figsize=(n_cols * 5, int(n_rows * len(level_keys) / 1.1)))
    for i in range(n_rows * n_cols - len(protein_intensities.index)):
        axarr[n_rows - 1, n_cols - 1 - i].remove()
    result_color_map = {value: f"C{i}" for i, value in enumerate(level_keys)}
    result_color_map.update(color_map if color_map is not None else {})
    if show_suptitle:
        fig.suptitle(pathway.replace("_", " ") + (TECHREP_SUFFIX if exp_has_techrep else ""), size=26)
    for protein, (pos, ax) in zip(protein_intensities.index, np.ndenumerate(axarr)):
        ax.scatter(protein_intensities.loc[protein],
                   [level_keys.index(c) for c in protein_intensities.columns.get_level_values(0)],
                   c=[result_color_map[c] for c in protein_intensities.columns.get_level_values(0)], edgecolors="none",
                   alpha=0.7)
        ax.set_title(protein)
        ax.set_ylim((-1, len(level_keys)))
        ax.set_yticks([i for i in range(len(level_keys))])
        level_keys_labels = [key.replace("_", " ") for key in level_keys]
        if len(level_keys_labels) == 0:
            level_keys_labels = level_keys
        ax.set_yticklabels(level_keys_labels)
        ax.set_xlabel(intensity_label)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    path, plot_name = get_path_and_name_from_kwargs(name="plots/{pathway}_no_labels", pathway=pathway, **kwargs)
    save_plot_func(fig, path, plot_name, save_pathway_analysis_results, **kwargs)

    if significances is not None:
        for protein, (pos, ax) in zip(protein_intensities.index, np.ndenumerate(axarr)):
            # adjust axis height based on number of significant differences
            to_annotate = significances.loc[protein]
            to_annotate = to_annotate[to_annotate <= threshold]
            xmin, xmax = ax.get_xbound()
            ax.set_xlim(right=xmax * (1 + to_annotate.shape[0] * 0.015))
            for i, (index, pval) in enumerate(to_annotate.items()):
                plot_annotate_line(ax, level_keys.index(index[0]), level_keys.index(index[1]),
                                   xmax * (1 + i * 0.015) - 0.005, pval, maxasterix=3)

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        path, plot_name = get_path_and_name_from_kwargs(name="plots/{pathway}", pathway=pathway, **kwargs)
        save_plot_func(fig, path, plot_name, save_pathway_analysis_results, **kwargs)
    return fig, axarr


@save_plot("boxplot")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_boxplot_results(
        protein_intensities: pd.DataFrame, intensity_label: str = "Intensity",
        plot: Optional[Tuple[plt.Figure, plt.Axes]] = None, vertical: bool = False, close_plots: str = "all", **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Boxplot of intensities. Saves the plot with prefix: {{name}}

    Parameters
    ----------
    protein_intensities
        DataFrame where each column are the intensities to boxplot, column names will be used as labels
    intensity_label
        label of the x axis of the plot
    plot
        figure to put plot
    vertical
        should a vertical boxplot be created
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    kwargs
        {kwargs}

    """
    # TODO give colors to the different groups
    if plot is None:
        if close_plots is not None:
            plt.close(close_plots)
        fig, ax = plt.subplots(figsize=(14, 1 + len(protein_intensities.columns) // 2))
    else:
        fig, ax = plot
    # indicate overall median with a line
    median_value = np.nanmedian(protein_intensities.values.flatten())
    line_kwargs = dict(color="black", alpha=0.5, linewidth=1)
    if vertical:
        ax.axhline(median_value, **line_kwargs)
    else:
        ax.axvline(median_value, **line_kwargs)
    # convert the data into a list of lists and filter nan values
    data = [
        protein_intensities.loc[~pd.isna(protein_intensities.loc[:, c]), c].tolist()
        for c in protein_intensities.columns
    ]
    labels = [label.replace("_", " ") for label in protein_intensities.columns]
    ax.boxplot(data, vert=vertical, labels=labels)
    if vertical:
        ax.set_ylabel(intensity_label, fontsize=15)
    else:
        ax.set_xlabel(intensity_label, fontsize=15)
    ax.tick_params(labelsize=15)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax


@save_plot("rel_std_{experiment_name}")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_relative_std_results(
        intensities: pd.DataFrame, experiment_name: str, intensity_label: str = "Intensity",
        show_suptitle: bool = True, bins=(10, 20, 30), cmap: dict = None, close_plots: str = "all",
        exp_has_techrep: bool = False, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Relative standard deviations of passed intensities with color marking based on the specified bins and color map.
    Save the plot with prefix: {{name}}

    Parameters
    ----------
    intensities
        DataFrame with experiment intensities to be plotted
    experiment_name
        name of the overall experiment
    intensity_label
        name of the intensities for the x label
    show_suptitle
        should the figure title be shown
    bins
        in which bins should the standard deviations be categorized
    cmap
        mapping for the digitized labels to a color
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    # TODO add percentage to absolute numbers
    # TODO see if code could be optimized
    if close_plots is not None:
        plt.close(close_plots)

    bins = np.array(bins)

    default_cm = {0: "navy", 1: "royalblue", 2: "skyblue", 3: "darkgray"}
    if cmap is not None:
        default_cm.update(cmap)

    if "Log_2" in intensity_label:
        relative_std_percent = intensities.std(axis=1) / intensities.mean(axis=1) * 100
    else:
        intensities = np.log2(intensities)
        relative_std_percent = intensities.std(axis=1) / intensities.mean(axis=1) * 100

    inds = np.digitize(relative_std_percent, bins).astype(int)
    plot_colors = pd.Series([default_cm.get(x, "black") for x in inds], index=relative_std_percent.index)
    color_counts = {color: (plot_colors == color).sum() for color in plot_colors.unique()}

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    ax.scatter(intensities.mean(axis=1), relative_std_percent, c=plot_colors, marker="o", s=(4 * 72. / fig.dpi) ** 2,
               alpha=0.7, edgecolors="none")

    experiment_name = experiment_name.replace("_", "  ")
    if show_suptitle:
        fig.suptitle(experiment_name + (TECHREP_SUFFIX if exp_has_techrep else ""))
    ax.set_xlabel(f"Mean {intensity_label}")
    ax.set_ylabel("Relative Standard deviation [%]")
    if "Log_2" not in intensity_label:
        ax.set_xscale('log')
    xmin, xmax = ax.get_xbound()
    cumulative_count = 0
    for i, bin_ in enumerate(bins):
        cumulative_count += color_counts.get(default_cm[i], 0)
        ax.axhline(bin_, color=default_cm[i])
        ax.text(xmin, bin_, cumulative_count)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax


@save_plot("detected_counts")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_detection_counts_results(
        counts: pd.DataFrame, intensity_label: str = "Intensity", show_suptitle: bool = True, close_plots: str = "all",
        exp_has_techrep: bool = False, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Saves the plot with prefix: {{name}}

    Parameters
    ----------
    counts
        DataFrame containing the counts to be plotted
    intensity_label
        label of the dataframe
    show_suptitle
        should the figure title be shown
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    if close_plots is not None:
        plt.close(close_plots)

    n_rows_experiment, n_cols_experiment = get_number_rows_cols_for_fig(counts.columns)
    fig, axarr = plt.subplots(n_rows_experiment, n_cols_experiment, squeeze=True,
                              figsize=(5 * n_cols_experiment, 3 * n_rows_experiment))
    for i in range(n_rows_experiment * n_cols_experiment - len(counts.columns)):
        axarr[n_rows_experiment - 1, n_cols_experiment - 1 - i].remove()

    if show_suptitle:
        fig.suptitle(f"Detection counts from {intensity_label}" + (TECHREP_SUFFIX if exp_has_techrep else ""))

    global_max = counts.max().max()
    for (pos, ax), col in zip(np.ndenumerate(axarr), counts.columns):
        col_data = counts.loc[:, col]
        col_data = col_data[~pd.isna(col_data)]
        col_name = col.replace("_", " ")

        ax.set_title(f"{col_name} \n total detected: {int(col_data.sum())}")
        ax.barh(col_data.index, col_data, color="skyblue")

        if len(col_data) in range(1, 9):
            fsize = 11
        elif len(col_data) in range(9, 13):
            fsize = 7
        else:
            fsize = 5

        for y, value in zip(col_data.index, col_data):
            ax.text(col_data.max() / 2, y, value,
                    verticalalignment='center', horizontalalignment='center', fontsize=fsize)

        ax.set_yticks(col_data.index)
        ax.set_yticklabels([f"detected in {i} replicates" for i in col_data.index], fontsize=fsize)
        ax.set_xlim(0, global_max + global_max*0.05)
        ax.set_xlabel("Counts", fontsize=fsize)
        ax.tick_params(axis="both", labelsize=fsize)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, axarr


@save_plot("{source}_correlation_heatmap")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_correlation_heatmap_results(
        correlations: pd.DataFrame, intensity_label: str = "Intensity", show_suptitle: bool = True,
        close_plots: str = "all", exp_has_techrep: bool = False, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Saves the plot with prefix: {{name}}

    Parameters
    ----------
    correlations
        DataFrame containing the correlations to be plotted
    intensity_label
        label of the dataframe
    show_suptitle
        should the figure title be shown
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    # TODO save csv
    if close_plots is not None:
        plt.close(close_plots)

    num_cols, a = correlations.shape
    assert num_cols == a, "Dataframe needs to be quadratic"

    mask = np.zeros_like(correlations).astype(bool)
    mask[np.triu_indices_from(mask)] = True

    wid_hei = 4 + 0.5 * num_cols
    fig, ax = plt.subplots(1, 1, figsize=(wid_hei, wid_hei))

    if show_suptitle:
        fig.suptitle(f"Correlation Heatmap {intensity_label}" + (TECHREP_SUFFIX if exp_has_techrep else ""))

    mesh = ax.pcolormesh(np.ma.masked_where(mask, correlations.values), cmap="coolwarm")

    ax.figure.colorbar(mesh, ax=ax)
    ax.invert_yaxis()

    # set x and y ticks
    ax.set_xticks(np.linspace(start=0.5, stop=num_cols - 0.5, num=num_cols))
    ax.set_xticklabels(correlations.index, rotation=90)
    ax.set_yticks(np.linspace(start=0.5, stop=num_cols - 0.5, num=num_cols))
    ax.set_yticklabels(correlations.index)

    # annotate values
    for x, col in enumerate(correlations.columns):
        for y, idx in enumerate(correlations.index):
            if not mask[y, x]:
                ax.text(x + 0.5, y + 0.5, f"{correlations.loc[idx, col]:.4f}", ha="center", va="center")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax


@save_plot("kde")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_kde_results(
        intensities: pd.DataFrame, quantile_range: Optional[np.array] = None, n_points: int = 1000,
        cmap: Union[str, colors.Colormap] = "viridis", plot: Optional[Tuple[plt.Figure, plt.Axes]] = None,
        intensity_label: str = "Intensity", close_plots: str = "all", **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    saves the plot with prefix: {{name}}
    
    Parameters
    ----------
    intensities
        protein intensities to be plotted
    quantile_range
        default is np.arange(0.05, 1, 0.05)
    n_points
        number of points to sample from distribution
    cmap
        color map to use
    plot
        figure to put plot
    intensity_label
        label to be put on x axis
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    kwargs
        {kwargs}

    Returns
    -------

    """
    if plot is not None:
        fig, ax = plot
    else:
        if close_plots is not None:
            plt.close(close_plots)
        fig, ax = plt.subplots(1, 1)

    if quantile_range is None:
        quantile_range = np.arange(0.05, 1, 0.05)

    for col in intensities.columns:
        intensity_quantiles = intensities.loc[:, col].quantile(quantile_range)
        kde_fit = gaussian_kde(intensities.loc[~pd.isna(intensities.loc[:, col]), col])

        x = np.linspace(intensities.loc[:, col].min() * 0.9, intensities.loc[:, col].max() * 1.1, n_points)
        y = kde_fit.evaluate(x)
        # Create a set of line segments so that we can color them individually
        # This creates the points as a N x 1 x 2 array so that we can stack points
        # together easily to get the segments. The segments array for line collection
        # needs to be (numlines) x (points per line) x 2 (for x and y)
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        #
        norm = QuantileNormalize(quantiles=intensity_quantiles)
        lc = LineCollection(segments, cmap=cmap, norm=norm, alpha=0.5)
        # Set the values used for colormapping
        lc.set_array(x)
        lc.set_linewidth(1)
        line = ax.add_collection(lc)

    ax.set_ylabel("Density", fontsize=15)
    ax.set_xlabel(intensity_label, fontsize=15)
    ax.tick_params(axis="both", labelsize=15)

    ax.autoscale_view()

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax


@save_plot("n_proteins_vs_quantile")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_n_proteins_vs_quantile_results(
        quantiles: pd.DataFrame, n_proteins: pd.Series, nstd: int = 1, cmap: Union[str, colors.Colormap] = "viridis",
        plot: Optional[Tuple[plt.Figure, plt.Axes]] = None, cbar_ax: Optional[plt.Axes] = None,
        intensity_label: str = "Intensity", fill_between: bool = False, close_plots: str = "all", **kwargs
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    """
    saves plot with prefix: {{name}}
    
    Parameters
    ----------
    quantiles
        quantiles to be plotted
    n_proteins
        number of identified proteins
    nstd
        how many standard deviations should the the line will between
    cmap
        color map to use
    plot
        figure to put plot
    cbar_ax
        axis for colorbar
    intensity_label
        label to put on x label
    fill_between
        should the area around the line be filled
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    kwargs
        {kwargs}

    Returns
    -------

    """
    if plot is not None:
        fig, ax = plot
    else:
        if close_plots is not None:
            plt.close(close_plots)
        fig, ax = plt.subplots(1, 1, figsize=(14, 7))

    if not isinstance(cmap, colors.Colormap):
        cmap = copy(cm.get_cmap(cmap))

    m = n_proteins.sort_values()
    for quant in quantiles.index:
        ax.scatter(quantiles.loc[quant, :], n_proteins, c=[cmap(quant)] * len(n_proteins), alpha=0.5, edgecolors="none")
        popt, pcov = curve_fit(linear, n_proteins, quantiles.loc[quant, :])
        fit = linear(m, *popt)

        ax.plot(fit, m, color=cmap(quant))

        if fill_between:
            perr = np.sqrt(np.diag(pcov))
            popt_up = popt + nstd * perr
            popt_dw = popt - nstd * perr
            fit_up = linear(m, *popt_up)
            fit_dw = linear(m, *popt_dw)
            ax.fill_betweenx(m, fit_dw, fit_up, alpha=.15, color=cmap(quant))

    if cbar_ax is None:
        cbar_ax = fig.add_axes([0.2, 0.00, 0.6, 0.02])  # [left, bottom, width, height]
    cb = ColorbarBase(cbar_ax, cmap=cmap, orientation="horizontal")
    cb.set_label("Quantile", fontsize=15)

    ax.set_ylabel("# detected proteins", fontsize=15)
    ax.set_xlabel(intensity_label, fontsize=15)
    ax.tick_params(axis="both", labelsize=15)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, (ax, cbar_ax)


@save_plot("normalization_overview")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_normalization_overview_results(
        quantiles, n_proteins, intensities, protein_intensities, height: int = 15, intensity_label: str = "Intensity",
        close_plots: str = "all", exp_has_techrep: bool = False, **kwargs
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes, plt.Axes, plt.Axes]]:
    """
    saves plot with prefix: {{name}}

    Parameters
    ----------
    quantiles
        quantiles for save_n_proteins_vs_quantile_results
    n_proteins
        n_proteins for save_n_proteins_vs_quantile_results
    intensities
        intensities for save_kde_results
    protein_intensities
        intensities for boxplot data
    height
        height of the figure
    intensity_label
        name of the experiment
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    if close_plots is not None:
        plt.close(close_plots)
    fig = plt.figure(figsize=(18, 18), constrained_layout=True)
    gs = fig.add_gridspec(height, 2)
    ax_density = fig.add_subplot(gs[0:height // 2, 0])
    ax_nprot = fig.add_subplot(gs[height // 2:height - 1, 0])
    ax_colorbar = fig.add_subplot(gs[height - 1, 0])
    ax_boxplot = fig.add_subplot(gs[0:height, 1])

    fig.suptitle(f"{intensity_label} overview" + (TECHREP_SUFFIX if exp_has_techrep else ""), size=28)
    # order the boxplot data after the number of identified peptides
    boxplot_data = protein_intensities.loc[:, n_proteins.sort_values(ascending=False).index[::-1]]

    plot_kwargs = dict(intensity_label=intensity_label)
    plot_kwargs.update(**kwargs)
    plot_kwargs["save_path"] = None
    plot_kwargs["close_plots"] = None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # ignore warning for mixing constrained layout with tight_layout
        save_kde_results(intensities=intensities,  plot=(fig, ax_density), **plot_kwargs)
        save_n_proteins_vs_quantile_results(quantiles=quantiles, n_proteins=n_proteins, plot=(fig, ax_nprot),
                                            cbar_ax=ax_colorbar, **plot_kwargs)
        save_boxplot_results(boxplot_data, plot=(fig, ax_boxplot), vertical=False, **plot_kwargs)
    ax_density.set_xlim(ax_nprot.get_xlim())

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, (ax_nprot, ax_density, ax_colorbar, ax_boxplot)


@save_plot("intensities_heatmap")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_intensities_heatmap_result(
        intensities: pd.DataFrame, cmap: Union[str, colors.Colormap] = "autumn_r", cmap_bad: str = "dimgray",
        cax: plt.Axes = None, plot: Optional[Tuple[plt.Figure, plt.Axes]] = None, vmax: Optional[float] = None,
        vmin: Optional[float] = None, intensity_label: str = "Intensity", show_suptitle: bool = True,
        close_plots: str = "all", exp_has_techrep: bool = False, **kwargs
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    """
    saves plot with prefix: {{name}}

    Parameters
    ----------
    intensities
        DataFrame containing protein intensities of samples
    cmap
        color map to use for heatmap coloring
    cmap_bad
        color for missing values
    cax
        axis for the color bar
    plot
        figure to put plot
    vmax
        passed to imshow
    vmin
        passed to imshow
    intensity_label
        name of the experiment
    show_suptitle
        should figure title be shown
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}
    """
    if plot is not None:
        fig, ax = plot
    else:
        if close_plots is not None:
            plt.close(close_plots)
        height = max(len(intensities.columns) * 0.5, 4)
        fig, ax = plt.subplots(figsize=(17, height))

    if not isinstance(cmap, colors.Colormap):
        cmap = copy(cm.get_cmap(cmap))
    if cmap_bad is not None:
        cmap.set_bad(color='dimgray')

    im = ax.imshow(intensities.values.T, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax, interpolation="none")
    if cax is None:
        cbar = ax.figure.colorbar(im, ax=ax)
    else:
        cbar = ax.figure.colorbar(im, cax=cax)

    if show_suptitle:
        fig.suptitle(f"Proteins detected or missing in {intensity_label}" + (TECHREP_SUFFIX if exp_has_techrep else ""),
                     fontsize=16)
    ax.set_xlabel("Total proteins detected", fontsize=12)

    y_lim = ax.get_ylim()
    ax.set_yticks(np.linspace(0, len(intensities.columns) - 1, len(intensities.columns)))
    labels = [sample.replace("_", " ") for sample in intensities.columns.values]
    ax.set_yticklabels(labels=labels)
    ax.set_ylim(*y_lim)
    ax.tick_params(axis="both", labelsize=10)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, (ax, cbar.ax)


@save_plot("detection_per_replicate")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_detected_proteins_per_replicate_results(
        all_heights: Dict[str, pd.Series], intensity_label: str = "Intensity", show_suptitle: bool = True,
        close_plots: str = "all", exp_has_techrep: bool = False, **kwargs
):
    """
    saves plot with prefix: {{name}}

    Parameters
    ----------
    all_heights
        mapping of sample to a pd.Series of heights
    intensity_label
        name of the experiment
    show_suptitle
        should the figure title be shown
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    if close_plots is not None:
        plt.close(close_plots)
    # determine number of rows and columns in the plot based on the number of experiments
    n_rows_experiment, n_cols_experiment = get_number_rows_cols_for_fig(all_heights.keys())
    fig, axarr = plt.subplots(n_rows_experiment, n_cols_experiment,
                              figsize=(5 * n_cols_experiment, 3 * n_rows_experiment))
    for i in range(n_rows_experiment * n_cols_experiment - len(all_heights.keys())):
        axarr[n_rows_experiment - 1, n_cols_experiment - 1 - i].remove()

    if show_suptitle:
        fig.suptitle(f"Number of detected proteins from {intensity_label}" + (
                     TECHREP_SUFFIX if exp_has_techrep else ""))

    global_max = max((ser.max() for ser in all_heights.values()))

    for experiment, (pos, ax) in zip(all_heights.keys(), np.ndenumerate(axarr)):
        experiment_heights = all_heights[experiment]
        mean_height = experiment_heights[1:].mean()
        y_pos = [x for x in range(len(experiment_heights))]
        ax.barh(y_pos, experiment_heights, color="skyblue")

        if len(experiment_heights) in range(1, 10):
            fsize = 11
        elif len(experiment_heights) in range(10, 16):
            fsize = 7
        else:
            fsize = 5

        for y, value in zip(y_pos, experiment_heights):
            ax.text(experiment_heights[0] / 2, y, value,
                    verticalalignment='center', horizontalalignment='center', fontsize=fsize)
        labels = experiment_heights.index.values
        labels = [label.replace((experiment + "_"), "").replace("_", " ") for label in labels]

        ax.set_title(experiment.replace("_", "  "))
        ax.axvline(mean_height, linestyle="--", color="black", alpha=0.6)
        ax.set_yticks([i for i in range(len(experiment_heights.index))])
        ax.set_yticklabels(labels, fontsize=fsize)
        ax.set_xlim(0, global_max + global_max*0.02)
        ax.set_xlabel("Counts", fontsize=fsize)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, axarr


@save_plot("intensity_histograms")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_intensity_histogram_results(
        hist_data: pd.DataFrame, intensity_label: str = "Intensity", show_suptitle: bool = False,
        compare_to_remaining: bool = False, legend: bool = True, n_bins: int = 25, show_mean: bool = True,
        histtype="bar", color=None, plot: Optional[Tuple[plt.Figure, plt.Axes]] = None, close_plots: str = "all",
        exp_has_techrep: bool = False, **kwargs
):
    """
    saves plot with prefix: {{name}}

    Parameters
    ----------
    hist_data
        data to be plotted
    intensity_label
        name of experiment
    show_suptitle
        should the figure title be shown
    compare_to_remaining
        should the sample be compared to the overall samples
    legend
        should the legend of the sample names be shown
    show_mean
        should the mean intensity be shown
    n_bins
        how many bins should the histograms have
    histtype
        passed to hist
    color
        passed to hist
    plot
        figure to put plot
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    if plot is not None:
        fig, axarr = plot
    else:
        if close_plots is not None:
            plt.close(close_plots)
        n_rows, n_cols = get_number_rows_cols_for_fig(hist_data.columns.get_level_values(0).unique())
        fig, axarr = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
        for i in range(n_rows * n_cols - len(hist_data.columns.get_level_values(0).unique())):
            axarr[n_rows - 1, n_cols - 1 - i].remove()

    if show_suptitle:
        fig.suptitle(f"{intensity_label} histograms" + (TECHREP_SUFFIX if exp_has_techrep else ""))

    counts = []
    for col in hist_data.columns:
        col_count = pd.cut(hist_data[col], bins=int(np.ceil(n_bins*0.9))).value_counts()
        counts.append(max(col_count))

    for col, (pos, ax) in zip(hist_data.columns.get_level_values(0).unique(), np.ndenumerate(axarr)):
        intensities = hist_data[col]
        try:
            labels = intensities.columns
        except AttributeError:
            labels = col
        labels = [label.replace("_", " ") for label in labels]

        if "Log_2" in intensity_label:
            bins = np.linspace(np.nanmin(intensities.values), np.nanmax(intensities.values), n_bins)
        else:
            bins = np.logspace(np.log2(np.nanmin(intensities.values)), np.log2(np.nanmax(intensities.values)), n_bins,
                               base=2)

        col_name = col.replace("_", "  ")
        ax.set_title(col_name)

        if not legend:
            ax.hist(intensities, bins=bins, histtype=histtype, color=color)
        elif legend and len(labels)<2:
            ax.hist(intensities, bins=bins, histtype=histtype, color=color)
        else:
            ax.hist(intensities, bins=bins, histtype=histtype, color=color, label=labels)

        if compare_to_remaining:
            remaining = hist_data.drop(col, axis=1)
            remaining = remaining.mean(axis=1)
            # remaining = remaining[intensities.notna()]
            ax.hist(remaining, bins=bins, histtype="step", alpha=0.5)
        if "Log_2" not in intensity_label:
            ax.set_xscale("log", base=2)
        ax.set_xlabel(intensity_label)
        ax.set_ylabel("Counts")
        ax.set_xlim(hist_data.min().min(), hist_data.max().max())
        if histtype == ["bar"]:
            ax.set_ylim(0, max(counts))

        if show_mean:
            means = []
            for col in intensities.columns:
                means.append(float(intensities[col].mean()))
            mean = float(sum(means)) / float(len(means))
            ax.axvline(mean, linestyle="--", color="black", alpha=0.6, linewidth=2.5,
                       label='mean: {:5.2f}'.format(mean))
        ax.legend(handlelength=1, handletextpad=0.8, loc='upper right', bbox_to_anchor=(0.98, 0.98), frameon=False)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, axarr


@save_plot("scatter_replicates_{full_name}")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_scatter_replicates_results(
        scatter_data: pd.DataFrame, intensity_label: str = "Intensity", show_suptitle: bool = True,
        close_plots: str = "all", exp_has_techrep: bool = False, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    saves plot with prefix: {{name}}

    Parameters
    ----------
    scatter_data
        data to create scatter plots
    intensity_label
        name of the experiment
    show_suptitle
        should the figure title be shown
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    if close_plots is not None:
        plt.close(close_plots)
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))

    if show_suptitle:
        ax.set_title(f"Scatter comparison of replicates using {intensity_label}" + (
                     TECHREP_SUFFIX if exp_has_techrep else ""))

    min_counts = scatter_data.min().min()
    df_corr = scatter_data.corr()

    for rep1, rep2 in combinations(scatter_data.columns, 2):
        x1 = scatter_data.loc[:, rep1]
        x2 = scatter_data.loc[:, rep2]
        r_2 = df_corr.loc[rep1, rep2] ** 2
        plot_mask = np.logical_or(x1.notna(), x2.notna())
        exp = r"$r^{2}$"
        rep1 = rep1.replace("_", "  ")
        rep2 = rep2.replace("_", "  ")
        ax.scatter(x1.fillna(min_counts)[plot_mask], x2.fillna(min_counts)[plot_mask],
                   label=f"{rep1}  vs  {rep2},  " + fr"{exp}: {r_2:.4f}",
                   alpha=0.5, s=40,  marker=".", edgecolors="none")
        ax.set_xlabel(f"{intensity_label} of x1")
        ax.set_ylabel(f"{intensity_label} of x2")

    fig.legend(frameon=False, bbox_to_anchor=(1.02, 0.5), loc="center left", title=r"$\bf{Sample\ x1\ vs\ Sample\ x2}$")
    if "Log_2" not in intensity_label:
        ax.set_xscale("log")
        ax.set_yscale("log")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax


@save_plot("rank_{full_name}")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_rank_results(
        rank_data: pd.Series, interesting_proteins, intensity_label: str = "Intensity", full_name="Experiment",
        close_plots: str = "all", exp_has_techrep: bool = False, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    saves plot with prefix: {{name}}

    Parameters
    ----------
    rank_data
        data for rank plot
    interesting_proteins
        mapping of pathway names to a list of proteins
    intensity_label
        name of the experiment
    full_name
        name of the sample/group plotted
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    if close_plots is not None:
        plt.close(close_plots)
    if interesting_proteins.values():
        all_pathway_proteins = set.union(*(set(x) for x in interesting_proteins.values()))
    else:
        all_pathway_proteins = set()
    # get the experiment intensities calculate mean intensity for the experiment and sort from highest to lowest
    # TODO apply filter for rare proteins before here?
    # protein ranks vs intensity
    # create dict to map each protein its respective rank and mean intensity
    rank_data = rank_data.dropna()
    dic = {idx: (i, value) for i, (idx, value) in enumerate(rank_data.items())}

    found_proteins = set(rank_data.index)
    # get all proteins that are not part of any pathway
    non_pathway_proteins = found_proteins - all_pathway_proteins
    # get all proteins that are part of any pathway

    # plot the non pathway proteins
    x = [dic[protein][0] for protein in non_pathway_proteins]  # rank
    x_percentage = [xi / len(rank_data) * 100 for xi in x]     # rank as percentage
    y = [dic[protein][1] for protein in non_pathway_proteins]  # intensity

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.scatter(x_percentage, y, c="darkgray", s=30, alpha=0.2, marker=".", label="no pathway", edgecolors="none")

    # plot all proteins of a specific pathway
    legend_text = []
    handles = []
    for i, (pathway, proteins) in enumerate(interesting_proteins.items()):
        proteins = set(proteins) & found_proteins
        if proteins:
            x = [dic[protein][0] for protein in proteins]  # the rank of each protein
            x_percentage_pathway = [xi / len(rank_data) * 100 for xi in x]  # rank as percentage
            y = [dic[protein][1] for protein in proteins]  # the intensity of each protein
            ax.scatter(x_percentage_pathway, y, c=f"C{i}", s=120, alpha=0.7, marker=".", edgecolors="none",
                       label=pathway.replace("_", " "))

            median_pathway_rank = int(np.median(x))
            median_pathway_rank_percentage = median_pathway_rank / len(rank_data) * 100
            median_intensity = rank_data.iloc[median_pathway_rank]
            #xmin, xmax = ax.get_xbound()
            #xm = (median_pathway_rank + abs(xmin)) / (abs(xmax) + abs(xmin))
            xmin, xmax = 0, 100
            xm = (median_pathway_rank_percentage + abs(xmin)) / (abs(xmax) + abs(xmin))
            ymin, ymax = ax.get_ybound()
            ym = (median_intensity - ymin) / (ymax - ymin)
            # plot the median rank and intensity at that rank
            ax.axvline(median_pathway_rank_percentage, ymax=ym, linestyle="--", color=f"C{i}", alpha=0.6)
            ax.axhline(median_intensity, xmax=xm, linestyle="--", color=f"C{i}", alpha=0.6)
            pathway_label = pathway.replace("_", " ")
            text = f"{pathway_label} : median rank: {median_pathway_rank / len(rank_data) * 100 :.1f}% ({len(x)})"
            legend_text.append(text)
            handle = mlines.Line2D([], [], color=f"C{i}", marker='.', markersize=10, label=text, linewidth=0)
            handles.append(handle)

    median_int_total = rank_data.median()
    legend_text.append(f"Median intensity: {median_int_total :.1f} {intensity_label}")
    handle = mlines.Line2D([], [], color="lightgray", marker='.', markersize=10, label=legend_text[-1], linewidth=0)
    handles.append(handle)

    exp_name = full_name.replace("_", " ")
    if "Log_2" not in intensity_label:
        ax.set_yscale("log")
    ax.set_xlabel("Protein rank [%]", size=10, labelpad=10)
    ax.set_ylabel(intensity_label, size=10, labelpad=10)
    fig.suptitle(f"{exp_name} mean" + (TECHREP_SUFFIX if exp_has_techrep else ""), weight="bold", size="14")
    fig.legend(labels=legend_text, handles=handles, bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=False)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    return fig, ax


@save_plot("scatter_comparison_{sample1}_vs_{sample2}")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_experiment_comparison_results(
        protein_intensities_sample1: pd.Series, protein_intensities_sample2: pd.Series,
        exclusive_sample1: pd.Series, exclusive_sample2: pd.Series, sample1: str, sample2: str,
        intensity_label: str = "Intensity", show_suptitle: bool = True,
        plot: Optional[Tuple[plt.Figure, plt.Axes]] = None, close_plots: str = "all", exp_has_techrep: bool = False,
        **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    saves plot with prefix: {{name}}

    Parameters
    ----------
    protein_intensities_sample1
        intensities of sample 1
    protein_intensities_sample2
        intensities of sample 2
    exclusive_sample1
        intensities exclusive to sample 1
    exclusive_sample2
        intensities exclusive to sample 2
    sample1
        name of sample 1
    sample2
        name of sample 2
    intensity_label
        name of experiment
    show_suptitle
        should the figure title be shown
    plot
        figure to put plot
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    # calculate r
    r = protein_intensities_sample1.corr(protein_intensities_sample2)

    if plot is not None:
        fig, ax = plot
    else:
        if close_plots is not None:
            plt.close(close_plots)
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))

    exp = r"$r^{2}$"
    sample1 = sample1.replace("_", " ")
    sample2 = sample2.replace("_", " ")

    ax.scatter(protein_intensities_sample1, protein_intensities_sample2, s=40, alpha=0.7, marker=".", edgecolors="none",
               label=f"{sample1}  vs  {sample2}, {exp}: {r ** 2:.4f}")
    ax.scatter(exclusive_sample1, [np.min(protein_intensities_sample2) * 0.95] * exclusive_sample1.shape[0],
               s=40, alpha=0.7, marker=".", edgecolors="none", label=f"unique for  {sample1}")
    ax.scatter([np.min(protein_intensities_sample1) * 0.95] * exclusive_sample2.shape[0], exclusive_sample2,
               s=40, alpha=0.7, marker=".", edgecolors="none", label=f"unique for  {sample2}")

    ax.set_xlabel(f"{intensity_label} of {sample1}", labelpad=10)
    ax.set_ylabel(f"{intensity_label} of {sample2}", labelpad=10)
    fig.legend(frameon=False, bbox_to_anchor=(1.02, 0.5), loc="center left")

    if "Log_2" not in intensity_label:
        ax.set_xscale("log")
        ax.set_yscale("log")
    if show_suptitle:
        fig.suptitle(f"Scatter comparison of groups using {intensity_label}" + (
                     TECHREP_SUFFIX if exp_has_techrep else ""))

    xmin, xmax = ax.get_xbound()
    ymin, ymax = ax.get_ybound()
    ax.set_xlim(min(xmin, ymin), max(xmax, ymax))
    ax.set_ylim(min(xmin, ymin), max(xmax, ymax))

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax


@save_plot("go_analysis")
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_go_analysis_results(
        heights: Dict[str, list], test_results: Dict[str, list], go_length: Dict[str, list],
        go_analysis_gene_names: list, show_suptitle: bool = True, intensity_label="Intensity", close_plots: str = "all",
        exp_has_techrep: bool = False, **kwargs
) -> Tuple[plt.Figure, plt.Axes]:
    """
    saves plot with prefix: {{name}}

    Parameters
    ----------
    heights
        mapping from samples to bar height
    test_results
        mapping from samples to p value
    go_length
        mapping of name to number of proteins of go terms
    go_analysis_gene_names
        names of the go term lists
    show_suptitle
        should the figure title be shown
    intensity_label
        name of the experiment
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    if close_plots is not None:
        plt.close(close_plots)

    go_analysis_gene_names = [go_term.replace("_", " \n").replace(" + ", " \n") for go_term in go_analysis_gene_names]

    heights_df = pd.DataFrame(data=heights, index=go_analysis_gene_names)
    test_results_df = pd.DataFrame(data=test_results, index=go_analysis_gene_names)
    # save go analysis data
    save_path, csv_name = get_path_and_name_from_kwargs(name="tables/go_analysis_data", **kwargs)
    save_csv_fn(save_path, csv_name, heights_df)
    save_path, csv_name = get_path_and_name_from_kwargs(name="tables/go_analysis_pvals", **kwargs)
    save_csv_fn(save_path, csv_name, test_results_df)
    # add extra col "Total" to data frame - easier data handeling later
    test_results_df.insert(0, "Total", [1] * len(go_analysis_gene_names), True)

    barwidth = 1 / (len(heights_df.columns) + 1)
    max_prot_count = heights_df.max().max()
    tick_y_pos = set()
    fig, ax = plt.subplots(1, 1, figsize=(10, int(len(heights_df.columns) * len(go_analysis_gene_names)/1.5)))
    #for every sample plot a single bar per go term selected in one subplot
    for i, sample in enumerate(heights_df.columns):
        width = 1 / (len(heights_df.columns) + 1) * i
        y_pos = np.arange(len(heights_df)) + width
        tick_y_pos.update(set(y_pos))
        ax.barh(y=y_pos, width=heights_df[sample], height=barwidth, align="center", edgecolor='white')
        #annotate each sample bar with the corresponding p valueis the p value is significant
        for pval, x, y in zip(test_results_df[sample], heights_df[sample], y_pos):
            if pval > 0.05:
                continue
            text = f"{pval:.4f}" if pval > 0.0005 else "< 0.0005"
            ax.annotate(f" p: {text}", xy=(x, y))

    for i, total_count in enumerate(heights_df["Total"]):
        ax.annotate(f"{total_count} / {go_length[i]}", xy=(total_count, list(tick_y_pos)[i]))
    #annotate each set of bars for a go term with the corresponding go term label = second layer of y axis labeling
    for i, go in enumerate(go_analysis_gene_names):
        text = go.split(".txt")[0].replace("_", " ")
        ax.annotate(text, xy=(-max_prot_count * 0.8, i + 0.5 - (1 / (2 * (len(heights_df.columns) + 1)))),
                    verticalalignment="center", annotation_clip=False)

    if show_suptitle:
        fig.suptitle(f"GO based analysis from {intensity_label}" + (TECHREP_SUFFIX if exp_has_techrep else ""))
    ax.set_yticks(sorted(list(tick_y_pos)))
    labels = [label.replace("_", " ") for label in heights_df.columns]
    ax.set_yticklabels(labels * len(heights_df))
    ax.set_xlim(0, max_prot_count + max_prot_count * 0.2)
    ax.set_xlabel('Number of detected proteins from GO list')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    return fig, ax


@save_plot("pathway_timecourse_{pathway}")
def save_pathway_timecourse_results():
    """
    Not Implemented at the moment
    """
    raise NotImplementedError
    """plt.close("all")
    n_rows, n_cols = get_number_rows_cols_for_fig(found_proteins)
    fig, axarr = plt.subplots(n_rows, n_cols, figsize=(n_cols * int(max_time / 5), 4 * n_rows))
    if show_suptitle:
        fig.suptitle(pathway)

    protein_minimum = self.all_intensities_dict[df_to_use].max().max()
    protein_maximum = self.all_intensities_dict[df_to_use].min().min()
    for protein, (pos, ax) in zip(found_proteins, np.ndenumerate(axarr)):
        ax.set_title(protein)
        ax.set_xlabel(f"Age [weeks]")
        ax.set_ylabel(intensity_label)
        for idx, experiment in enumerate(level_keys):
            protein_intensities = self.all_tree_dict[df_to_use][experiment].aggregate(None, index=protein)
            mask = protein_intensities > 0
            protein_minimum = min(protein_minimum, protein_intensities[mask].min())
            protein_maximum = max(protein_maximum, protein_intensities[mask].max())
            ax.scatter([x_values[experiment]] * sum(mask), protein_intensities[mask],
                       label=f"{groups[experiment]}", color=group_colors[groups[experiment]])
    # adjust labels based on overall min and max of the pathway
    for protein, (pos, ax) in zip(found_proteins, np.ndenumerate(axarr)):
        ax.set_ylim(bottom=protein_minimum * 0.99, top=protein_maximum * 1.01)
        ax.set_xlim(left=0, right=max_time + 1)
    handles, labels = next(np.ndenumerate(axarr))[1].get_legend_handles_labels()
    fig.legend(handles, labels, bbox_to_anchor=(1.04, 0.5), loc="center left")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    return fig, axarr"""


@save_plot("plots/venn_bar_{ex}")
@save_venn_to_txt({"named_sets": "txts/set_bar"})
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_bar_venn(
        named_sets: Dict[str, set], ex: str, show_suptitle: bool = True, close_plots: str = "all",
        exp_has_techrep: bool = False, **kwargs
) -> Optional[Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]]:
    """
    saves plot with prefix: {{name}}

    Parameters
    ----------
    named_sets
        a mapping of samples to protein names
    ex
        figure title
    show_suptitle
        should the figure title be shown
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    if close_plots is not None:
        plt.close(close_plots)
    if len(named_sets) > 6:
        warnings.warn(f"Skipping bar-venn for {ex} because it has more than 6 experiments")
        return

    # create a mapping from name to a y coordinate
    y_mappings = {name: i for i, name in enumerate(named_sets)}
    # get all the heights and other info required for the plot
    heights = []
    x = []
    ys = []
    for i, (intersected, unioned, result) in enumerate(venn_names(named_sets)):
        heights.append(len(result))
        x.append(i)
        ys.append([y_mappings[x] for x in intersected])

    # initial figure setup
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(1 * len(heights), 7))
    if show_suptitle:
        fig.suptitle(ex.replace("_", " ") + (TECHREP_SUFFIX if exp_has_techrep else ""), fontsize=17, weight="bold")
    # create the bar plot
    ax1.bar(x, heights, color="skyblue")
    # add text to the bar plot
    for x_level, height in zip(x, heights):
        ax1.text(x_level, max(heights) / 2, height, verticalalignment='center', horizontalalignment='center')
    ax1.set_ylabel("Number of proteins")

    labels = [sample.split("_", 1)[1] if "_" in sample else sample for sample in y_mappings]
    labels = [sample.replace("_", " ") for sample in labels]

    # create the line plots
    for x_level, y in zip(x, ys):
        # we just want to draw a straight line every time so we repeat x as often as needed
        ax2.plot([x_level] * len(y), y, linestyle="-", color="gray", marker=".")
    # replace the yticks with the names of the samples
    ax2.set_yticks([i for i in range(len(labels))])
    ax2.set_yticklabels(labels)
    ax2.set_ylabel("Sample name")
    ax2.set_xlabel("Number of comparison")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, (ax1, ax2)


@save_plot("plots/venn_replicate_{ex}")
@save_venn_to_txt({"named_sets": "txts/set"})
@format_docstrings(kwargs=_get_path_and_name_kwargs_doc)
def save_venn(
        named_sets: Dict[str, set], ex: str, show_suptitle: bool = True, title_font_size=20, set_label_font_size=16,
        subset_label_font_size=14, close_plots: str = "all", exp_has_techrep: bool = False, **kwargs
) -> Optional[Tuple[plt.Figure, plt.Axes]]:
    """
    Creates Venn Diagrams from passed data. saves plot with prefix: {{name}}

    Parameters
    ----------
    named_sets
        a mapping of samples to protein names
    ex
        title for the plot
    show_suptitle
        should the figure title be shown
    title_font_size
        font size of the title
    set_label_font_size
        font size of sets
    subset_label_font_size
        font size of subsets
    close_plots
        which plots should be closed when creating the plot, if None no plots will be closed
    exp_has_techrep
        whether technical replicates were aggregated for the plot
    kwargs
        {kwargs}

    """
    if close_plots is not None:
        plt.close(close_plots)
    fig, ax = plt.subplots(1, 1, figsize=(14, 7))
    if show_suptitle:
        fig.suptitle(ex.replace("_", " ") + (TECHREP_SUFFIX if exp_has_techrep else ""),
                     fontsize=title_font_size, weight="bold")

    # create venn diagram based on size of set
    sets = named_sets.values()
    set_names = named_sets.keys()
    split_names = False
    for name in set_names:
        if "_" in name:
            split_names = True
    if split_names:
        set_names = [sample.split("_", 1)[1] for sample in set_names]
    set_names = [name.replace("_", " ") for name in set_names]

    if len(sets) < 2:
        warnings.warn(f"Could not create venn diagram for {ex} because it has less than 2 replicates")
        return
    elif len(sets) == 2:
        venn = venn2(subsets=sets, set_labels=set_names, ax=ax)
    elif len(sets) == 3:
        venn = venn3(subsets=sets, set_labels=set_names, ax=ax)
    else:
        warnings.warn(f"Could not create venn diagram for {ex}"
                      f" because it has more than 3 replicates ({len(sets)})")
        return

    # if a figure was created, do some further configuration
    for text in venn.set_labels:
        try:
            text.set_fontsize(set_label_font_size)
        except AttributeError:
            pass
    handles = []
    labels = []
    for text, patch in zip(venn.subset_labels, venn.patches):
        try:
            handles.append(patch)
            labels.append(text.get_text())
            text.set_fontsize(subset_label_font_size)
        except AttributeError:
            pass
    fig.legend(handles, labels, bbox_to_anchor=(1.02, 0.5), loc="center left")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, ax
