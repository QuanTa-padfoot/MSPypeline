import logging
from mspypeline.core import MSPInitializer
from mspypeline.core.MSPPlots import BasePlotter

class SpectroPlotter(BasePlotter):
    """
    SpectroPlotter is a child class of the :class:`BasePlotter` and inherits all functionality to get data and
    generate plots.
    """
    def __init__(
        self,
        start_dir: str,
        reader_data: dict,
        intensity_df_name: str = "proteins",
        interesting_proteins: dict = None,
        go_analysis_gene_names: dict = None,
        configs: dict = None,
        required_reader="spectroReader",
        intensity_entries=(("raw", "Intensity ", "Intensity"), ("lfq", "LFQ intensity ", "LFQ intensity"), ("ibaq", "iBAQ ", "iBAQ intensity")),
        loglevel=logging.DEBUG
    ):
        """
        Parameters
        ----------
        start_dir
            location to save results
        reader_data
            mapping to provide input data
        intensity_df_name
            name/key to input data
        interesting_proteins
            mapping with pathway proteins to analyze
        go_analysis_gene_names
            mapping with go terms to analyze
        configs
            mapping of configuration
        required_reader
            name of the file reader
        intensity_entries
            tuple of (key in all_tree_dict, prefix in data, name in plot). See :meth:`add_intensity_column`.
        loglevel
            level of the logger
        """
        
        super().__init__(
            start_dir,
            reader_data,
            intensity_df_name,
            interesting_proteins,
            go_analysis_gene_names,
            configs,
            required_reader,
            intensity_entries,
            loglevel
        )

    @classmethod
    def from_MSPInitializer(cls, mspinit_instance: MSPInitializer, **kwargs):
        default_kwargs = dict(
            intensity_entries = (("raw", "Intensity ", "Intensity"), ("lfq", "LFQ intensity ", "LFQ intensity"),
                                ("ibaq", "iBAQ", "iBAQ intensity")),
            intensity_df_name="proteins",
            required_reader="spectroReader"
        )
        default_kwargs.update(**kwargs)
        return super().from_MSPInitializer(mspinit_instance, **default_kwargs)

    @classmethod
    def from_file_reader(cls, reader_instance, **kwargs):
        default_kwargs = dict(
            intensity_df_name="proteins",
            intensity_entries=(("raw", "Intensity ", "Intensity"), ("lfq", "LFQ intensity ", "LFQ intensity"),
                               ("ibaq", "iBAQ ", "iBAQ intensity")),
        )
        default_kwargs.update(**kwargs)
        return super().from_file_reader(reader_instance, **kwargs)
        