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
                                ("ibaq", "iBAQ ", "iBAQ intensity")),
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
        
    def read_peptide_data(self, dir_peptide_data_folder = None, time_level = None):
        """Read and preprocess the .csv file containing peptide data from Spectronaut. The data is first log2-transformed by default.
        Rows with multiple protein names are split to separate rows, one for each protein. The time points are detected from the sample names.
        By default (i.e., if time_level is None), the second-to-last level in analysis_design is used for the detection.

        Parameters
        ----------
        dir_peptide_data_folder
            directory of the folder containing the peptide data, which must be in .csv or .tsv format. If None, MSPypeline will search in a folder named `peptide_data` in the current directory
            The folder `peptide_data` should contain only 1 .csv or .tsv file to avoid reading the wrong file
        time_level
            Which level in analysis_design the time points or doses were declared. If None, the second-to-last level will be used

        Returns
        --------
        peptide_df
            Dataframe containing the expression value of the peptides.
        timepoints_dict
            Dictionary specifying the time points or dose for each sample.
        time_numerical_dict
            Dictionary specifying the numerical values of the time points or dose for each sample
        """
        import numpy as np
        import os
        import pandas as pd
        import re
        peptide_dir=[]
        for file in os.listdir(dir_peptide_data_folder):
            if file.endswith(".xls"):
                filename, ext = os.path.splitext(file)
                new_filename = filename + '.csv'
                old_filedir = os.path.join(dir_peptide_data_folder, file)
                new_filedir = os.path.join(dir_peptide_data_folder, new_filename)
                os.replace(old_filedir, new_filedir)
                peptide_dir.append(new_filedir)
            elif file.endswith(".csv") or file.endswith(".tsv"):
                peptide_dir.append(os.path.join(dir_peptide_data_folder, file))
        try:
            peptide_dir = peptide_dir[0]
        except (IndexError):
            peptide_dir = []
        for cur_separator in ['\t', ';', ',']:
            try:
                peptide_df = pd.read_csv(peptide_dir, delimiter=cur_separator, dtype=str)
                if ('PG.ProteinGroups' in peptide_df.columns) and ('PG.Genes' in peptide_df.columns) and ('PEP.StrippedSequence' in peptide_df.columns) and ('EG.PrecursorId' in peptide_df.columns):
                    break
            except (ValueError):
                pass
        peptide_df.columns = [col.replace(" ","") for col in peptide_df.columns]
        # remove "PG.IsSingleHit", ".EG.Quantity", and ".Qvalue" columns
        peptide_df.drop(peptide_df.filter(regex='PG.IsSingleHit|PG.Quantity|EG.Qvalue').columns, axis=1, inplace=True)
        # rename columns
        all_quant_cols = peptide_df.filter(regex="EG.TotalQuantity").columns.to_list()
        all_prefixes = [s.split(".EG")[0] for s in all_quant_cols]
        all_prefixes = [s.split(".raw")[0] for s in all_prefixes]
        all_sample_name = [s.split("]")[1] for s in all_prefixes]
        rename_dict = {all_quant_cols[i]: all_sample_name[i] for i in range(len(all_quant_cols))} 
        peptide_df.rename(columns=rename_dict, inplace= True)
        
        sample_mapping = os.path.join(os.path.dirname(dir_peptide_data_folder), "config/sample_mapping.txt")
        try:
            with open(sample_mapping, "r") as f:
                next(f)  # skip the title line
                for line in f.readlines():
                    sample_name = line.split('\t')
                    old_col = peptide_df.filter(regex=sample_name[0]).columns.to_list()
                    rename_col_dict = {col: col.replace(sample_name[0], sample_name[1][:-1]) for col in old_col}
                    peptide_df.rename(columns=rename_col_dict, inplace=True)
                f.close()
        except FileNotFoundError:
            pass

        timepoints_dict, time_numerical_dict = {}, {}
        return {"peptide_df": peptide_df,"timepoints_dict": timepoints_dict,"time_numerical_dict": time_numerical_dict}
