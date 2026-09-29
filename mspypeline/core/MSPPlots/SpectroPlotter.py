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
        if not peptide_dir:
            raise FileNotFoundError(f"No .csv, .tsv or .xls peptide data found in {dir_peptide_data_folder}")
        peptide_dir = peptide_dir[0]
        required_cols = ('PG.ProteinGroups', 'PG.Genes', 'PEP.StrippedSequence', 'EG.PrecursorId')
        # only sniff the separator on the first rows, a peptide report of a larger experiment is far too big to be
        # read once per candidate separator
        separator, raw_columns = None, None
        for cur_separator in ['\t', ';', ',']:
            try:
                candidate = pd.read_csv(peptide_dir, delimiter=cur_separator, dtype=str, nrows=2)
            except (ValueError, pd.errors.ParserError):
                continue
            if all(col in [c.replace(" ", "") for c in candidate.columns] for col in required_cols):
                separator, raw_columns = cur_separator, list(candidate.columns)
                break
        if separator is None:
            # without a clear message this used to surface as an UnboundLocalError further down
            raise ValueError(f"Could not read {peptide_dir} as a Spectronaut peptide report: none of the tested "
                             f"separators (tab, ';', ',') produced the required columns {required_cols}")
        # "PG.IsSingleHit", "PG.Quantity" and "EG.Qvalue" are not used, so they are never read in the first place.
        # A report of a real experiment has one of those per run, i.e. they make up the bulk of the file.
        drop_pattern = re.compile('PG.IsSingleHit|PG.Quantity|EG.Qvalue')
        use_columns = [col for col in raw_columns if not drop_pattern.search(col.replace(" ", ""))]
        quant_pattern = re.compile(r"EG\.TotalQuantity")
        # read in chunks and narrow the quantity columns to float right away, keeping 200k+ rows of a peptide
        # report as python strings needs several GB of memory
        chunks = []
        for chunk in pd.read_csv(peptide_dir, delimiter=separator, dtype=str, usecols=use_columns,
                                 chunksize=50000):
            chunk.columns = [col.replace(" ", "") for col in chunk.columns]
            quant_cols = [col for col in chunk.columns if quant_pattern.search(col)]
            # Spectronaut writes "Filtered" for runs without a quantity and may use ',' as the decimal separator
            chunk[quant_cols] = chunk[quant_cols].replace(',', '.', regex=True).apply(pd.to_numeric,
                                                                                     errors="coerce").astype("float32")
            chunks.append(chunk)
        if not chunks:
            raise ValueError(f"{peptide_dir} does not contain any peptide data")
        peptide_df = pd.concat(chunks, ignore_index=True) if len(chunks) > 1 else chunks[0]
        del chunks
        # rename columns
        all_quant_cols = peptide_df.filter(regex=r"EG\.TotalQuantity").columns.to_list()
        all_prefixes = [s.split(".EG")[0] for s in all_quant_cols]
        all_prefixes = [s.split(".raw")[0] for s in all_prefixes]
        # Spectronaut prefixes the run with "[<n>]"; not every export has it, so only strip it when it is there
        all_sample_name = [s.split("]")[1] if "]" in s else s for s in all_prefixes]
        rename_dict = {all_quant_cols[i]: all_sample_name[i] for i in range(len(all_quant_cols))}
        peptide_df.rename(columns=rename_dict, inplace= True)
        
        sample_mapping = os.path.join(os.path.dirname(dir_peptide_data_folder), "config/sample_mapping.txt")
        try:
            with open(sample_mapping, "r") as f:
                next(f)  # skip the title line
                for line in f.readlines():
                    sample_name = line.split('\t')
                    if len(sample_name) < 2:
                        continue
                    old_name, new_name = sample_name[0].strip(), sample_name[1].strip()
                    # the raw file names may contain regex metacharacters (e.g. "."), so match them literally
                    old_col = peptide_df.filter(regex=re.escape(old_name)).columns.to_list()
                    rename_col_dict = {col: col.replace(old_name, new_name) for col in old_col}
                    peptide_df.rename(columns=rename_col_dict, inplace=True)
                f.close()
        except FileNotFoundError:
            pass

        timepoints_dict, time_numerical_dict = {}, {}
        return {"peptide_df": peptide_df,"timepoints_dict": timepoints_dict,"time_numerical_dict": time_numerical_dict}
