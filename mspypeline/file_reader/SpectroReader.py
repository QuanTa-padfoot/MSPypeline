import os
import pandas as pd
import logging
from collections import defaultdict
from mspypeline.core.MSPPlots import SpectroPlotter
from mspypeline.helpers import dict_depth
from mspypeline.file_reader import BaseReader, MissingFilesException

class SpectroReader(BaseReader):
    """
    | A child class of the :class:`~BaseReader`.
    | The SpectroReader preprocesses data from Spectronaut file into the internal data format to provide the correct input
      for the plotters. Required files to start the SpectroReader is the xls files from Spectronaut.
    | use_imputed is a class variable that controls whether to use the imputed values from Spectronaut.
    """
    name = "spectroReader"
    required_files = ['.xls, .tsv, or .csv file']
    plotter = SpectroPlotter
    use_imputed = False

    # separators tried (in order) when the input file's delimiter is not known upfront
    _candidate_separators = [",", "\t", ";"]
    _required_columns = ('PG.ProteinGroups', 'PG.Genes', 'PG.ProteinDescriptions')

    def __init__(self, start_dir: str,
                reader_config: dict,
                index_col: str = "PG.Genes",
                loglevel: int = logging.DEBUG):
        """
        Parameters
        ----------
        start_dir
            location where the directory/txt folder to the data can be found.
        reader_config
            mapping of the file reader configuration (as e.g. given in the config.yml file)
        index_col
            with which identification type should detected proteins in the *proteins.txt* file be handled.
            If provided in the reader_config will be taken from there.
        loglevel
            level of the logger
        """
        super().__init__(start_dir, reader_config, loglevel=loglevel)
        self.data_dir = self.start_dir
        self.index_col = index_col

        try:
            file_dir = self.ext_change(self.data_dir)[0]
        except IndexError:
            raise MissingFilesException("Could not find all ")

        # cache the resolved file path + detected separator so preprocess_proteins (called lazily
        # when full_data["proteins"] is first accessed) can reuse them instead of re-scanning the
        # directory and re-sniffing the separator from scratch
        self.file_dir = file_dir
        self.separator = None
        formatted_proteins_txt_columns = None
        for cur_separator in self._candidate_separators:
            try:
                df = pd.read_csv(file_dir, sep=cur_separator, nrows=2)
            except (pd.errors.ParserError, UnicodeDecodeError, ValueError) as e:
                self.logger.warning("SpectroReader cannot open file with (%s) separator: %s", cur_separator, e)
                continue
            if all(col in df.columns for col in self._required_columns):
                self.separator = cur_separator
                df = self.preprocess_df(df)
                header = [element for element in df.columns if ".Quantity" in element]
                formatted_proteins_txt_columns, self.analysis_design = self.format_spektrocols(header)
                self.intensity_column_names = formatted_proteins_txt_columns
                break
        if self.separator is None:
            raise MissingFilesException("Could not find all ")

        if not self.reader_config.get("all_replicates", False):
            self.reader_config["all_replicates"] = formatted_proteins_txt_columns
        if not self.reader_config.get("analysis_design", False):
            self.reader_config["analysis_design"] = self.analysis_design
            self.reader_config["levels"] = dict_depth(self.analysis_design)
            self.reader_config["level_names"] = [x for x in range(self.reader_config["levels"])]

    def preprocess_proteins(self):
        """
        Indices are set to the value given in initialization (PG.Genes per default). Quantity columns are extracted.
        If selected imputed values are set to 0 and duplicate indices are renamed.
        """
        df = pd.read_csv(self.file_dir, sep=self.separator)
        df = self.preprocess_df(df)
        use_index = self.format_double_indx(df[self.index_col])
        df.set_index(use_index, drop=False, inplace=True)
        missing_map = df.filter(regex=(".IsIdentified")).replace({"Filtered": False, "True": True, "False": False})
        missing_map.fillna(False, inplace = True)
        df_result = None
        for intensity in ['.Quantity', '.LFQ', '.IBAQ']:
            quant_cols = [col for col in df.columns if intensity in col]
            if quant_cols != []:
                df1 = df.filter(regex=(intensity))
                df1 = df1.replace({"Filtered": float(0)}).fillna(0)
                use_cols = [col.split(' ')[1] for col in df1.columns]
                if intensity == '.Quantity':
                    use_cols = ["Intensity " + col for col in use_cols]
                elif intensity == '.IBAQ':
                    use_cols = ["iBAQ " + col for col in use_cols]
                elif intensity == '.LFQ':
                    use_cols = ["LFQ intensity " + col for col in use_cols]
                use_cols = [col.split('.PG')[0] for col in use_cols]
                use_cols = [col.split('.raw')[0] for col in use_cols]
                df1.columns = use_cols
                if self.use_imputed is False:
                    df1 = pd.DataFrame(df1.values * missing_map.values, columns=df1.columns, index=df1.index)
                if df_result is None:
                    df_result = df1
                else:
                    df_result = df_result.join(df1)
        df_result.index = df.index.fillna("nan")
        return df_result

    def ext_change(self, data_dir):
        list_files=[]
        for file in os.listdir(data_dir):
            if file.endswith(".xls"):
                filename, ext = os.path.splitext(file)
                new_filename = filename + '.csv'
                old_filedir = os.path.join(data_dir, file)
                new_filedir = os.path.join(data_dir, new_filename)
                os.replace(old_filedir, new_filedir)
                list_files.append(new_filedir)
            elif file.endswith(".csv") or file.endswith(".tsv"):
                list_files.append(os.path.join(data_dir, file))
        return(list_files)

    def format_double_indx(self, indx):
        """
        Used to rename duplicate indices like so: dupl_indx, dupl_indx => dupl_indx_1, dupl_indx_2

        Returns
        -------
        Series
            Pandas Series containing the renamed indices.
        """
        s = pd.Series([str(ind) for ind in indx])
        is_dup = s.duplicated(keep=False)
        occurrence = (s.groupby(s).cumcount() + 1).astype(str)
        return s.where(~is_dup, s + "_" + occurrence)

    def format_spektrocols(self, cols):
        """
        Reformats the column naming from Spektronaut to make it compatible with the setup of msypypeline.
        Columns are split at `.raw` (old Spectronaut version) or `.PG` (newer version) to remove the .Quantity tag. They are then split on the underscore and
        are reordered and put together, so they have the same number of "components". This is necessary due to
        the way analysis_design works in mspypeline.
        """
        processed_cols = []
        analysis_design = {}
        for col in cols:
            cur_col = col.split("] ")[1]
            cur_col = cur_col.split(".PG")[0]
            cur_col = cur_col.split(".raw")[0]
            processed_cols.append(cur_col)
            self.dictizeString(cur_col, cur_col, analysis_design)
        return processed_cols, analysis_design


    def dictizeString(self, string, final_value, dictionary):
        """
        Helper function to generate the analysis design dictionary.
        """
        parts = string.split('_', 1)
        if len(parts) > 1:
            branch = dictionary.setdefault(parts[0], {})
            self.dictizeString(parts[1], final_value, branch)
        else:
            if parts[0] in dictionary:
                # This part is True when there are duplicate sample names. Best to report this as an error
                raise ValueError(f"Found duplicate sample name for [{final_value}], please rename one(s) of these samples differently.")
            else:
                dictionary[parts[0]] = final_value

    def preprocess_df(self, df):
        '''
        Initial processing of the dataframe if numbers were imported as string by read_csv. Also process rows with
        multiple protein names.

        Returns df after processing
        '''
        # process rows with multiple protein names
        dupl_row = df.index[df["PG.ProteinGroups"].str.contains(";", na=False)]
        cols_to_clean = [c for c in df.columns
                 if c not in ["PG.Genes", "PG.ProteinGroups", "PG.ProteinDescriptions"]]

        def take_first_if_split(x):
            if isinstance(x, str) and ";" in x:
                return x.split(";", 1)[0]  # first value only
            return x  # leave as-is (handles no ";" or non-strings)

        df.loc[dupl_row, cols_to_clean] = df.loc[dupl_row, cols_to_clean].applymap(take_first_if_split)

        df["PG.Genes"] = df["PG.Genes"].fillna(df["PG.ProteinGroups"])
        df.index = range(len(df.index))

        # convert non-numeric intensities to numeric:
        value_col = [col for col in df.columns if '.Quantity' in col or '.IBAQ' in col or '.iBAQ' in col]
        df[value_col] = df[value_col].replace(',','.', regex=True)
        df[value_col] = df[value_col].apply(pd.to_numeric, errors = "coerce")

        # rename columns according to sample_mapping.txt
        try:
            sample_mapping = os.path.join(self.start_dir, "config", "sample_mapping.txt")
            mapping = pd.read_csv(sample_mapping, sep="\t", header=0, dtype=str)
            for old_name, new_name in zip(mapping.iloc[:, 0], mapping.iloc[:, 1]):
                old_col = df.filter(regex=old_name).columns.to_list()
                rename_col_dict = {col: col.replace(old_name, new_name) for col in old_col}
                df.rename(columns=rename_col_dict, inplace=True)
        except FileNotFoundError:
            self.logger.debug("File sample_mapping.txt not found in the config folder")
        return df
