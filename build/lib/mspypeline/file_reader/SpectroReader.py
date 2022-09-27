from msilib.schema import Error
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
      for the plotters. Required files to start the SpectroReader is the **proteins.txt** file from Spectronaut.
    | use_imputed and format_doubleindx are class variables that control whether to use the imputed values from Spectronaut
      and whether to rename duplicate row indices.
    """
    name = "spectroReader"
    protein_txt = "proteins.txt"
    required_files = [protein_txt]
    plotter = SpectroPlotter
    use_imputed = False
    format_doubleindx = False

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
            file_dir = os.path.join(self.data_dir, "proteins.txt")
            separators = [",", "\t", ";"]
            df = pd.DataFrame()
            while(df.shape[1] <= 1) & bool(separators):
                try:
                    df = pd.read_csv(file_dir, sep=separators.pop(0))
                except:
                    print("Unable to open file with given seperator")
            df = df.filter(regex=(".Quantity"))
            formatted_proteins_txt_columns, self.analysis_design = self.format_spektrocols(df.columns)
            self.intensity_column_names = formatted_proteins_txt_columns
        except FileNotFoundError:
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
        file_dir = os.path.join(self.data_dir, SpectroReader.protein_txt)
        separators = [",", "\t", ";"]
        df = pd.DataFrame()
        while(df.shape[1] <= 1) & bool(separators):
            try:
                df = pd.read_csv(file_dir, sep=separators.pop(0))
            except:
                print("Unable to open file with" + (separators[0]) + "seperator")
        use_index = df[self.index_col]
        missing_map  = df.filter(regex=(".IsIdentified"))
        df = df.filter(regex=(".Quantity"))
        use_cols = ["Intensity " + col for col in self.intensity_column_names]
        df.columns = use_cols
        if self.use_imputed is False:
            df = pd.DataFrame(df.values * missing_map.values, columns=df.columns, index=df.index)
        if self.format_double_indx:
            use_index = self.format_double_indx(use_index)
        df.set_index(use_index, drop=False, inplace=True)
        df.index = df.index.fillna("nan")
        return df

    def format_double_indx(self, indx):
        """
        Used to rename duplicate indices like so: dupl_indx, dupl_indx => dupl_indx_1, dupl_indx_2
        
        Returns
        -------
        Series
            Pandas Series containing the renamed indices.
        """
        doubled_indx = set([str(ind) for ind in indx[indx.duplicated()]])
        new_indx = []
        cnt = defaultdict(lambda: 1)
        for _ind in list(indx):
            ind = str(_ind)
            if ind in doubled_indx:
                    new_indx.append(ind + "_" + str(cnt[ind]))
                    cnt[ind] += 1            
            else:
                new_indx.append(ind)
        return pd.Series(new_indx)

    def format_spektrocols(self, cols):
        """
        Reformats the column naming from Spektronaut to make it compatible with the setup of msypypeline.
        Columns are split at the dot to remove the .Quantity tag. They are then split on the underscore and
        are reordered and put together, so they have the same number of "components". This is necessary due to 
        the way analysis_design works in mspypeline.
        """
        processed_cols = []
        analysis_design = {}
        for col in cols:
            cur_col = col.split(".")[0]
            cur_col = cur_col.split("_")[2:]
            if len(cur_col)==3:
                new_col = "_".join(cur_col)
            elif len(cur_col)==5:
                new_col = cur_col[0] + cur_col[1] + "_" + cur_col[2] + "_" + cur_col[3] + cur_col[4]
            elif len(cur_col)==4:
                if("rep" in cur_col[3]):
                    new_col = cur_col[0] + "_" + cur_col[1] + "_" + cur_col[2] + cur_col[3]
                else:
                    new_col = cur_col[0] + cur_col[1] + "_" + cur_col[2] + "_" + cur_col[3]
            processed_cols.append(new_col)
            self.dictizeString(new_col, new_col, analysis_design)
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
                dictionary[parts[0]] += 1
            else:
                dictionary[parts[0]] = final_value
