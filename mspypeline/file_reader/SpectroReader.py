from fileinput import filename
from msilib.schema import Error
import os
import re
import pandas as pd
import logging
from collections import defaultdict
from mspypeline.core.MSPPlots import SpectroPlotter
from mspypeline.helpers import dict_depth
from mspypeline.file_reader import BaseReader, MissingFilesException
from itertools import groupby

class SpectroReader(BaseReader):
    """
    | A child class of the :class:`~BaseReader`.
    | The SpectroReader preprocesses data from Spectronaut file into the internal data format to provide the correct input
      for the plotters. Required files to start the SpectroReader is the xls files from Spectronaut.
    | use_imputed and format_doubleindx are class variables that control whether to use the imputed values from Spectronaut
      and whether to rename duplicate row indices.
    """
    name = "spectroReader"
    required_files = ['.xls or .csv file']
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
            try:
                file_dir = self.ext_change(self.data_dir)[0]
            except (IndexError):
                file_dir = []
            df = pd.DataFrame()
            with open(file_dir, "r") as df_head:
                head = [next(df_head) for x in range(4)]
            breaker = False
            for cur_separator in ["\t", ";", ","]:
                h = head[0].split(cur_separator)
                if (len(h)>3):
                    for line in head[1:4]:
                        l = line.split(cur_separator)
                        l = l[3:]
                        if bool(any("," in element for element in l)):
                            try:
                                df = pd.read_csv(file_dir, sep=cur_separator, decimal=",")
                                breaker = True
                                break
                            except:
                                print(f"Unable to open file with ({cur_separator}) separator and (,) decimal point")
                    if not breaker:
                        try:
                            df = pd.read_csv(file_dir, sep=cur_separator, decimal=".")
                        except:
                            print(f"Unable to open file with ({cur_separator}) separator and (.) decimal point")
            if (df.shape[1] > 3):
                for column in ["PG.ProteinGroups", "PG.Genes", "PG.ProteinDescriptions"]:
                    for row in df.index:
                        df[column][row] = df[column][row].replace("!", ";")
                        df[column][row] = df[column][row].replace("@", ",")
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
        try:
            file_dir = self.ext_change(self.data_dir)[0]
        except IndexError:
            file_dir = []
        separators = [",", "\t", ";"]
        df = pd.DataFrame()
        with open(file_dir, "r") as df_head:
            head = [next(df_head) for x in range(4)]
        breaker = False
        for cur_separator in ["\t", ";", ","]:
            h = head[0].split(cur_separator)
            if (len(h) > 3):
                for line in head[1:4]:
                    l = line.split(cur_separator)
                    l = l[3:]
                    if bool(any("," in element for element in l)):
                        try:
                            df = pd.read_csv(file_dir, sep=cur_separator, decimal=",")
                            breaker = True
                            break
                        except:
                            print(f"Unable to open file with ({cur_separator}) separator and (,) decimal point")
                if not breaker:
                    try:
                        df = pd.read_csv(file_dir, sep=cur_separator, decimal=".")
                    except:
                        print(f"Unable to open file with ({cur_separator}) separator and (.) decimal point")
        if (df.shape[1] > 3):
            for column in ["PG.ProteinGroups", "PG.Genes", "PG.ProteinDescriptions"]:
                for row in df.index:
                    df.loc[row, column] = df.loc[row, column].replace("!", ";")
                    df.loc[row, column] = df.loc[row, column].replace("@", ",")
        use_index = df[self.index_col]
        quant_cols = [col for col in df.columns if '.Quantity' in col]
        filt_cols = [col for col in df.columns if '.IsIdentified' in col]
        quantfilt_cols = quant_cols + filt_cols
        keyf = lambda text: text.split(".")[0]
        grouped_cols = [list(items) for gr, items in groupby(sorted(quantfilt_cols), key=keyf)]
        for i in range(0, len(grouped_cols)):
            df.loc[pd.isna(df[grouped_cols[i][1]]) == True, grouped_cols[i][0]]=False

        missing_map  = df.filter(regex=(".IsIdentified")).replace({"Filtered": False, "True": True, "False": False})
        df = df.filter(regex=(".Quantity"))
        df = df.replace({"Filtered": float(0)}).fillna(0)
        
        use_cols = ["Intensity " + col for col in self.intensity_column_names]
        df.columns = use_cols
        if self.use_imputed is False:
            df = pd.DataFrame(df.values * missing_map.values, columns=df.columns, index=df.index)
            
        if self.format_double_indx:
            use_index = self.format_double_indx(use_index)
        df.set_index(use_index, drop=False, inplace=True)
        df.index = df.index.fillna("nan")
        return df

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
            elif file.endswith(".csv"):
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
            cur_col = col.split(" ")[1]
            cur_col = cur_col.split(".")[0]
            cur_col = cur_col.split("_")
            print(cur_col)
            if len(cur_col)<=2:
                new_col = cur_col[0] + "_" + cur_col[1]
            elif len(cur_col) ==3:
                new_col = cur_col[0] + "_" + cur_col[1] + "_" + cur_col[2]
            else:
                cur_col = cur_col[2:]
                print(cur_col)
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
