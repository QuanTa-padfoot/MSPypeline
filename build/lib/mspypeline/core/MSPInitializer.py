import os
from copy import deepcopy
from inspect import isclass
from typing import Tuple, Dict, Type, Optional
try:
    from ruamel_yaml import YAML
except ModuleNotFoundError:
    from ruamel.yaml import YAML
import logging

from mspypeline.helpers import get_logger
from mspypeline import path_package, path_package_config, version
from mspypeline.file_reader import MissingFilesException, BaseReader


class MSPInitializer:
    """
    | An initializer class which is responsible for creating the directory to save the default YAML configuration
      file as well as reading and saving the specified settings.
    | The initializer also operates as a means of passing stored configurations to the plotter classes.
    """
    # set all file names that are required
    yml_file_name_tmp = "config_tmp.yml"
    yml_file_name = "config.yml"
    default_yml_name = "ms_analysis_default.yml"
    go_path = "go_terms"
    pathway_path = "pathways"
    possible_gos = sorted([x for x in os.listdir(os.path.join(path_package_config, go_path))
                           if x.endswith(".txt")])
    possible_pathways = sorted([x for x in os.listdir(os.path.join(path_package_config, pathway_path))
                                if x.endswith(".txt")])

    def __init__(self, path: str, file_path_yml: Optional[str] = None, loglevel=logging.DEBUG):
        """
        Parameters
        ----------
        path
            location where the directory/txt folder to the data can be found.
        file_path_yml
            path to the yaml config file
        loglevel
            level of the logger
        """
        self.logger = get_logger(self.__class__.__name__, loglevel=loglevel)
        # create a yaml file reader
        self.yaml = YAML()
        self.yaml.indent(mapping=2, sequence=4, offset=2)
        # self.yaml.indent(offset=2)
        self.yaml.default_flow_style = False
        self.yaml.width = 4096

        # attributes that change upon changing the starting dir
        #: configurations for the run. also saved configurations for the reader under the respective reader name
        self.configs = {}
        self.reader_data = {}

        self.interesting_proteins, self.go_analysis_gene_names = None, None

        # properties
        self._start_dir = None
        self._file_path_yaml = None

        # list to store all selectable terms; custom and provided
        self.list_full_gos = []
        self.list_full_pathways = []

        # set the specified dirs
        self.start_dir = path
        if file_path_yml is not None:
            self.file_path_yaml = file_path_yml

    @property
    def start_dir(self):
        return self._start_dir

    @start_dir.setter
    def start_dir(self, start_dir):
        start_dir = os.path.normpath(start_dir)
        # make sure to be on the right level and set starting dir
        if os.path.split(start_dir)[1] == "txt":
            self.logger.debug("Removing txt ending from path")
            self._start_dir = os.path.split(start_dir)[0]
        else:
            self._start_dir = start_dir
        self.logger.info(f"Starting dir: {self.start_dir}")
        # set all attributes back None that where file specific
        self.configs = {}
        self.reader_data = {}
        self.file_path_yaml = "file"
        # see if any custom lists can be found
        self.list_full_gos = []
        self.list_full_gos += MSPInitializer.possible_gos
        try:
            self.list_full_gos += [x for x in os.listdir(os.path.join(self._start_dir, "go_terms"))
                                   if os.path.isfile(os.path.join(self._start_dir, "go_terms", x))]
        except FileNotFoundError:
            pass

        self.list_full_pathways = []
        self.list_full_pathways += MSPInitializer.possible_pathways
        try:
            self.list_full_pathways += [x for x in os.listdir(os.path.join(self._start_dir, "pathways"))
                                        if os.path.isfile(os.path.join(self._start_dir, "pathways", x))]
        except FileNotFoundError:
            pass

    @property
    def path_config(self):
        return os.path.join(self.start_dir, "config")

    @property
    def file_path_yaml(self):
        """
        Setting the yaml file path will set the configurations of the class to the ones specified in the file.

        Note
        -----
        The value can be set to either:

        - "default"
        - "file"
        - a path to a yml file

        Raises
        ------
        ValueError
            if no valid value was provided
        FileNotFoundError
            if the file specified by the file_path_yml was not found

        """
        return self._file_path_yaml

    @file_path_yaml.setter
    def file_path_yaml(self, file_path_yml: str):
        if file_path_yml.lower() == "default":
            self._file_path_yaml = self.get_default_yml_path()
        elif file_path_yml.lower() == "file":
            if self.has_yml_file():
                self._file_path_yaml = os.path.join(self.start_dir, "config", MSPInitializer.yml_file_name)
            else:
                self._file_path_yaml = self.get_default_yml_path()
        elif file_path_yml.lower().endswith(('.yml', '.yaml')):
            self._file_path_yaml = os.path.normpath(file_path_yml)
        else:
            raise ValueError(f"Invalid value provided for yaml file: {file_path_yml}")
        self.logger.debug("yml file location: %s", self._file_path_yaml)

        # load the config from the yml file
        self.logger.info("loading yml file")
        with open(self.file_path_yaml) as f:
            self.configs = self.yaml.load(f)
        self.logger.debug(f"Config file contents: {self.configs}")

    def init_config(self):
        """
        | Creates the directory to save the configuration file if not present, updates and saves the configuration.
        | The function is usually applied to ensure that configs are provided to the initializer in order to avoid
          problems initializing the file reader with :meth:`read_data`.
        """
        os.makedirs(self.path_config, exist_ok=True)
        self.update_config_file()

    def has_yml_file(self) -> bool:
        if not os.path.isdir(self.start_dir):
            return False
        if "config" in os.listdir(self.start_dir):
            self.logger.debug("Found config dir")
            config_dir = os.path.join(self.start_dir, "config")
            if MSPInitializer.yml_file_name in os.listdir(config_dir):
                self.logger.debug("Found config.yml file in config dir")
                return True
        return False

    def get_default_yml_path(self) -> str:
        self.logger.debug("Loading default yml file from: %s, since 'default' or no (valid) file was selected",
                          path_package)
        return os.path.join(path_package_config, MSPInitializer.default_yml_name)

    def init_interest_from_txt(self) -> Tuple[Dict[str, list], Dict[str, list]]:
        dict_pathway = {}
        dict_go = {}
        for pathway in self.configs.get("pathways"):
            name, proteins = self.read_config_txt_file(pathway)
            dict_pathway[name] = proteins

        for go in self.configs.get("go_terms"):
            name, proteins = self.read_config_txt_file(go, False)
            dict_go[name] = proteins
        return dict_pathway, dict_go

    def read_config_txt_file(self, file, is_pathway: bool = True) -> Tuple[str, list]:
        path_full = os.path.join(self.start_dir, "pathways" if is_pathway else "go_terms", file)
        if not os.path.isfile(path_full):
            path_full = os.path.join(path_package_config, "pathways" if is_pathway else "go_terms", file)
            if not os.path.isfile(path_full):
                raise FileNotFoundError(f"The selected file: {file} cannot be found.")
        name = file.replace(".txt", "")
        with open(path_full) as f:
            proteins = []
            for line in f:
                proteins.append(line.strip())
        return name, proteins

    def update_config_file(self):
        # store the config file as tmp
        self.logger.debug("Updating yml settings file")
        yml_file_loc_tmp = os.path.join(self.path_config, MSPInitializer.yml_file_name_tmp)
        with open(yml_file_loc_tmp, "w") as outfile:
            self.yaml.dump(self.configs, outfile)

        # delete non tmp if exists
        yml_file_loc = os.path.join(self.path_config, MSPInitializer.yml_file_name)
        if MSPInitializer.yml_file_name in os.listdir(self.path_config):
            os.remove(yml_file_loc)

        # rename to non tmp
        os.rename(yml_file_loc_tmp, yml_file_loc)

    def read_data(self):
        """
        | Initiates the file reader by providing the directory and the configs to the reader.
        | The configs for the reader are taken from the configs from the name of the reader as key. E.g. mqreader.
        | In turn a :class:`~mspypeline.helpers.Utils.DataDict` is generated to provide the mapping to the input data
          (*reader_data*) for the further analysis with the :ref:`mspypeline plotters <plotters>`.
        """
        for Reader in BaseReader.__subclasses__():
            Reader: Type[BaseReader]  # for IDE hints
            try:
                reader = Reader(self.start_dir, self.configs.get(Reader.name, {}))
                self.configs[str(Reader.name)] = deepcopy(reader.reader_config)
                self.reader_data[Reader.name] = reader.full_data

            except MissingFilesException:
                self.logger.debug("No files found for reader: %s", Reader.name)

        # read all proteins and receptors of interest from the config dir
        self.logger.info("Reading pathway and GO list of interest")
        self.interesting_proteins, self.go_analysis_gene_names = self.init_interest_from_txt()
        self.update_config_file()
