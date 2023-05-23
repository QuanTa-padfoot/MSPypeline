import argparse
import os
import tkinter as tk
import tkinter.ttk as ttk
import tkinter.tix as tix
from tkinter import filedialog
import logging
from threading import Thread
from typing import Optional, Iterable
import webbrowser
try:
    from ruamel_yaml import YAML
except ModuleNotFoundError:
    from ruamel.yaml import YAML
from mspypeline import create_app
from mspypeline.core import MSPInitializer
from mspypeline.helpers import get_logger
from mspypeline.modules import default_normalizers
from mspypeline.file_reader import BaseReader, MQReader


class UIHandler:
    """
    | Used to take the mapping of arguments provided by the :class:`~MSPParser` through the command line to organize
      through which point of entry the data analysis should be performed.
    | The UIHandler can initialize the creation of the GUI, perform the whole analysis instantly by running the
      pipeline and creating all plots according to the configs or (not yet available) host the mspypeline on a flask
      server.
    """
    def __init__(self, file_dir, yml_file=None, gui=False, host_flask=False, selected_reader=MQReader.MQReader,
                 loglevel=logging.DEBUG, configs: dict = None):
        """
        Parameters
        ----------
        file_dir
            location where the directory/txt folder to the data can be found.
        yml_file
            path to the yaml config file
        gui
            should a GUI be compiled
        host_flask
            currently not implemented (should a flsk server be hosted)
        selected_reader
            instance of an :class:`~BaseReader` used to process data to an internal format for the plotter.
        loglevel
            level of the logger
        configs
            mapping of configuration
        """
        base_config = {
            "has_techrep": False,
            "has_groups": False,
        }
        if configs is None:
            configs = {}
        base_config.update(**configs)
        configs = base_config

        if gui and host_flask:
            raise ValueError("Can only specify one of host_flask and gui")
        if gui:
            app = MSPGUI(file_dir=file_dir, yml_file=yml_file, loglevel=loglevel, configs=configs)
            app.mainloop()
        elif host_flask:
            # TODO pass arguments to create app
            app = create_app()
            app.run(debug=True)
        else:
            mspinit = MSPInitializer(file_dir, yml_file, loglevel=loglevel)
            mspinit.init_config()
            mspinit.configs.update(configs)
            mspinit.read_data()
            # create plotter from initializer
            mspplots = selected_reader.plotter.from_MSPInitializer(mspinit)
            # create all plots and other results
            mspplots.create_results()
    

class MSPGUI(tk.Tk):
    """
    | This class is called from the :class:`~UIHandler` to create the GUI.
    """    
    def __init__(self, file_dir, yml_file=None, loglevel=logging.DEBUG, configs: dict = None):
        """
        Parameters
        ----------
        file_dir
            location where the directory/txt folder to the data can be found.
        yml_file
            path to the yaml config file
        loglevel
            level of the logger
        configs
            mapping containing the configurations
        """
        super().__init__()
        self.yaml_options = ["default"]
        self.reader_options = {reader.name: reader for reader in BaseReader.__subclasses__()}
        self.selected_reader = MQReader.MQReader
        self.normalize_options = ["None"] + list(default_normalizers.keys())
        self.mspinit = MSPInitializer(file_dir, yml_file, loglevel=loglevel)
        self.logger = get_logger(self.__class__.__name__, loglevel=loglevel)

        #icon = tk.PhotoImage(file = (os.path.dirname(os.path.abspath(__file__)) + '/GUI/mspypeline_ico.png'))
        #self.iconphoto(False, icon)

        
        #self.tk.call('source', (os.path.dirname(os.path.abspath(__file__)) + '/GUI/forest-light.tcl'))
        #self.tk.call('source', (os.path.dirname(os.path.abspath(__file__)) + '/GUI/forest-dark.tcl'))
        self.tk.call('source', (os.path.dirname(os.path.abspath(__file__)) + '/GUI/azure.tcl'))
        
        #style.theme_use('forest-light')
        #style.theme_use('forest-dark')
        self.tk.call("set_theme", "light")
        style = ttk.Style()
        style.configure('my.TButton',font=('Helvetica', 12, 'bold'))

        self.number_of_plots = 0
        self.plots_per_section = 0

        self.plot_settings = {}
        self.intensity_options = ["lfq_log2", "raw_log2", "ibaq_log2"]
        #,"lfq_normalized_log2", "raw_normalized_log2", "ibaq_normalized_log2]

        self.title("mspypeline")

        path_label = ttk.Label(self, text="Dir to analyze", font="Helvetica 10 bold").grid(row=0, column=0, sticky=tk.W, padx = 5)

        yaml_label = ttk.Label(self, text="Yaml file", font="Helvetica 10 bold").grid(row=2, column=0, sticky=tk.W, padx = 5)

        reader_label = ttk.Label(self, text="File reader", font="Helvetica 10 bold").grid(row=2, column=1, sticky=tk.W, padx = 5)

        self.dir_text = tk.StringVar(value=file_dir)
        dir_button = ttk.Button(self, textvariable=self.dir_text,
                                width=len(max(self.mspinit.list_full_gos, key=len)),
                                command=lambda: browsefunc(filedialog.askdirectory, self.dir_text, fn_params={
                                   "title": "Please select a directory with result files"}))
        create_tool_tip(dir_button, "Select a directory to analyze")
        dir_button.grid(row=1, column=0, sticky=tk.W, padx=20, columnspan=3)

        self.yaml_text = tk.StringVar()
        self.yaml_button = ttk.OptionMenu(self, self.yaml_text, *self.yaml_options)
        create_tool_tip(self.yaml_button, "Leave blank for first analysis"
                                            "\nFor re-analysis, load previously created Yaml file")
        self.yaml_button.grid(row=3, column=0, sticky=tk.W, padx=20)

        self.reader_text = tk.StringVar(value="mqreader")
        self.reader_button = ttk.OptionMenu(self, self.reader_text, list(self.reader_options.keys())[0], *self.reader_options.keys())
        self.reader_button.grid(row=3, column=1, sticky=tk.W, padx=20)

        self.replicate_var = tk.IntVar(value=1)
        replicate_button = ttk.Checkbutton(self, text="Does the file have technical replicates?",
                                          variable=self.replicate_var)
        replicate_button.grid(row=3, column=2, padx=15, pady=10, sticky=tk.W)
        create_tool_tip(replicate_button, "If selected, the samples of the last level are averaged")
        

        go_proteins_label = ttk.Label(self, text="Go analysis lists")
        create_tool_tip(go_proteins_label, "For a full list of Proteins see the Documentation."
                                           "\nCustom Gene Lists -> GO Terms")
        go_proteins_label.grid(row=5, column=0, sticky=tk.W, padx=5)
        # Button for submitting customized GO gene lists
        GO_button = ttk.Button(self, text='Upload file', width=10,
                               command=lambda: self.upload_file(list_type='GO'))
        create_tool_tip(GO_button, "Upload your GO gene list")
        GO_button.grid(row=5, column=1, sticky=tk.W)

        pathways_label = ttk.Label(self, text="Pathway analysis lists")
        create_tool_tip(pathways_label, "For a full list of Proteins see the Documentation."
                                           "\nCustom Gene Lists -> Pathways")
        pathways_label.grid(row=7, column=0, sticky=tk.W, padx=5)
        # Button for submitting customized pathway gene lists
        pathways_button = ttk.Button(self, text='Upload file',width=10,
                                command=lambda: self.upload_file(list_type='pathways'))
        create_tool_tip(pathways_button, "Upload your pathways gene list")
        pathways_button.grid(row=7, column=1, sticky=tk.W)

        design_label = ttk.Label(self, text="Sample names")
        create_tool_tip(design_label, "Inferred sample names of the experiment")
        design_label.grid(row=9, column=0, sticky=tk.W, padx=5)

        self.go_term_list = tk.Listbox(self, selectmode="multiple", height=5,
                                       width=len(max(self.mspinit.list_full_gos, key=len)))
        self.go_term_list.configure(exportselection=False)
        for x in self.mspinit.list_full_gos:
            self.go_term_list.insert("end", x)

        
        #scrollbar = ttk.Scrollbar(orient=tk.VERTICAL, command=self.go_term_list.yview)
        #self.go_term_list.config(yscrollcommand=scrollbar.set)
        #scrollbar.grid(row=4, column=0, columnspan=3, sticky=tk.W, padx=20, )

        self.go_term_list.grid(row=6, column=0, columnspan=3,sticky=tk.W, padx=20)

        self.pathway_list = tk.Listbox(self, selectmode="multiple", height=5,
                                       width=len(max(self.mspinit.list_full_gos, key=len)))
        self.pathway_list.configure(exportselection=False)
        for x in self.mspinit.list_full_pathways:
            self.pathway_list.insert("end", x)

        self.pathway_list.grid(row=8, column=0, columnspan=3, sticky=tk.W, padx=20)

        self.experiments_list = tk.Listbox(self, height=5, width=30)
        self.experiments_list.grid(row=10, column=0, columnspan=2, sticky=tk.W, padx=20)

        report_button = ttk.Button(self, text="Create Report",
                                  command=lambda: self.report_button())
        report_button.grid(row=11, column=0, padx=20, pady=20)
        create_tool_tip(report_button, "MaxQuant report for quality control")

        #plot_label = ttk.Label(self, text="Plot selection", font='Helvetica 10 bold').grid(row=14, column=0, sticky = tk.W)

        #intensity_label = ttk.Label(self, text="Intensities").grid(row=7, column=1)

        #levels_label = ttk.Label(self, text="Levels").grid(row=7, column=2)
        
        # Button to reorder the level
        reorder_level_button = ttk.Button(self, text="Reorder level", command=lambda:self.reorder_level())
        reorder_level_button.grid(row=9, column=2, sticky=tk.W, padx=20)
        create_tool_tip(reorder_level_button,"Reorder level of the samples and save on sample_mapping.txt, currently only available for SpectroReader")

        self.p_val_var = tk.IntVar(value=1)
        pval_button = ttk.Checkbutton(self, text="Use adjusted p value", variable=self.p_val_var)
        pval_button.grid(row=10, column=2, sticky=tk.W, padx=15)
        create_tool_tip(pval_button,"Select to focus on regulated proteins, deselect to focus on affected pathways and processes")

        update_button = ttk.Button(self, text="Update", command=lambda: self.update_button(), width=10)
        update_button.grid(row=11, column=1)
        create_tool_tip(update_button, "Press if Yaml or sample_mapping.txt files were changed")

        start_button = ttk.Button(self, text="Start", style='my.TButton',
                                 command=lambda: self.start_button())
        start_button.grid(row=12, column=3, sticky=tk.SE, padx=20, pady=10)

        documentation_link = ttk.Label(self, text="Documentation link", font="Helvetica 10 underline",
                                      cursor="hand2")
        documentation_link.bind("<ButtonRelease-1>",
                                lambda _: webbrowser.open_new("https://mspypeline.readthedocs.io/en/latest/"))
        documentation_link.grid(row=11, column=2)

        #self.running_text = tk.StringVar(value="Please press Start")
        #self.running_label = ttk.Label(self, textvariable=self.running_text).grid(row=12, column=3, sticky=tk.NE, padx=30)

        style_switch = ttk.Checkbutton(self, text='Dark theme', style='Switch.TCheckbutton', command=lambda: self.style_handler())
        style_switch.grid(row=12, column=0)


        self.heading_length = 12

        tabControl = ttk.Notebook(self, height=450)
        tab1 = ttk.Frame(tabControl)
        tab2 = ttk.Frame(tabControl)
        tab3 = ttk.Frame(tabControl)
        tabControl.add(tab1, text ='Normalization')
        tabControl.add(tab2, text ='Outlier Detection')
        tabControl.add(tab3, text ='Statistical Inference')
        tabControl.grid(row = 1, column=3, rowspan=self.heading_length, sticky=tk.NE, padx=20)

        ##Section Normalization

        self.plots_per_section = 0
        
        #ttk.Label(tab1, text="Normalization plots", font="Helvetica 10 bold").grid(
        #    row=self.heading_length + self.number_of_plots, column=0, sticky=tk.W, padx=5)
        #self.number_of_plots += 1
        norm_method_label = ttk.Label(tab1, text="Choose a Normalization Method:", font="Helvetica 10 bold")
        norm_method_label.grid(row=self.heading_length + self.number_of_plots+3, column=0, pady=20)
        create_tool_tip(norm_method_label, "For more information about normalization visit the documentation.\n"
                        "Can be found in Workflow -> Data Preprocessing -> Normalization options.")
        self.plots_per_section +=2
        self.number_of_plots += 1

        self.plot_row("Normalization overview", "normalization_overview_all_normalizers",
                      "Displays original data and data after different normalization methods", tab=tab1)
        self.plots_per_section +=1
        self.plot_row("Heatmap overview", "heatmap_overview_all_normalizers",
                      "Displays intensity heatmap of original data and after different normalization methods", tab=tab1)
        self.plots_per_section +=1

        self.plot_intermediate_row("Choose a Normalization Method", tab = tab1)


        ##Section Outlier detection
        self.plots_per_section = 0
        #ttk.Label(tab2, text="Outlier detection / Comparisons", font="Helvetica 10 bold").grid(
        #    row=self.heading_length + self.number_of_plots, column=0, sticky=tk.W, padx=5)
        #self.number_of_plots += 1
        self.plot_row("Detection counts", "detection_counts",
                      "How many proteins were detected how frequently in the samples of a group?", tab = tab2)
        self.plots_per_section +=1
        self.plot_row("Number of detected proteins", "detected_proteins_per_replicate",
                      "How many proteins were detected in each of my samples and in total for each group?", tab = tab2)
        self.plots_per_section +=1
        self.plot_row("Venn diagrams", "venn_results",
                      "How large is the intersection of detected proteins of my samples in each group?\n How many proteins are uniquely detected in a sample?", tab = tab2)
        self.plots_per_section +=1
        self.plot_row("Group diagrams", "venn_groups",
                      "How large is the intersection of detected proteins between different groups?\n How many proteins are uniquely detected in a group?", tab = tab2)
        self.plots_per_section +=1
        self.plot_row("PCA overview", "pca_overview",
                      "How similar are my samples? Do samples cluster together?", tab = tab2)
        self.plots_per_section +=1
        self.plot_row("Intensity histogram", "intensity_histograms",
                      "How does the intensity profile of my samples look? How similar are the intensity profiles?", tab = tab2)
        self.plots_per_section +=1
        self.plot_row("Relative std", "relative_std",
                      "What is the relative standard deviation of the samples of a group?", tab = tab2)
        self.plots_per_section +=1
        self.plot_row("Scatter replicates", "scatter_replicates",
                      "How well do the overall protein intensities of the samples of each group correlate?", tab = tab2)
        self.plots_per_section +=1
        self.plot_row("Experiment comparison", "experiment_comparison",
                      "How well do the overall protein intensities of different groups correlate?", tab = tab2)
        self.plots_per_section +=1
        self.plot_row("Rank", "rank",
                      "Where do my proteins of interest rank in intensity compared to all other proteins?", tab = tab2)

        ##Section Statistical Inference
        
        #ttk.Label(tab3, text="Statistical inference", font="Helvetica 10 bold").grid(
        #    row=self.heading_length + self.number_of_plots, column=0, sticky=tk.W, padx=5)
        #self.number_of_plots += 1
        self.plots_per_section = 0
        self.plot_row("Pathway Analysis", "pathway_analysis",
                      "What is the intensity of my proteins of interest, and is it significantly different in one group versus the other?", tab = tab3)
        #self.plot_row("Pathway Timecourse", "pathway_timecourse")
        self.plot_row("Go analysis", "go_analysis",
                      "Are the proteins of a group enriched for the selected GO terms?", tab = tab3)
        self.plot_row("Volcano plot (R)", "r_volcano",
                      "Which proteins are significantly higher or lower in intensity comparing two groups?\n Which proteins are detected only in one group and not in the other?", tab = tab3)
        # button for selecting samples to plot volcanoes
        volcano_label = ttk.Label(tab3, text="Select samples for volcano plot:", font="Helvetica 10 bold")
        volcano_label.grid(row=self.heading_length + self.number_of_plots, column=0, pady=20)
        self.customize_sample_button(plot_text="volcano", tab=tab3)        
        
        self.plot_row("Time-course intensities (R)", "r_timecourse",
                      "What is the dynamics of the protein level across several condition?\n Genes to be plotted are detected from selected GO and Pathway gene lists",
                      tab=tab3)
        # button for customizing sample selection for timecourse FC
        norm_method_label = ttk.Label(tab3, text="Settings for plotting timecourse:", font="Helvetica 10 bold")
        norm_method_label.grid(row=self.heading_length + self.number_of_plots, column=0, pady=20)
        self.customize_sample_button(plot_text="timecourse", tab=tab3)
        total_length = self.heading_length + self.number_of_plots
        total_length = self.heading_length + self.number_of_plots

        
        # add all tracing to the variables
        self.dir_text.trace("w", self.dir_setter)
        self.yaml_text.trace("w", self.yaml_path_setter)
        self.reader_text.trace("w", self.reader_setter)
        # make the GUI resizable
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        
        self.protocol("WM_DELETE_WINDOW", self._exit)
        self.mspinit.configs.update(configs)
        self.update_yaml_options()
        self.resizable(False, False) 

    def upload_file(self, list_type: str):
        '''To upload custom gene list for GO and pathways analysis'''
        filenames = filedialog.askopenfilenames(initialdir=self.dir_text,filetypes=[('Text Files', '*.txt')])
        targetPath = os.path.realpath(__file__)
        targetPath = targetPath.split(sep='mspypeline\core')[0]
        if list_type == 'GO':
            targetPath = targetPath + 'mspypeline\config\go_terms'
        elif list_type == 'pathways':
            targetPath = targetPath + 'mspypeline\config\pathways'
        for filedir in filenames:
            if filedir is not None:
                filename = filedir.replace('\\', '/')
                filename_nodir = filename.split('/')
                filename_nodir = filename_nodir[len(filename_nodir) - 1]
                target_filedir = os.path.join(targetPath, filename_nodir)
                target_filedir = target_filedir.replace('\\', '/')

                # Create a list of gene names, remove duplicates and sort the list
                gene_list = []
                for line in open(filename, 'r').readlines():
                    if line not in gene_list and line != '\n':
                        gene_list.append(line)

                # Write gene_list to the new file
                with open(target_filedir, 'w') as target:
                    for line in gene_list:
                        target.write(line)
    
    def style_handler(self):
        if self.tk.call("ttk::style", "theme", "use") == "azure-dark":
            # Set light theme
            self.tk.call("set_theme", "light")
            style = ttk.Style()
            style.configure('my.TButton',foreground="black", font=('Helvetica', 12, 'bold'))
        else:
            # Set dark theme
            style = ttk.Style()
            self.tk.call("set_theme", "dark")
            style.configure('my.TButton',foreground="white", font=('Helvetica', 12, 'bold'))

    def _exit(self):
        self.quit()
        self.destroy()

    def dir_setter(self, *args):
        start_dir = self.dir_text.get()
        if os.path.split(start_dir)[1] == "txt":
            start_dir = os.path.split(start_dir)[0]
        self.mspinit.start_dir = start_dir
        self.update_yaml_options()

    def update_yaml_options(self):
        if self.mspinit.has_yml_file():
            self.yaml_text.set("file")
            self.yaml_options = ["default", "file"]
            # next steps should also be done hy what the update button does?
        else:
            self.yaml_text.set("default")
            self.yaml_options = ["default"]
        self.yaml_button["menu"].delete(0, "end")
        for op in self.yaml_options:
            self.yaml_button["menu"].add_command(label=op, command=tk._setit(self.yaml_text, op))

    def update_listboxes(self):
        # delete all experiments then add from file (try sample_mapping.txt first, if sample_mapping is not present then add from config)
        self.experiments_list.delete(0, "end")
        mapping_txt = os.path.join(self.mspinit.start_dir,"config/sample_mapping.txt")
        try:
            with open(mapping_txt, "r") as f:
                next(f)
                for line in f.readlines():
                    op = line.split('\t')[1]
                    op = op[:-1]
                    self.experiments_list.insert("end", op)
                f.close()
        except FileNotFoundError:
            for op in self.mspinit.configs.get(self.selected_reader.name, {}).get("all_replicates", []):
                self.experiments_list.insert("end", op)
        # clear selection then select from configs
        for i, pathway in enumerate(self.mspinit.list_full_pathways):
            self.pathway_list.select_clear(i)
        if self.mspinit.configs.get("pathways"):
            for pathway in self.mspinit.configs.get("pathways"):
                try:
                    self.pathway_list.select_set(self.mspinit.list_full_pathways.index(pathway))
                except ValueError:
                    self.logger.warning("Selected pathway file %s not found", pathway)
        # clear selection then select from configs
        for i, go in enumerate(self.mspinit.list_full_gos):
            self.go_term_list.select_clear(i)
        if self.mspinit.configs.get("go_terms"):
            for go in self.mspinit.configs.get("go_terms"):
                try:
                    self.go_term_list.select_set(self.mspinit.list_full_gos.index(go))
                except ValueError:
                    self.logger.warning("Selected go term file %s not found", go)

    def yaml_path_setter(self, *args):
        self.mspinit.file_path_yaml = self.yaml_text.get()
        # get the reader class by the saved name
        self.selected_reader = self.reader_options.get(self.mspinit.configs.get("selected_reader", "mqreader"))
        reader_settings = self.mspinit.configs.get(self.selected_reader.name, {})
        self.reader_text.set(self.selected_reader.name)
        level_names = reader_settings.get("level_names", [])
        level_names = {i: name for i, name in enumerate(level_names)}
        levels = reader_settings.get("levels", 0)
        for plot_name in self.selected_reader.plotter.possible_plots:
            plot_settings_name = plot_name + "_settings"
            plot_settings = self.mspinit.configs.get(plot_settings_name, {})
            var_name = plot_name.replace("plot_", "") + "_var"
            int_name = var_name.replace("var", "int")
            levels_name = var_name.replace("var", "levels")
            # update all settings in the GUI
            self.plot_settings[int_name].set(plot_settings.get("create_plot", False))

            plot_intensities = plot_settings.get("dfs_to_use", [])
            self.plot_settings[var_name].update_selection(plot_intensities)

            plot_levels = plot_settings.get("levels", [])
            selected_levels = self.plot_settings[levels_name]
            selected_levels.update_options([level_names.get(l, l) for l in range(levels)])
            selected_levels.update_selection([level_names.get(pl, pl) for pl in plot_levels])
        self.replicate_var.set(self.mspinit.configs.get("has_techrep", True))
        self.p_val_var.set(self.mspinit.configs.get("plot_r_volcano_settings", {}).get("adj_pval", False))
        self.normalizer_text.set(self.mspinit.configs.get("selected_normalizer", "None"))
        self.update_listboxes()

    def reader_setter(self, *args):
        self.selected_reader = self.reader_options[self.reader_text.get()]

    def update_button(self):
        self.mspinit.configs["has_techrep"] = bool(self.replicate_var.get())
        self.mspinit.configs["selected_reader"] = str(self.reader_text.get())
        self.mspinit.configs["selected_normalizer"] = str(self.normalizer_text.get())
        yamlReader = YAML()
        with open(self.mspinit.get_default_yml_path()) as f:
            defaultConfigs = yamlReader.load(f)
        reader_settings = self.mspinit.configs.get(self.selected_reader.name, {})
        level_names = reader_settings.get("level_names", [])
        level_names = {name: i for i, name in enumerate(level_names)}
        for plot_name in self.selected_reader.plotter.possible_plots:
            plot_settings = plot_name + "_settings"
            var_name = plot_name.replace("plot_", "") + "_var"
            int_name = var_name.replace("var", "int")
            levels_name = var_name.replace("var", "levels")
            selected_settings = {
                "create_plot": bool(self.plot_settings[int_name].get()),
                "dfs_to_use": [k for k, v in self.plot_settings[var_name].get_selection().items() if v],
                "levels": [level_names.get(k, k) for k, v in self.plot_settings[levels_name].get_selection().items() if v]
            }
            additional_settings = self.mspinit.configs.get(plot_settings, {})
            additional_settings = {k: v for k, v in additional_settings.items()
                                   if k != "create_plot" and k != "dfs_to_use" and k != "levels"}
            selected_settings.update(additional_settings)
            for k, v in selected_settings.items():
                self.mspinit.configs[plot_settings][k] = v
        gos = self.go_term_list.curselection()
        gos = [self.mspinit.list_full_gos[int(go)] for go in gos]
        pathways = self.pathway_list.curselection()
        pathways = [self.mspinit.list_full_pathways[int(pathway)] for pathway in pathways]
        self.mspinit.configs["plot_r_volcano_settings"]["adj_pval"] = bool(self.p_val_var.get())
        if not gos:
            gos = defaultConfigs["go_terms"]
        self.mspinit.configs["go_terms"] = gos
        if not pathways:
            pathways = defaultConfigs["pathways"]
        self.mspinit.configs["pathways"] = pathways
        self.mspinit.init_config()
        self.mspinit.read_data()
        self.update_listboxes()
        self.update_yaml_options()

    def start_button(self):
        self.warningbox = WarningBox('Status Update', 'Creating Plots')
        self.warningbox.wait_visibility()
        x = self.winfo_x() + self.winfo_width()//2 - self.warningbox.winfo_width()//2
        y = self.winfo_y() + self.winfo_height()//2 - self.warningbox.winfo_height()//2
        self.warningbox.geometry(f"+{x}+{y}")
        self.start_mspypeline_thread(None)
        #self.running_text.set("Please press Start")

    def start_ops(self):
        try:
            self.update()
            self.update_button()
            mspplots = self.selected_reader.plotter.from_MSPInitializer(self.mspinit)
            mspplots.create_results()
        except Exception as err:
            self.err = type(err)
            return
        else:
            self.err = None
    
    def start_mspypeline_thread(self, event):
        global start_thread
        start_thread = Thread(target=self.start_ops)
        self.warningbox.progressbar.start()
        start_thread.daemon = True
        start_thread.start()
        self.after(20, self.check_mspypeline_thread)
    
    def check_mspypeline_thread(self):
        if start_thread.is_alive():
            self.after(20, self.check_mspypeline_thread)
        elif self.err == None:
            self.popup_window('Status Update', 'Tasks completed')
            #self.warningbox.updateInfo('Status Update', 'Tasks completed')
        elif self.err == KeyError:
            self.popup_window(title='Status Update', message=('File could not be read with selected reader\nError code: ' + self.err.__name__), error=True)
            #self.warningbox.updateInfo('Status Update', 'File could not be read with selected reader')
        ### ADD HERE ERROR TYPES FOR PROMPT DISPLAY ###
        elif self.err == NotImplementedError:
            self.popup_window(title='Status Update', message=('A requested task is currently not implemented\nError code: ' + self.err.__name__), error = True)
        else:
            self.popup_window(title='Status Update', message=('An error occurred, please check Terminal\nError code: ' + self.err.__name__), error=True)
            #self.warningbox.updateInfo('Status Update', 'An error occured, please check Terminal')
    
    def report_button(self):
        #self.running_text.set("Please press Start")
        self.warningbox = WarningBox('Status Update', 'Creating Report')
        self.warningbox.wait_visibility()
        x = self.winfo_x() + self.winfo_width()//2 - self.warningbox.winfo_width()//2
        y = self.winfo_y() + self.winfo_height()//2 - self.warningbox.winfo_height()//2
        self.warningbox.geometry(f"+{x}+{y}")
        self.report_mspypeline_thread(None)
    
    def report_ops(self):
        try:
            self.update()
            self.update_button()
            mspplots = self.selected_reader.plotter.from_MSPInitializer(self.mspinit)
            mspplots.create_report()
        except Exception as err:
            self.err = type(err)
            return
        else:
            self.err = None

    def report_mspypeline_thread(self, event):
        global report_thread
        report_thread = Thread(target=self.report_ops)
        report_thread.daemon = True
        self.warningbox.progressbar.start()
        report_thread.start()
        self.after(20, self.check_report_thread)
    
    def check_report_thread(self):
        if report_thread.is_alive():
            self.after(20, self.check_report_thread)
        elif self.err == None:
            self.popup_window('Status Update', 'Report completed')
            #self.warningbox.updateInfo('Status Update', 'Tasks completed')
        elif self.err == KeyError:
            self.popup_window(title='Status Update', message=('File could not be read with selected reader\nError code: ' + self.err.__name__), error=True)
            #self.warningbox.updateInfo('Status Update', 'File could not be read with selected reader')
        ### ADD HERE ERROR TYPES FOR PROMPT DISPLAY ###
        else:
            self.popup_window(title='Status Update', message=('An error occurred, please check Terminal\nError code: ' + self.err.__name__), error=True)
            #self.warningbox.updateInfo('Status Update', 'An error occured, please check Terminal')
        
    def popup_window(self, title='', message='', error=False):
        self.warningbox.progressbar.stop()
        self.warningbox.destroy()
        if error:
            tk.messagebox.showerror(title=title, message=message)
        else:
            tk.messagebox.showinfo(title=title, message=message)

    def plot_intermediate_row(self, text: str, tab = None):
        if tab == None:
            col = 0
            row = self.heading_length + self.number_of_plots
            int_var = tk.IntVar(value=1)
            self.normalizer_text = tk.StringVar(value="None")
            self.normalizer_button = ttk.OptionMenu(self, self.normalizer_text, *self.normalize_options)
            self.normalizer_button.grid(row=row, column=col)

            self.number_of_plots += 1
        else:
            col = 1
            row = self.heading_length + self.number_of_plots
            int_var = tk.IntVar(value=1)
            self.normalizer_text = tk.StringVar(value="None")
            self.normalizer_button = ttk.OptionMenu(tab, self.normalizer_text, *self.normalize_options)
            self.normalizer_button.grid(row=row, column=col, sticky=tk.W, padx = 8)

            self.number_of_plots += 1

    def customize_sample_button(self, plot_text: str = None, tab = None):
        '''Create a button for customizing plots (currently available for volcano and timecourse)

        need to identify plot_text, can be either "volcano" or "timecourse".'''
        if tab == None:
            col = 0
            row = self.heading_length + self.number_of_plots
            if plot_text == "timecourse":
                customize = ttk.Button(self, text="Customize", command= lambda: self.customize_timecourse())
                customize.grid(row=row, column=col, sticky=tk.W, padx=5)
            elif plot_text == "volcano":
                customize = ttk.Button(self, text="Select sample", command= lambda: self.customize_volcano())
                customize.grid(row=row, column=col, sticky=tk.W, padx=5)

            self.number_of_plots += 1
        else:
            col = 1
            row = self.heading_length + self.number_of_plots
            if plot_text == "timecourse":
                customize = ttk.Button(tab, text="Customize", command= lambda: self.customize_timecourse())
                customize.grid(row=row, column=col, sticky=tk.W, padx=5)
            elif plot_text == "volcano":
                customize = ttk.Button(tab, text="Select sample", command= lambda: self.customize_volcano())
                customize.grid(row=row, column=col, sticky=tk.W, padx=5)

            self.number_of_plots += 1

    def customize_volcano(self):
        '''Popup window to select samples for plotting volcano.
        The window finds all samples at a given level, and one can choose them to be either "downregulated" or "upregulated".
        The top option, DEFAULT, means that MSPypeline will plot all sample combinations.
        Nevertheless, this does NOT mean the desired upregulated-downregulated sample pair will always be plotted.
        If default is selected, all other selected samples will be ignored.
        This function serves to avoid plotting unnecessary combinations.
        Selected settings are passed to configs, which is then picked up by BasePlotter to make volcano plots.'''
        window = tk.Toplevel()
        window.geometry("500x400")
        window.title("Select samples for plotting volcanoes")

        # get the level for volcano plot
        all_level = self.plot_settings.get("r_volcano_levels", []).get_selection()
        selected_level = [level for level in all_level.keys() if all_level[level]]
        if selected_level == [] or len(selected_level) > 1:
            note = tk.Label(window, text="No level or more than 1 level was selected. To select samples, only 1 level is allowed at a time")
            note.grid(row=0, column=0)
            note1 = tk.Label(window, text="Please close this window and try again after checking the level selection")
            note1.grid(row=1, column=0)
        else:
            selected_level = selected_level[0]
            try:
                all_replicates = []
                with open(os.path.join(self.mspinit.start_dir,"config/sample_mapping.txt"), 'r') as f:
                    next(f)
                    for line in f.readlines():
                        sample = line.split('\t')[1]
                        all_replicates.append(sample[:-1])
                    f.close()
            except FileNotFoundError:
                all_replicates = self.mspinit.configs.get(self.selected_reader.name, {}).get("all_replicates", [])
            max_level = len(all_replicates[0].split("_")) - 1
            all_sample = ["default"]
            if max_level == selected_level:
                note = tk.Label(window,
                                text="Volcano plot is not possible at the lowest level")
                note.grid(row=0, column=0)
                note1 = tk.Label(window,
                                 text="Please close this window and try again after checking the level selection")
                note1.grid(row=1, column=0)
            else:
                for rep in all_replicates:
                    rep_split = rep.split("_")
                    del rep_split[-(max_level-selected_level):]
                    rep_split = "_".join(rep_split)
                    if rep_split not in all_sample:
                        all_sample.append(rep_split)

                # layout for the popup windown
                sample1_label = tk.Label(window, text= "Select all 1st samples to compare\n (downregulated)")
                sample1_label.grid(row=0,column=0)
                sample1_list = tk.Listbox(window, selectmode="multiple", height=15, width=len(max(all_sample, key=len))+2)
                sample1_list.configure(exportselection=False)
                sample1_list.grid(row=1, column=0, sticky=tk.W, padx=20)

                filler_label = tk.Label(window, text= "  ")
                filler_label.grid(row=0, column=1)

                sample2_label = tk.Label(window, text="Select all 2nd samples to compare\n (upregulated)")
                sample2_label.grid(row=0, column=2)
                sample2_list = tk.Listbox(window, selectmode="multiple", height=15, width=len(max(all_sample, key=len)) + 2)
                sample2_list.configure(exportselection=False)
                sample2_list.grid(row=1, column=2, sticky=tk.W, padx=20)

                # Add a default option, if chosen then the all combinations will be plotted

                # Add all samples to both list
                for x in all_sample:
                    sample1_list.insert("end", x)
                    sample2_list.insert("end", x)

                # Add previous selection from configs
                s1 = self.mspinit.configs["plot_r_volcano_settings"].get("sample1_list")
                if s1 == "default":
                    sample1_list.select_set(0)
                    sample2_list.select_set(0)
                else:
                    try:
                        for sample in s1:
                            sample_index = all_sample.index(sample)
                            sample1_list.select_set(sample_index)
                        for sample in self.mspinit.configs["plot_r_volcano_settings"].get("sample2_list"):
                            sample_index = all_sample.index(sample)
                            sample2_list.select_set(sample_index)
                    except (ValueError, TypeError):
                        sample1_list.select_set(0)
                        sample2_list.select_set(0)
                def update_volcano_settings():
                    '''Triggered when press the OK button. This function saves the selected samples and close the popup window.
                    If default was chosen, all other selected samples will be ignored.
                    If nothing was chosen, the default option will be applied'''
                    sample1, sample2 = [], []
                    for i in sample1_list.curselection():
                        sample1.append(sample1_list.get(i))
                    for i in sample2_list.curselection():
                        sample2.append(sample2_list.get(i))
                    if "default" in sample1 or "default" in sample2:
                        sample1, sample2 = "default", "default"
                        print("The \"default\" option was selected, other selected samples will be ignored")
                    elif sample1 == [] or sample2 == []:
                        sample1, sample2 = "default", "default"
                        print("No sample was selected, the \"default\" option will be applied")
                    self.mspinit.configs["plot_r_volcano_settings"].update({
                        "sample1_list": sample1,
                        "sample2_list": sample2
                        })
                    print("Sample selection for volcano updated!")
                    window.destroy()

                def reset_selection():
                    '''Triggered when press the Reset button. All selected samples will be cleared and the default is selected'''
                    sample1_list.selection_clear(0, 'end')
                    sample2_list.selection_clear(0, 'end')
                    sample1_list.select_set(0)
                    sample2_list.select_set(0)

                resetButton = tk.Button(window, text="Reset",
                                    command=lambda: reset_selection())
                resetButton.grid(row=2, column=1, padx=5, sticky=tk.W)
                okButton = tk.Button(window, text="OK",
                                 command=lambda: update_volcano_settings())
                okButton.grid(row=2, column=3, padx=5, sticky=tk.W)
        window.mainloop()
            
    def customize_timecourse(self):
        '''Popup window to select samples for plotting the timecourse.
        Selected settings are passed to configs, which are then picked up by BasePlotter to make plots.'''
        window = tk.Toplevel()
        window.geometry("500x400")
        window.title("Selecting samples for plotting timecourse Fold Change")
        # get sample names from sample_mapping.txt; if not present then take it from config
        try:
            all_replicates = []
            with open(os.path.join(self.mspinit.start_dir, "config/sample_mapping.txt"), 'r') as f:
                next(f)
                for line in f.readlines():
                    sample = line.split('\t')[1]
                    all_replicates.append(sample[:-1])
                f.close()
        except FileNotFoundError:
            all_replicates = self.mspinit.configs.get(self.selected_reader.name, {}).get("all_replicates", [])

        # Remove the last two levels in sample display
        all_sample = []
        for i in range(len(all_replicates)):
            sample = all_replicates[i].split("_")
            del sample[-1]
            del sample[-1]
            new_sample = "_".join(sample)
            if new_sample not in all_sample and new_sample != "":
                all_sample.append(new_sample)

        if all_sample == []:
            note = tk.Label(window, text="Samples not found! Check that the directory to the data is correct")
            note.grid(row=0,column=0)
            note1 = tk.Label(window, text="In addition, samples must have at least 3 levels, \nwith the second to last one displaying timepoints")
            note1.grid(row=1, column=0)
        else:

            plot_errorbar = tk.IntVar(value=int(self.mspinit.configs["plot_r_timecourse_settings"].get("plot_errorbar")))
            errorbar_option = tk.Checkbutton(window, text="Plot errorbars", variable= plot_errorbar)
            errorbar_option.grid(row=0, column=0, sticky= tk.W)

            align_yaxis = tk.IntVar(value=int(self.mspinit.configs["plot_r_timecourse_settings"].get("align_yaxis")))
            align_yaxis_option = tk.Checkbutton(window, text="Align y axes", variable= align_yaxis)
            align_yaxis_option.grid(row=1, column=0, sticky= tk.W)

            normalizing_sample_label = tk.Label(window, text="Select a sample to normalize below for plotting fold change\nOr choose None for plotting protein intensities")
            normalizing_sample_label.grid(row=2, column=0)
            normalizing_sample = tk.Listbox(window, selectmode="single", height=5,
                                       width=len(max(all_sample, key=len))+2)
            normalizing_sample.configure(exportselection=False)
            normalizing_sample.grid(row=3, column=0, columnspan=2, sticky=tk.W, padx=20)

            time_match_norm = tk.IntVar(
                value=int(self.mspinit.configs["plot_r_timecourse_settings"].get("matching_time_normalization")))
            time_match_norm_option = tk.Checkbutton(window, text="Normalize by each time point instead of the first time point\nInvalid if does not select sample to normalize", variable=time_match_norm)
            time_match_norm_option.grid(row=4, column=0, sticky = tk.W)

            sample_list_label = tk.Label(window, text="Select all samples to plot below")
            sample_list_label.grid(row=5, column=0, sticky = tk.W)
            sample_list = tk.Listbox(window, selectmode="multiple", height=5,
                                        width=len(max(all_sample, key=len))+2)
            sample_list.configure(exportselection=False)

            sample_list.grid(row=6, column=0, columnspan=2, sticky=tk.W, padx=20)

            # Add a no normalizer option, if chosen then the intensities will be plotted instead of fold change
            normalizing_sample.insert("end", "None")
            # Add all samples to both list
            for x in all_sample:
                normalizing_sample.insert("end", x)
                sample_list.insert("end", x)

            # Add previous selection from configs to normalizing_sample and sample_list. If both are not present then pass
            try:
                s = self.mspinit.configs["plot_r_timecourse_settings"].get("sample_to_normalize")
                if s == 'None':
                    sample_index = 0
                else:
                    sample_index = all_sample.index(s)+1  # add 1 since a None option was added at the beginning of normalizing_sample
                normalizing_sample.select_set(sample_index)
                for sample in self.mspinit.configs["plot_r_timecourse_settings"].get("samples_to_plot"):
                    sample_index = all_sample.index(sample)
                    sample_list.select_set(sample_index)
            except ValueError:
                normalizing_sample.select_set(0)

            def update_timecourse_settings():
                samples_to_plot = []
                errorbar = plot_errorbar.get()
                timenorm = time_match_norm.get()
                alignyaxis = align_yaxis.get()
                for i in sample_list.curselection():
                    samples_to_plot.append(sample_list.get(i))
                self.mspinit.configs["plot_r_timecourse_settings"].update({
                    "plot_errorbar": bool(errorbar),
                    "align_yaxis": bool(alignyaxis),
                    "sample_to_normalize": normalizing_sample.get(normalizing_sample.curselection()),
                    "matching_time_normalization": bool(timenorm),
                    "samples_to_plot": samples_to_plot
                })
                print("Settings for plotting timecourse updated!")
                window.destroy()


            okButton = tk.Button(window, text="OK",
                                command= lambda: update_timecourse_settings())
            okButton.grid(row=7, column=1, padx=5, sticky=tk.W)

        window.mainloop()
        
    def plot_row(self, text: str, plot_name: str, plot_tool_tip: str = None, tab = None):
        row = self.heading_length + self.number_of_plots
        col = 0
        row = self.heading_length + self.number_of_plots
        int_var = tk.IntVar(value=1)
        if tab == None:
            checkbutton = ttk.Checkbutton(self, text=text, variable=int_var)
            checkbutton.grid(row=row, column=col, sticky=tk.W, padx=20)
            if plot_tool_tip:
                create_tool_tip(checkbutton, plot_tool_tip)
            intensity_list = MultiSelectOptionMenu(self, self.intensity_options, "Select Intensities")
            intensity_list.grid(row=row, column=0, sticky=tk.W, padx=5)

            level_list = MultiSelectOptionMenu(self, button_text="Select Levels")
            level_list.grid(row=row, column=0, sticky=tk.W, padx=5)

            self.number_of_plots += 1
            self.plot_settings.update({
                f"{plot_name}_int": int_var,
                f"{plot_name}_var": intensity_list,
                f"{plot_name}_levels": level_list
            })
        else:
            checkbutton = ttk.Checkbutton(tab, text=text, variable=int_var)
            checkbutton.grid(row=row, column=col, sticky=tk.W, padx=20)
            if plot_tool_tip:
                create_tool_tip(checkbutton, plot_tool_tip)
            intensity_list = MultiSelectOptionMenu(tab, self.intensity_options, "Select Intensities")
            intensity_list.grid(row=row, column=1+col, sticky=tk.W, padx=5)

            level_list = MultiSelectOptionMenu(tab, button_text="Select Levels")
            level_list.grid(row=row, column=2+col, sticky=tk.W, padx=5)

            self.number_of_plots += 1
            self.plot_settings.update({
                f"{plot_name}_int": int_var,
                f"{plot_name}_var": intensity_list,
                f"{plot_name}_levels": level_list
            })
    
    def reorder_level(self):
        """Reorder the level of the samples. Settings are saved in the text file sample_mapping.txt
        Called by clicking the reorder_level_button"""
        window = tk.Toplevel()
        window.title("Reordering level")
        mapping_txt = os.path.join(self.mspinit.start_dir,"config/sample_mapping.txt")

        try:
            with open(mapping_txt, "r") as f:
                next(f)
                original_name_example = f.readline().split('\t')[0]
                f.close()
        except FileNotFoundError:
            original_name_example = self.mspinit.configs.get(self.selected_reader.name, {}).get("all_replicates", [])[0]
        original_level_dict = {}
        n = 0
        for element in original_name_example.split('_'):
            original_level_dict[element] = n
            n += 1
        instruction_note = ttk.Label(window, text="You can change the level in analysis_design by typing the new name of a sample.\n"
                                                  "Other samples' names will be changed accordingly.\n"
                                                  "E.g. Name as in data file: Abc_Def_Ghi\n             Type new name: Ghi_Abc_Def\n   ")
        instruction_note.grid(row=0, column=0)
        orig_name_note = ttk.Label(window, text=f"Name as in the data file: {original_name_example}")
        orig_name_note.grid(row=1, column=0)
        new_name_note = ttk.Label(window, text=f"Enter the new name below")
        new_name_note.grid(row=2, column=0)

        inputtxt = tk.Entry(window)
        inputtxt.grid(row=3, column=0, padx=5, sticky=tk.W, pady=20)

        remove_config_note = ttk.Label(window, text="After pressing Confirm, mspypeline will recompile the config.yml file\n"
                                                    "This might take a while")
        remove_config_note.grid(row=4, column=0)

        def confirm_new_level():
            new_name_example = inputtxt.get().split('_')
            new_level_order = []
            if len(new_name_example) != n:
                error_note = tk.Label(window,
                                      text="Cannot update the new level. Please check the spelling and try again")
                error_note.grid(row=6, column=0)
                return
            try:
                for element in new_name_example:
                    new_level_order.append(original_level_dict[element])
                try:
                    with open(mapping_txt, 'r') as f:
                        all_replicates = []
                        next(f)
                        for line in f.readlines():
                            all_replicates.append(line.split('\t')[0])
                        f.close()
                except FileNotFoundError:
                    all_replicates = self.mspinit.configs.get(self.selected_reader.name, {}).get("all_replicates", [])
                with open(mapping_txt, 'w') as f:
                    f.write('Name as in data file\tNew name\n')
                    for name in all_replicates:
                        old_name = name.split('_')
                        new_name = [old_name[n] for n in new_level_order]
                        f.write(f"{name}\t{'_'.join(new_name)}\n")
                    f.close()
                # delete the old config file and compile it again
                os.remove(os.path.join(self.mspinit.start_dir, "config/config.yml"))
                self.update_button()
                print('New level order updated in sample_mapping.txt. Please delete the config.yml file and run Update again')
                window.destroy()
            except KeyError:
                error_note = tk.Label(window, text="Cannot update the new level. Please check the spelling and try again")
                error_note.grid(row=6, column=0)
        confirm_button = ttk.Button(window, text="Confirm", command=lambda: confirm_new_level())
        confirm_button.grid(row=5, column=0)

        window.geometry("500x300")
        window.mainloop()

class WarningBox(tk.Toplevel):
    def __init__(self, title='', message=''):
        tk.Toplevel.__init__(self)
        self.geometry('300x100')
        self.title(title)
        self.messageLabel = tk.Label(self, text=message)
        self.messageLabel.pack(expand=True, fill=tk.BOTH)
        self.progressbar = ttk.Progressbar(
            self,
            orient='horizontal',
            mode='indeterminate',
            length=150
        )
        self.progressbar.pack(anchor=tk.S, pady=(0,20))
        self.progressbar.start()

    def updateInfo(self, title='', message=''):
        self.progressbar.destroy()
        self.title(title)
        self.messageLabel.configure(text=message)
        exitbutton = ttk.Button(self, text='OK', command=lambda: self.destroy())
        exitbutton.pack(anchor=tk.SE, pady=(0,10), padx=(0,20))

    def errorbox(self, title='', message=''):
        self.progressbar.destroy()
        self.title(title)
        self.errorlogo = tk.PhotoImage(self, file = "Your_image.png")
        self.messageLabel.configure(text=message)
        exitbutton = ttk.Button(self, text='OK', command=lambda: self.destroy())
        exitbutton.pack(anchor=tk.SE, pady=(0,10), padx=(0,20))


class MultiSelectOptionMenu(tk.Frame):
    def __init__(self, parent, choices: Optional[Iterable] = None, button_text: str = "Default text"):
        super().__init__(parent)
        menubutton = ttk.Menubutton(self, text=button_text, style='custom.TMenubutton')
        self.menu = tk.Menu(menubutton, tearoff=False)
        menubutton.configure(menu=self.menu)
        menubutton.pack(pady=3, padx=3)
        self.choices_dict = {}
        self.choices = choices if choices is not None else ()
        self.update_options()

    def update_options(self, choices: Optional[Iterable] = None):
        if choices is not None:
            self.choices = choices
            self.choices_dict.clear()
            self.menu.delete(0, "end")
        for choice in self.choices:
            self.choices_dict[choice] = tk.BooleanVar(value=False)
            self.menu.add_checkbutton(label=choice, variable=self.choices_dict[choice], onvalue=True, offvalue=False)

    def update_selection(self, choices: Iterable):
        for choice in self.choices:
            if choice in choices:
                self.choices_dict[choice].set(True)
            else:
                self.choices_dict[choice].set(False)

    def get_selection(self):
        return {k: v.get() for k, v in self.choices_dict.items()}


class ToolTip:
    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None

    def showtip(self, text: str):
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 57
        y = y + cy + self.widget.winfo_rooty() + 27
        self.tipwindow = tk.Toplevel(self.widget)
        self.tipwindow.wm_overrideredirect(1)
        self.tipwindow.wm_geometry("+%d+%d" % (x, y))
        label = ttk.Label(self.tipwindow, text=text, justify=tk.LEFT,
                         background="#ffffe0", foreground='black',relief=tk.SOLID, borderwidth=1,
                         font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        if self.tipwindow:
            self.tipwindow.destroy()


def create_tool_tip(widget, text):
    tool_tip = ToolTip(widget)
    widget.bind('<Enter>', lambda x: tool_tip.showtip(text))
    widget.bind('<Leave>', lambda x: tool_tip.hidetip())


class MSPParser(argparse.ArgumentParser):
    """
    | Uses the ``argparse`` module to provide a parser for command line options, arguments and sub-commands through
      which the analysis can be started or the GUI can be called (see :ref:`get-started`).
    """
    def __init__(self):
        super().__init__(description="A pipeline to analyze result files from a MaxQuant report. "
                                     "The required result files are in the txt directory.")
        self.add_argument(
            '--dir',
            dest='file_dir',
            action='store',
            default="",
            help="Path to directory of analysis which should be analyzed."
                 "If not set the program will open a prompt and ask to select one."
        )
        self.add_argument(
            '--yml-file',
            dest='yml_file',
            action='store',
            default=None,
            help="Path to yml file which should be used for analysis, or 'default' / 'file'."
        )
        self.add_argument(
            "--loglevel",
            dest="loglevel",
            action="store",
            default=logging.WARNING,
            help="Logging level of analysis. Should be from options (lowest to highest): DEBUG < INFO < WARNING < ERROR. "
                 "The higher the logging level the fewer messages are shown. Default: WARNING"
        )
        self.add_argument(
            "--has-techrep",
            dest="has_techrep",
            default=False,
            const=True,
            nargs="?",
            help="If you have replicates of each experiment specify this"
                 "Replicates need to be enumerated in a way that numbers are added at the end of the name"
                 "If no replicates are in data set no venn diagrams will be generated"
        )
        self.add_argument(
            "--gui",
            default=False,
            const=True,
            nargs="?",
            help="specify this if a gui should be opened"
        )
        self.add_argument(
            "--host-flask",
            dest="host_flask",
            default=False,
            const=True,
            nargs="?",
            help="specify this if a flask server should be hosted"
        )
        self.args = self.parse_args()
        self.args_dict = vars(self.args)
        move_to_config = ["has_techrep"]
        self.args_dict["configs"] = {k: self.args_dict.pop(k) for k in move_to_config}


def browsefunc(fn, var, fn_params: dict = None):
    if fn_params is None:
        fn_params = {}
    filename = fn(**fn_params)
    var.set(filename)
