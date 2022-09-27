import logging
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from mspypeline.core import MSPInitializer
from mspypeline.core.MSPPlots import BasePlotter
from mspypeline.plotting_backend import matplotlib_plots

class MQReader:  # TODO currently circular dependency
    pass


class MaxQuantPlotter(BasePlotter):
    """
    MaxQuant Plotter is a child class of the :class:`BasePlotter` and inherits all functionality to get data and
    generate plots.
    """
    def __init__(
            self,
            start_dir: str,
            reader_data: dict,
            intensity_df_name: str = "proteinGroups",
            interesting_proteins: dict = None,
            go_analysis_gene_names: dict = None,
            configs: dict = None,
            required_reader="mqreader",
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
            intensity_entries=(("raw", "Intensity ", "Intensity"), ("lfq", "LFQ intensity ", "LFQ intensity"),
                               ("ibaq", "iBAQ ", "iBAQ intensity")),
            intensity_df_name="proteinGroups",
            required_reader="mqreader"
        )
        default_kwargs.update(**kwargs)
        return super().from_MSPInitializer(mspinit_instance, **default_kwargs)

    @classmethod
    def from_file_reader(cls, reader_instance: MQReader, **kwargs):
        default_kwargs = dict(
            intensity_df_name="proteinGroups",
            intensity_entries=(("raw", "Intensity ", "Intensity"), ("lfq", "LFQ intensity ", "LFQ intensity"),
                               ("ibaq", "iBAQ ", "iBAQ intensity")),
        )
        default_kwargs.update(**kwargs)
        return super().from_file_reader(reader_instance, **default_kwargs)

    def create_report(self, target_dir: str = None):
        """
        Creates a MaxQuantReport.pdf, which can be used as :ref:`quality control <plotters>`.

        | For overview of plots see :ref:`analysis options <plotters>`
        | For exemplary plot see :ref:`gallery <mqreport>`

        Parameters
        ----------
        target_dir
            directory where report will be written

        """
        def bar_from_counts(ax, counts, compare_counts=None, title=None, relative=False, yscale=None,
                            ylabel = None, bar_kwargs=None):
            if ylabel is not None:
                ax.set_ylabel(ylabel)
            elif relative:
                ax.set_ylabel("Relative counts")
                counts = counts / counts.sum()
            else:
                ax.set_ylabel("Counts")
            if title is not None:
                ax.set_title(title)
            if bar_kwargs is None:
                bar_kwargs = {}
            bar_container = ax.bar([x for x in range(len(counts))], counts.values, **bar_kwargs)

            if compare_counts is not None:
                if relative:
                    compare_counts = compare_counts / compare_counts.sum()
                for bar, height in zip(bar_container, compare_counts):
                    bar_x = bar.get_x()
                    bar_w = bar.get_width()
                    ax.plot((bar_x, bar_x, bar_x + bar_w, bar_x + bar_w),
                            (0, height, height, 0), color="black")

            ax.set_xticks([i for i in range(len(counts))])
            ax.set_xticklabels(counts.index)
            if yscale is not None:
                if isinstance(yscale, str):
                    ax.set_yscale(yscale)
                elif isinstance(yscale, dict):
                    ax.set_yscale(**yscale)
            return bar_container

        def hist2d_with_hist(xdata, ydata, title=None, xlabel=None, ylabel=None, max_x=145):
            fig = plt.figure(figsize=(14, 7))
            if title is not None:
                fig.suptitle(title)
            spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[2, 1], height_ratios=[1, 2])

            ax2dhist = fig.add_subplot(spec[1, 0])
            ax1dhistvert = fig.add_subplot(spec[0, 0])
            ax1dhisthor = fig.add_subplot(spec[1, 1])

            h, xedges, yedges, image = ax2dhist.hist2d(xdata, ydata,
                                                       bins=100, range=((0, max_x), (0, 2)))  # TODO find bins
            ax2dhist.set_xlabel(xlabel)
            ax2dhist.set_ylabel(ylabel)

            ax1dhistvert.hist(xdata, bins=xedges)
            ax1dhistvert.set_ylabel("Peptide Counts")

            ax1dhisthor.hist(ydata, bins=yedges, orientation="horizontal")
            ax1dhisthor.set_xlabel("Peptide Counts")

            ax1dhistvert.set_xlim(*ax2dhist.get_xlim())
            ax1dhisthor.set_ylim(*ax2dhist.get_ylim())

            fig.tight_layout(rect=[0, 0.03, 1, 0.95])

            return fig, (ax2dhist, ax1dhistvert, ax1dhisthor)

        def get_plot_data_from_hist(data, density=False, n_bins=16):
            d_min, d_max = np.nanmin(data.values), np.nanmax(data.values)
            bins = np.linspace(d_min, d_max, n_bins)

            if data.shape[0] * data.shape[1] > 1e7:
                y = np.zeros(bins.shape[0] - 1)
                for i in range(0, data.shape[0] // 5000 + 1):
                    for j in range(0, data.shape[1] // 50 + 1):
                        y_del, x = np.histogram(data.iloc[
                                                i * 5000: (i + 1) * 5000, j * 50: (j + 1) * 50].values.flatten(),
                                                bins=bins)
                        y += y_del
                if density:
                    db = np.array(np.diff(bins), float)
                    y = y / db / y.sum()

            else:
                y, x = np.histogram(data.values.flatten(), bins=bins, density=density)
            y = np.concatenate(([0], np.repeat(y, 2), [0]))
            x = np.repeat(x, 2)
            return x, y, bins

        import matplotlib.cm as cm
        cmap = cm.get_cmap("jet")

        prefix = "Intensity "
        group_iter = None
        plot_colors = {}

        self.logger.info("Reading files")

        try:
            self.logger.debug("Reading parameters")
            parameters = self.required_reader_data['parameters']
        except KeyError:
            self.logger.warning("Did not find parameters")
            parameters = None
        try:
            self.logger.debug("Reading summary")
            summary = self.required_reader_data['summary']
            summary.columns = [sum_col.upper() for sum_col in summary.columns]
        except KeyError:
            self.logger.warning("Did not find summary")
            summary = None
        try:
            self.logger.debug("Reading peptides")
            peptides = self.required_reader_data["peptides"]
            peptides_prefix_columns = [x for x in peptides.columns if x.startswith(prefix)]
            peptides_intensities = peptides[peptides_prefix_columns].replace({0: np.nan})
            peptides_intensities.columns = pd.MultiIndex.from_arrays(
                [["Grouped Intensity"] * len(peptides_intensities.columns), peptides_intensities.columns],
                names=("agg", "sample")
            )

            last_aa = pd.concat([peptides["Last amino acid"].rename(col)[peptides[col].notna()]
                                 for col in peptides.columns if col.startswith("Experiment")], axis=1)
            last_aa_counts = last_aa.apply(pd.Series.value_counts)
            last_aa_counts = last_aa_counts.fillna(0).rename(lambda x: x.replace("Experiment ", ""), axis=1)

            before_aa = pd.concat([peptides["Amino acid before"].rename(col)[peptides[col].notna()]
                                   for col in peptides.columns if col.startswith("Experiment")], axis=1)
            before_aa_counts = before_aa.apply(pd.Series.value_counts)
            before_aa_counts = before_aa_counts.fillna(0).rename(lambda x: x.replace("Experiment ", ""), axis=1)
        except (KeyError, ValueError):
            self.logger.warning("Did not find peptides")
            peptides = None
        try:
            self.logger.debug("Reading proteinGroups")
            prot_groups = self.required_reader_data["proteinGroups"]
            prot_groups_prefix_columns = [x for x in prot_groups.columns if x.startswith(prefix)]
            prot_groups_colors = [x.replace(prefix, "") for x in prot_groups_prefix_columns]
            plot_colors.update({col: cmap(i/len(prot_groups_colors)) for i, col in enumerate(prot_groups_colors)})
            prot_groups_intensities = prot_groups[prot_groups_prefix_columns].replace({0: np.nan})
            prot_groups_intensities.columns = pd.MultiIndex.from_arrays(
                [["Grouped Intensity"] * len(prot_groups_intensities.columns), prot_groups_intensities.columns],
                names=("agg", "sample")
            )
            has_lfq = str(any([x.startswith("LFQ") for x in prot_groups.columns]))
            has_ibaq = str(any([x.startswith("iBAQ") for x in prot_groups.columns]))
        except KeyError:
            self.logger.warning("Did not find proteinGroups")
            prot_groups = None
            has_lfq = "File is missing"
            has_ibaq = "File is missing"
        try:
            contaminants = self.required_reader_data["contaminants"]
        except KeyError:
            self.logger.warning("No contaminants found")
            contaminants = None
        try:
            self.logger.debug("Reading evidence")
            evidence = self.required_reader_data["evidence"]
            mz = evidence.pivot(index=None, columns="Experiment", values="m/z")
            plot_colors.update({col: cmap(i/len(mz.columns)) for i, col in enumerate(mz.columns)})
            charge = evidence.pivot(index=None, columns="Experiment", values="Charge")
            charge = charge.apply(pd.Series.value_counts)
            charge.index = charge.index.astype(int)
            missed_cleavages = evidence.pivot(index=None, columns="Experiment", values="Missed cleavages")
            missed_cleavages = missed_cleavages.apply(pd.Series.value_counts)
            missed_cleavages.index = missed_cleavages.index.astype(int)
            retention_length = evidence.pivot(index=None, columns="Experiment", values="Retention length")
            retention_time = evidence.pivot(index=None, columns="Experiment", values="Retention time")
        except KeyError:
            self.logger.warning("Did not find evidence")
            evidence = None
        try:
            self.logger.debug("Reading msScans")
            ms_scans = self.required_reader_data["msScans"]
            ms_scan_groups = ms_scans.groupby("Raw file")
            group_iter = ms_scan_groups.groups
        except KeyError:
            self.logger.warning("Did not find msScans")
            ms_scans = None
        try:
            self.logger.debug("Reading msmsScans")
            msms_scans = self.required_reader_data["msmsScans"]
            msms_scan_groups = msms_scans.groupby("Raw file")
            group_iter = msms_scan_groups.groups
        except KeyError:
            self.logger.warning("Did not find msmsScans")
            msms_scans = None

        self.logger.info("Creating plots")
        target_dir = target_dir if target_dir is not None else self.start_dir
        with PdfPages(os.path.join(target_dir, "MaxQuantReport.pdf")) as pdf:
            self.logger.debug("Creating start page")
            fig = plt.figure(figsize=(14, 7))
            text_conf = dict(transform=fig.transFigure, size=24, ha="center")
            fig.text(0.5, 0.92, "MaxQuant report", **text_conf)
            text_conf.update({"size": 20})
            fig.text(0.5, 0.85, "parameter.txt info", **text_conf)
            text_conf.pop("size")
            if parameters is not None:
                fig.text(0.5, 0.8, f"Version: {parameters.get('Version', 'missing')}"
                         f"run at: {parameters.get('Date of writing', 'missing')}", **text_conf)
                fig.text(0.5, 0.75, f"Fasta File: {os.path.split(parameters.get('Fasta file', os.path.join('missing', 'file')))[1]}, "
                         f"Match between runs: {parameters.get('Match between runs', 'missing')}", **text_conf)
                fig.text(0.5, 0.7, "Min. to Max. peptide length for unspecific search: "
                         f"{parameters.get('Min. peptide length for unspecific search', 'missing')} to "
                         f"{parameters.get('Max. peptide length for unspecific search', 'missing')}", **text_conf)
            else:
                fig.text(0.5, 0.8, "Missing", **text_conf)

            text_conf.update({"size": 20})
            fig.text(0.5, 0.65, "summary.txt info", **text_conf)
            text_conf.pop("size")
            if summary is not None:
                fig.text(0.5, 0.6, f"Used Enzyme: {summary.loc[0, 'ENZYME']}", **text_conf)
                fig.text(0.5, 0.55, f"Variable modifications: {summary.loc[0, 'VARIABLE MODIFICATIONS']}", **text_conf)
                fig.text(0.5, 0.5, f"Mass Standard Deviation: mean {summary.loc[:, 'MASS STANDARD DEVIATION [PPM]'].mean():.5f} ppm, max {summary.loc[:, 'MASS STANDARD DEVIATION [PPM]'].max():.5f} ppm", **text_conf)
            else:
                fig.text(0.5, 0.6, "Missing", **text_conf)

            if prot_groups is not None:
                fig.text(0.5, 0.45, f"Identified proteins (without contaminants): {prot_groups.shape[0]}", **text_conf)
            if peptides is not None:
                fig.text(0.5, 0.4, f"Identified peptides (without contaminants): {peptides.shape[0]}", **text_conf)
            fig.text(0.5, 0.35, f"Has LFQ intensities: {has_lfq}", **text_conf)
            fig.text(0.5, 0.3, f"Has iBAQ: {has_ibaq}", **text_conf)

            pdf.savefig()
            plt.close(fig)
            # ######

            # figure
            if peptides is not None:
                self.logger.debug("Creating peptide overview")
                fig, axarr = plt.subplots(3, 1, figsize=(14, 7))
                bar_from_counts(axarr[0], peptides["Missed cleavages"].value_counts(), title="Missed Cleavages", relative=True)
                bar_from_counts(axarr[1], peptides["Amino acid before"].value_counts(), title="Amino acid before", yscale="log")
                bar_from_counts(axarr[2], peptides["Last amino acid"].value_counts(), title="Last amino acid", yscale="log")

                fig.tight_layout(rect=[0, 0.03, 1, 0.95])

                pdf.savefig()
                plt.close(fig)
            # ######

            # figure stuff
            self.logger.debug("Creating technical overview")
            fig, axarr = plt.subplots(2, 1, figsize=(14, 7))
            if peptides is not None:
                bar_from_counts(axarr[0], peptides["Charges"].str.split(";").explode().value_counts().sort_index(),
                                title="Peptide Charges")

            if evidence is not None:
                axarr[1].hist(evidence["m/z"])
                axarr[1].set_xlabel("m/z")
                axarr[1].set_ylabel("Counts")
                axarr[1].set_title("Peptide m/z")

            fig.tight_layout(rect=[0, 0.3, 1, 0.95])

            pdf.savefig()
            plt.close(fig)
            # ###########

            # hist with peptide m/z from evidence["m.s"]
            self.logger.debug("Creating identified proteins and peptides per sample")
            fig, axarr = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
            # hist with identified proteins and hist with identified peptides, shared axis
            if prot_groups is not None:
                identified_proteins = (prot_groups_intensities["Grouped Intensity"] > 0).sum()
                identified_proteins = identified_proteins.rename(lambda x: x.replace("Intensity ", "").replace("_", " ")
                                                                 , axis=0)
                bar_from_counts(axarr[0], identified_proteins, title="Identified proteins")
            # proteins from proteinGroups, peptides from peptides file per sample
            if peptides is not None:
                identified_peptides = (peptides_intensities["Grouped Intensity"] > 0).sum()
                identified_peptides = identified_peptides.rename(lambda x: x.replace("Intensity ", "").replace("_", " ")
                                                                 , axis=0)
                bar_from_counts(axarr[1], identified_peptides, title="Identified peptides")
                axarr[1].xaxis.set_tick_params(rotation=90)

            fig.tight_layout(rect=[0, 0.03, 1, 0.95])

            pdf.savefig()
            plt.close(fig)
            # #####################

            # Page with stuff
            if summary is not None:
                self.logger.debug("Creating scan overview")
                fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(14, 7))

                axarr[0].set_title("MS scans")
                axarr[0].bar(range(summary.shape[0]), summary["MS"])
                axarr[0].set_ylabel("Count")

                axarr[1].set_title("MS/MS scans")
                axarr[1].bar(range(summary.shape[0]), summary["MS/MS"])
                axarr[1].set_ylabel("Count")

                axarr[2].set_title("MS/MS identified [%]")
                axarr[2].bar(range(summary.shape[0]), summary["MS/MS IDENTIFIED [%]"])
                axarr[2].set_ylabel("Percent")
                axarr[2].set_xticks(range(summary.shape[0]))
                labels = [sample.replace("_", " ") for sample in summary["EXPERIMENT"]]
                axarr[2].set_xticklabels(labels, rotation=90)

                fig.tight_layout(rect=[0, 0.03, 1, 0.95])

                pdf.savefig()
                plt.close(fig)
            # ##################

            # page with stuff
            if contaminants is not None:
                self.logger.debug("Creating overview of share of contamination intensity from total")
                df_contaminants_int = contaminants[prot_groups_prefix_columns]
                df_no_contaminants_int = prot_groups[prot_groups_prefix_columns]
                sum_int = pd.DataFrame({
                    "contaminants": df_contaminants_int.sum(axis=0),
                    "no_contaminants": df_no_contaminants_int.sum(axis=0)})
                sum_int["total"] = sum_int["contaminants"] + sum_int["no_contaminants"]
                sum_int["percent_cont"] = (sum_int["contaminants"] / sum_int["total"])*100
                sum_int["percent_no_cont"] = (sum_int["no_contaminants"] / sum_int["total"])*100
                fig, ax = plt.subplots(1, 1, sharex=True, figsize=(14, 7))

                labels = [label.replace("Intensity ", "").replace("_", " ") for label in sum_int.index.values]

                ax.bar(x=labels, height=sum_int["percent_cont"], width=0.8)
                ax.set_title("Intensity of contaminants from total intensity", size=14, pad=20)
                ax.set_ylabel("Percent", size=13)
                ax.set_xticks(range(len(labels)))
                ax.set_xticklabels(labels, rotation=90)
                ax.axhline(5, linestyle="-", linewidth=2, color="red", alpha=0.6)

                fig.tight_layout()

                pdf.savefig(figure=fig)
                plt.close(fig)
            # ############

            # page with stuff
            if prot_groups is not None:
                self.logger.debug("Creating overall intensity histograms")
                fig, axarr = plt.subplots(2, 1, sharex=True, figsize=(7, 7))

                # stacked histogram of log2 intensities
                colors = prot_groups_intensities["Grouped Intensity"].rename(lambda x: x.replace("Intensity ", ""), axis=1).columns
                colors = [plot_colors[c] for c in colors]
                matplotlib_plots.save_intensity_histogram_results(prot_groups_intensities, n_bins=11, histtype="barstacked",
                                                                  plot=(fig, axarr[0]), color=colors, show_mean=False,
                                                                  legend=False)
                # overlayed histogram of log2 intensities
                matplotlib_plots.save_intensity_histogram_results(prot_groups_intensities, n_bins=11, histtype="step",
                                                                  plot=(fig, axarr[1]), color=colors, show_mean=False,
                                                                  legend=False)
                fig.legend(bbox_to_anchor=(1.02, 0.5), loc="center left")

                fig.tight_layout()

                pdf.savefig()
                plt.close(fig)
            # ############

            # page with stuff
            # two histograms with heatmap
            # retention time vs retention length
            # from evidence["Retention time"], evidence["Retention length"]
            if evidence is not None:
                self.logger.debug("Creating overall retention time vs retention length")

                fig, ax = hist2d_with_hist(title="Overall Retention time vs Retention length",
                                           xdata=evidence["Retention time"], ydata=evidence["Retention length"],
                                           xlabel="Retention time [min]", ylabel="Retention length [min]",
                                           max_x=evidence["Retention time"].max())

                pdf.savefig(figure=fig)
                plt.close(fig)
            # ##############

            # individual comparison
            if evidence is not None and peptides is not None:
                self.logger.debug("Creating individual experiment comparison")
                charge_flat = charge.sum(axis=1)
                missed_cleavages_flat = missed_cleavages.sum(axis=1)
                before_aa_counts_flat = before_aa_counts.sum(axis=1)
                last_aa_counts_flat = last_aa_counts.sum(axis=1)

                mz_x, mz_y, mz_bins = get_plot_data_from_hist(mz, n_bins=15, density=True)

                for experiment in mz.columns:
                    plot_color = plot_colors[experiment]
                    fig, axarr = plt.subplots(3, 2, figsize=(14, 7))
                    fig.suptitle(experiment.replace("_", " "))

                    axarr[0, 0].hist(mz[experiment], density=True, color=plot_color, bins=mz_bins)
                    axarr[0, 0].plot(mz_x, mz_y, color="black")
                    # axarr[0, 0].hist(mz.drop(experiment, axis=1).values.flatten(), histtype="step", density=True, color="black", bins=bins, linewidth=2)
                    # axarr[0, 0].hist(mz_flat, histtype="step", density=True, color="black", bins=bins, linewidth=2)
                    axarr[0, 0].set_xlabel("m/z")
                    axarr[0, 0].set_ylabel("Density")
                    axarr[0, 0].set_title("Peptide m/z")

                    bar_from_counts(axarr[0, 1], charge[experiment],
                                    compare_counts=charge_flat,
                                    relative=True,
                                    title="Peptide charges", bar_kwargs={"color": plot_color})
                    axarr[0, 1].set_xlabel("Peptide charge")

                    bar_from_counts(axarr[1, 0], missed_cleavages[experiment],
                                    compare_counts=missed_cleavages_flat,
                                    relative=True,
                                    title="Number of missed cleavages", bar_kwargs={"color": plot_color})
                    axarr[1, 0].set_xlabel("Missed cleavages")

                    # TODO this might be missing
                    bar_from_counts(axarr[1, 1], before_aa_counts[experiment],
                                    compare_counts=before_aa_counts_flat,
                                    relative=True,
                                    bar_kwargs={"color": plot_color})
                    axarr[1, 1].set_title("Amino acid before")

                    bar_from_counts(axarr[2, 0], last_aa_counts[experiment],
                                    compare_counts=last_aa_counts_flat,
                                    relative=True,
                                    bar_kwargs={"color": plot_color})
                    axarr[2, 0].set_title("Last amino acid")

                    axarr[2, 1].remove()
                    fig.tight_layout()

                    pdf.savefig()
                    plt.close(fig)
            # ###############

            # Intensity histograms of individual samples compared to remaining
            if prot_groups is not None:
                self.logger.debug("Creating individual intensity histograms")
                log2_intensities = np.log2(prot_groups_intensities["Grouped Intensity"])
                log2_intensities = log2_intensities.rename(lambda x: x.replace("Intensity ", ""), axis=1)

                b, h, bins = get_plot_data_from_hist(log2_intensities, density=True, n_bins=16)

                n_figures = int(np.ceil(len(log2_intensities.columns) / 9))
                for n_figure in range(n_figures):
                    fig, axarr = plt.subplots(3, 3, figsize=(15, 15))
                    for i, (pos, ax) in enumerate(np.ndenumerate(axarr)):
                        idx = n_figure * 9 + i
                        try:
                            experiment = log2_intensities.columns[idx]
                        except IndexError:
                            break
                        ax.hist(log2_intensities.loc[:, experiment], bins=bins, density=True,
                                color=plot_colors[experiment])
                        ax.plot(b, h, color="black")
                        ax.set_title(experiment.replace("_", " "))
                        ax.set_xlabel("Intensity")
                        ax.set_ylabel("Density")

                    if n_figure == (n_figures - 1):
                        n_empty = n_figures * 9 - len(log2_intensities.columns)
                        for i in range(1, n_empty + 1):
                            axarr.flat[-i].remove()

                    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

                    pdf.savefig(fig)
                    plt.close(fig)
            # ################

            # Retention time of individuals samples vs remaining
            if evidence is not None:
                self.logger.debug("Creating individual retention time histograms")
                b, h, bins = get_plot_data_from_hist(retention_time, density=True, n_bins=25)

                n_figures = int(np.ceil(len(retention_time.columns) / 9))
                for n_figure in range(n_figures):
                    fig, axarr = plt.subplots(3, 3, figsize=(15, 15))
                    for i, (pos, ax) in enumerate(np.ndenumerate(axarr)):
                        idx = n_figure * 9 + i
                        try:
                            experiment = retention_time.columns[idx]
                        except IndexError:
                            break
                        ax.hist(retention_time.loc[:, experiment], bins=bins, density=True,
                                color=plot_colors[experiment])
                        ax.plot(b, h, color="black")
                        ax.set_title(experiment.replace("_", " "))
                        ax.set_xlabel("Retention time")
                        ax.set_ylabel("Density")

                    if n_figure == (n_figures - 1):
                        n_empty = n_figures * 9 - len(retention_time.columns)
                        for i in range(1, n_empty + 1):
                            axarr.flat[-i].remove()

                    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

                    pdf.savefig(fig)
                    plt.close(fig)

            # retention time vs retention length individual
            if evidence is not None:
                self.logger.debug("Creating individual retention time vs retention length")
                for experiment in retention_length.columns:
                    fig, ax = hist2d_with_hist(title=experiment.replace("_", " "), xdata=retention_time[experiment],
                                               ydata=retention_length[experiment], xlabel="Retention time [min]",
                                               ylabel="Retention length [min]", max_x=retention_time.max().max())

                    pdf.savefig(figure=fig)
                    plt.close(fig)
            # ##########

            # total ion current vs retention length
            import matplotlib.ticker as ticker

            @ticker.FuncFormatter
            def scientific_formatter(x, pos):
                if x != 0:
                    return f"{x:.1E}"
                else:
                    return "0"

            if group_iter is not None:
                self.logger.debug("Creating MS scan and MSMS scan overview")
                for n_plot in range(int(np.ceil(len(group_iter) / 4))):
                    fig = plt.figure(figsize=(14, 7))
                    outer = fig.add_gridspec(2, 2, wspace=0.2, hspace=0.4)

                    for i in range(4):
                        inner = outer[i].subgridspec(2, 1, wspace=0.1, hspace=0.0)

                        group_counter = 4 * n_plot + i
                        try:
                            if msms_scans is not None:
                                group_name = list(msms_scan_groups.groups.keys())[group_counter]
                            elif ms_scans is not None:
                                group_name = list(ms_scan_groups.groups.keys())[group_counter]
                            else:
                                raise ValueError("Logic error")
                        except IndexError:
                            break

                        # msms plot
                        ax_msms: plt.Axes = plt.subplot(inner[1])
                        ax_msms.text(0.1, 0.9, 'MSMS', horizontalalignment='center',
                                     verticalalignment='center', transform=ax_msms.transAxes)
                        if msms_scans is not None:
                            df = msms_scan_groups.get_group(group_name)
                            ax_msms.plot(df["Retention time"], df["Total ion current"], color="black", linewidth=0.2)
                        ax_msms.yaxis.set_major_formatter(scientific_formatter)
                        ax_msms.set_xlabel("Retention time")
                        ax_msms.set_ylabel("Total ion current")
                        fig.add_subplot(ax_msms)

                        # ms plot
                        # get the axis with shared x axis
                        ax_ms: plt.Axes = plt.subplot(inner[0], sharex=ax_msms)
                        # add the text
                        ax_ms.text(0.1, 0.9, 'MS', horizontalalignment='center',
                                   verticalalignment='center', transform=ax_ms.transAxes)
                        # disable the axis ticks
                        ax_ms.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
                        ax_ms.set_title(group_name)
                        if ms_scans is not None:
                            df = ms_scan_groups.get_group(group_name)
                            ax_ms.plot(df["Retention time"], df["Total ion current"], color="black", linewidth=0.2)
                        ax_ms.yaxis.set_major_formatter(scientific_formatter)
                        fig.add_subplot(ax_ms)

                    pdf.savefig()
                    plt.close(fig)

            self.logger.info("Done creating report")
