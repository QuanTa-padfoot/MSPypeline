Custom Gene Lists
=====================

Pathway and GO analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

| Pathway and GO analysis gene set files are lists of genes that are somewhat associated with the respective pathway
  or GO term. Upon data preparation by the :class:`MQReader`, measured proteins are indexed by their **gene name**
  originating from the FASTA file (header) to which the protein was mapped. Thus, measured protein intensities can be
  analysed using functional gene sets like those that are incorporated by the ``mspypeline`` package.
| Since some analysis methods deploy such pathway and GO term lists, the provided files can be used to get an idea and
  generate an example for a potential analysis (see :ref:`gallery`).

  * GO lists are required for the GO analysis plot using :meth:`~mspypeline.BasePlotter.plot_go_analysis`, in
    which an enrichment analysis for each selected GO Term file is created.
  * Pathway lists serve multiple plotting options including the rank plot using
    :meth:`~mspypeline.BasePlotter.plot_rank`, the pathway analysis plot using
    :meth:`~mspypeline.BasePlotter.plot_pathway_analysis` and the volcano plot using
    :meth:`~mspypeline.BasePlotter.plot_r_volcano`.

| Pathway and GO term lists are configured globally for a data analysis, as a single selection of gene lists, and not
  individually for each plot. If one or more gene lists are selected in the GUI or :ref:`configs file <default-yaml>`,
  they will be used for any of the plots listed above that make use of them, whether that is a GO enrichment plot or
  one of the pathway-based plots.
| To change the choice of gene lists for an analysis, the previously selected lists have to be un-checked
  in the gene list selection box in the GUI or the *"gene_lists:"* argument in the configs file has
  to be edited manually.

.. tip::
    Any desired pathway and GO list can be manually provided to ``mspypeline`` by the user. The file simply has to:

    1. follow the *one-column-txt-format* that can be seen in the exemplary files listed below,
    2. be stored in one of these two locations:

        - saved in the *.../mspypeline/config/gene_lists* directory, where all the other files are stored (files
          saved here are available for all experiments and from the GUI).
        - saved in a *gene_lists* directory in the same location where the txt folder of the experiment
          data is stored (files saved here are available for the particular experiment and are callable when
          ``mspypeline`` is used as a :ref:`python module <python-quickstart>` or when the list is specified in the
          :ref:`configs file <default-yaml>`.

.. attention::
    * All pathway and GO analysis gene files provided here are based on the **HUMAN** genome/proteome.
    * All pathway and GO analysis gene files are retrieved from the open source
      `GSEA Molecular Signature Data Base <https://www.gsea-msigdb.org/gsea/msigdb/index.jsp>`__ (22. Feb. 2021)

.. _pathway-proteins:
.. _go-term-proteins:

Gene Lists (Pathways & GO Terms)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. toctree::
   :glob:

   gene_lists/*