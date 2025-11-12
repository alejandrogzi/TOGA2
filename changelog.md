## v2.0.6
* New TOGA2 mode added: `integrate` (early access functionality)
* Query gene inference improved:
    * Gene inference step is moved from `loss_summary` to a separate workflow step now placed before gene loss summary step;
    * **All** query genes now get their names after progenitor genes in the reference, with the following prefixes:
        * intact orthologous loci do not get any prefix;
        * inactivated orthologous loci get `missing_` prefix if they were inferred from `Missing` projections alone; otherwise, the prefix is `lost_`;
        * paralogous loci get the `paralog_` prefix;
        * retrogene loci (annotated based on `Intact/Fully Intact` processed pseudogene projections) get the `retro_` prefix
* Updated loss summary files:
    * `loss_summary` step now goes after `gene_inference` to accommodate for transcript/gene status changes 
    after gene inference step;
    * top-level file `loss_summary.tsv` now contains loss statuses only for projections appearing in the final output files (`query_annotation.bed`, `query_annotation.with_utrs.bed`, UCSC browser files). It still contains loss statuses for **all** reference transcripts and genes. For loss statuses for all anotated projections, including rejected items see `meta/loss_summary_extended.tsv`
* '#paralog' postfix added for paralogous projections in the final output files (`query_annotation.bed`, `query_annotation.with_utrs.bed`, UCSC browser files)
* '#paralog' and "#retro" suffixes also added to respective projections' names in `nucleotide.fa.gz`, `protein.fa.gz`, `loss_summary.tsv`
* Entries for fragmented projections get numerical postfixes based on their order in the restored query sequence. Numbers are 1-based and are separated from the base name with dollar sign ('$')
* A GTF format copy is produced for final query annotation file (either `query_annotation.bed` or `query_annotation.with_utrs.bed`)
* Folder `nextflow_configs` contains example config files and executor script for parallel steps performed with Nextflow
* Bug fixes:
    * Third-party binaries are now sought in the `bin/` directory first.
    * Removed hardcoded instances of `bedToBigBed` and `ixIxx` in `src/python/modules/make_ucsc_report.py`
    * Error-exit if all batches for a given step failed prior to ok-file check
    * Projections discarded at gene tree filtering step now removed from the final output files
    * Added post-gene-tree orthology resolution step to the main logging channel
    * 'missing_' query gene inference;
    * Parially Intact consistently removed from accepted retrogene/trusted second-level ortholog statuses
