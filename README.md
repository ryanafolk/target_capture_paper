# Triple domain target capture
This repository contains code, data, and figures for the paper entitled, "Multidomain triple target capture for navigating complex symbioses." 

## Base level repository
`kraken2_analysis_targetcapturepaper_16S.R` is the community analysis script for bacteria. `kraken2_analysis_targetcapturepaper_ITS.R` is the community analysis script for fungi. The base level contains scripts and some miscellaneous data. The mock community files (`mock_community.csv`, `mock_community_ITS.csv`, and `mock_community_amplicon.csv` contain relative abundance data spreadsheets representing the mock community labeled organism abundance.

## Biom files
Two directories, `biom_files_amplicon` and `biom_files_targetcapture`, contain community data in the biom format.

## Figures
Six directories, `figs_ampliconsequencing`, `figs_capture_success`, `figs_combined`, `figs_mockcommunity_ampliconsequencing/16S`, `figs_mockcommunity_targetcapture`, and `figs_targetcapture`, contain raw R figure output, which were assembled into multipanel figures for the paper. 

## metadata
This directory contains metadata files used for aggregating samples and reporting QC. `NEON-SiteMap-Table.csv` relates the location of NEON sites for the samples. `amplicon_original_metadata.tsv` describes sample metadata, such as plant host and site, for the amplicon sequencing data. `amplicon_resequencing_metadata.tsv` contains similar data for amplicon samples that were resequenced using target capture. `capture_qc.csv` contains QC information for target capture across the nearly 3,000 NEON project samples. `intbio_metadata_draft.csv` similarly contains metadata information for the NEON samples. `capture_paper_functionalcomparison.tsv` lists samples used for the functional target capture gene recovery figure.

## kraken_scripts
This directory contains six shell scripts representing the Kraken and Braken scripts used to generate biom files. These are numbered and should be run in order.

## functional_annotation
This directory contains scripts for the functional locus recovery analysis.

## baits
Bait design files. `Combined-masked-clust100_used.fas` contains the baits that were synthesized for ITS and symbiosis islands. `centroids_largesubunit.0.9id_0.9cov.fasta`, `centroids_5.8S.0.9id_0.9cov.fasta`, and `centroids_smallsubunit.0.9id_0.9cov.fasta` are the references from which ITS probes were synthesized, according to the strategy describe in the paper; it is recommended to use UNITE as an assembly reference and not this file. `symbiosis_islands.fasta` is the reference from which symbiosis island probes were synthesized, and may be used as an assembly reference.

Note on sample names: NEON samples are of the form: 4letterSiteCode-plantSerialNumber-replicateNumber-sampleType. Sample codes are described in `./metadata/NEON-SiteMap-Table.csv`. Sample types are described in the text, with Rh = "rhizosphere," Ro = "root," and No = "nodule."
