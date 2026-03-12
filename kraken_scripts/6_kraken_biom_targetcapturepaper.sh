#!/bin/bash

# This script creates a standard .biom format for community analysis
# This version is designed to combine results from multiple sequencing batches, which were classified separately

#call kraken biom
kraken-biom ./intbio_targetcapturepaper_report/*.k2_bracken_classes.report --fmt json --max C --min G -o bracken_classes.targetcapturepaper.biom

#call kraken biom
kraken-biom ./intbio_targetcapturepaper_report/*.k2_bracken_orders.report --fmt json --max C --min G -o bracken_orders.targetcapturepaper.biom

#call kraken biom
kraken-biom ./intbio_targetcapturepaper_report/*.k2_bracken_families.report --fmt json --max C --min G -o bracken_families.targetcapturepaper.biom

#call kraken biom
kraken-biom ./intbio_targetcapturepaper_report/*.k2_bracken_genuses.report --fmt json --max C --min G -o bracken_genera.targetcapturepaper.biom
