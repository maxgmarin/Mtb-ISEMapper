# Mtb-ISEMapper

The `MTb-ISEMapper` pipeline analyzes paired-end short-read WGS data for M. tuberculosis isolates.

At a high-level the pipeline has the following steps:
1) Subset the input data for reads that align to the IS6110 reference sequence. (IS6110 reads)
2) Align the IS6110 reads to the H37Rv genome. Note this version of the H37Rv genome must have all IS6110 elements masked.
3) Calculate alignment depth for IS6110 reads across the H37Rv genome.
4) Analyze depth to identify putative IS6110 insertion locations.

## Installation:

<Installation instructions will go here>

## Usage:
<Usage instructions will go here>

## FAQ:

<FAQ will go here>
