# Modules

## [Clonosets](functions.md#clonosets)

Manipulations with groups of clonosets:

* Creates clonosets DataFrames by searching for typical filenames in directories
* Detects clonoset formats
* Filters out non-target clonosets (same sample_id but low % of reads)
* Pools clonosets into one

## [Stats](functions.md#stats)

* Basic clonoset properties, like clone/read/umi counts functional or with OOF/Stops
* CDR3 amino acid properties: N-counts, physico-chemical properties, Kidera Factors
* Diversity statistics: observed diversity, (normalized) Shannon-Wiener, chao1
* Convergence estimate
* V/D/J/C-gene frequencies or VJ-combinations
* All calculations are parallelized

## [Intersections](functions.md#intersections)

* Finds parwise intersecting clonotypes between clonosets
* Intersection metrics: F, F2, D
* Count tables for clonotypes (similarity groups of clonotypes)
* Intersect clusters with clonosets
* TCRnet integration

## [Clustering](functions.md#clustering)

This module implements different immune repertoire clustering analyses:

* basic clustering: using hamming distance similarity in CDR3 regions and/or same V/J-segments
* tcr-dist clustering: using distance metrics from tcrdist software
* ALICE: find expanded clonotypes by analyzing the probability of neighbour generation with OLGA algrorithm
* split clusters with community detection algorithms (Louvain, Leiden)
* easy and modular customisation of cluster analysis
* output graphs and metadata to Cytoscape format

<br>Graph representation with NetworkX library.

## [Diffexp](functions.md#diffexp)

Finds differentially expressing clonotypes/clusters of clonotypes in CFSE-assays or similar experiments.

## [Clone Filter](functions.md#clone_filter)

Easy filtering of clonosets by one Filter object, integrated with other analysis procedures.
<br>Filtering includes following features:

* counting by reads/UMIs/clonotypes
* use top N clonotypes (tails mixing included for randomly mixing clonotypes with same counts)
* randomly downsample to N UMIs/reads (you can specify seed, highly recommended for reproducibility)
* remove low count clonotypes
* filter out non-functional(OOF,Stop in CDR3)/functional clonotypes
* white/black list of clonotypes
* recount frequencies (by reads/UMIs)
* convert to vdjtools-like format
* combine (pool) clonotypes with similar features: CDR3/V/J

## [IO module](functions.md#io)

* reads and understands clonosets of following formats: MiXCR 3/4, vdjtools, Adaptive Biosciences
* tsv, .gz, .zip

## [MiXCR module](functions.md#mixcr)

As MiXCR is the leading software for generating clonoset tables from raw FastQ files this module helps to run MiXCR 4.3+ batch analyses with SLURM queue manager.
<br>Easy accumulation of most sensible processing data from json-reports of MiXCR into one table.




