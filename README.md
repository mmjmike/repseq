# RepSeq

## Purpose

The purpose of this Python library is to make the postanalysis of TCR/BCR repertoire sequencing data easy and modular.

## Main modules

### Clonosets

Manipulations with groups of clonosets:
- creating clonosets DataFrames by searching for typical filenames in directories
- detecting clonoset formats
- filtering out non-target clonosets (same sample_id, but with low % of reads)
- pooling clonosets into one

### Stats

- Basic clonoset properties, like clone/read/umi counts functional or with OOF/Stops
- CDR3 amino acid properties: N-counts, physico-chemical properties, Kidera Factors
- Diversity statistics: observed diversity, (normalized) Shannon-Wiener, chao1
- Convergence estimate
- V/D/J/C-gene frequencies or VJ-combinations

All calculations are parallelized

### Intersections

- Find parwise intersecting clonotypes between clonosets
- Intersection metrics: F, F2, D
- Count tables for clonotypes (similarity groups of clonotypes)
- Intersect clusters with clonosets
- TCRnet integration

### Clustering

This module implements different immune repertoire clustering analyses:
- basic clustering: using hamming distance similarity in CDR3 regions and/or same V/J-segments
- tcr-dist clustering: using distance metrics from tcrdist software
- ALICE: find expanded clonotypes by analyzing the probability of neighbour generation with OLGA algrorithm
- split clusters with community detection algorithms (Louvain, Leiden)
- easy and modular customisation of cluster analysis
- output graphs and metadata to Cytoscape format

Graph representation with NetworkX library

### DiffExp

- find differentially expressing clonotypes/clusters of clonotypes in CFSE-assays or similar experiments

### Clone Filter

Easy filtering of clonosets by one Filter object, integrated with other analysis procedures

Filtering includes following features:
- counting by reads/UMIs/clonotypes
- use top N clonotypes (tails mixing included for randomly mixing clonotypes with same counts)
- randomly downsample to N UMIs/reads (you can specify `seed`, highly recommended for reproducibility)
- remove low count clonotypes
- filter out non-functional(OOF,Stop in CDR3)/functional clonotypes
- white/black list of clonotypes
- recount frequencies (by reads/UMIs)
- convert to vdjtools-like format
- combine (pool) clonotypes with similar features: CDR3/V/J 

### IO-module

- reads and understands clonosets of following formats: MiXCR 3/4, vdjtools, Adaptive Biosciences
- tsv, .gz, .zip

### MiXCR

As MiXCR is the leading software for generating clonoset tables from raw FastQ files this module helps to run MiXCR 4.3+ batch analises with SLURM queue manager.

Easy accumulation of most sensible processing data from json-reports of MiXCR into one table