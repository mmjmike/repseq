# IO module

This module contains functions for input-output procedures

## Clonoset

Lightweight object wrapper around a pandas clonoset table.

::: clonoset.Clonoset
    options:
          show_root_toc_entry: false

## standardize_to_vdjtools_columns

::: clonoset.standardize_to_vdjtools_columns
    options:
          show_root_toc_entry: false

## read_ngsik_metadata

::: io.read_ngsik_metadata
    options:
          show_root_toc_entry: false

## read_yaml_metadata

Deprecated alias for `read_ngsik_metadata`.

::: io.read_yaml_metadata
    options:
          show_root_toc_entry: false

## load_olga_models

Load OLGA probability-generation and sequence-generation models from an
organism model folder. OLGA is an optional dependency and can be installed
with `pip install repseq[pgen]`.

```py
from repseq.io import load_olga_models

human_olga_model_folder = "/home/mmyshkin/soft/OLGA/olga/default_models/human_T_beta/"
hum_pgen_model, hum_seq_gen_model = load_olga_models(human_olga_model_folder)
```

By default, the function loads `model_params.txt`, `model_marginals.txt`,
`V_gene_CDR3_anchors.csv`, and `J_gene_CDR3_anchors.csv`. Each filename can be
changed independently:

```py
hum_pgen_model, hum_seq_gen_model = load_olga_models(
    human_olga_model_folder,
    params_filename="custom_model_params.txt",
    marginals_filename="custom_model_marginals.txt",
    v_anchor_filename="custom_V_gene_CDR3_anchors.csv",
    j_anchor_filename="custom_J_gene_CDR3_anchors.csv",
)
```

::: io.load_olga_models
    options:
          show_root_toc_entry: false

## read_clonoset

::: io.read_clonoset
    options:
          show_root_toc_entry: false

## save_to_vdjtools

::: io.save_to_vdjtools
    options:
          show_root_toc_entry: false

## save_to_airr

::: io.save_to_airr
    options:
          show_root_toc_entry: false

## read_json_report

::: io.read_json_report
    options:
          show_root_toc_entry: false
