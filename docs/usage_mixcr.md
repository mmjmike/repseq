# Usage: working with MiXCR

MiXCR is the leading software for generating clonoset tables from raw FastQ files. [MiXCR module](functions.md#mixcr) allows to run MiXCR 4.3+ batch analyses with SLURM queue manager.

!!! note "Setting up the environment"
    Before getting started, make sure that main_repseq environment is chosen. Otherwise, check the [installation guide](installation.md)

``` py
from repseq import mixcr as mx
from repseq import slurm
from repseq import io as repseqio

# 1) Create sample_df from dataset metadata in .yaml format (if it's in a tabular format, use external libraries such as Pandas). 
# Remove unnesessary columns if needed.
# If METADATA_FILENAME is absent, it is set to `metadata.yaml` by default. If your dataset does not have metadata, create the dataframe manually. 
# The neccessary columns are: R1, R2, sample_id, where R1 and R2 contain paths (using full paths is advised) to respective raw files, 
# and sample_id are arbitrary unique identificators
sample_df = repseqio.repseqio.read_yaml_metadata(RAW_DATA_DIR, METADATA_FILENAME)
metadata = sample_df.prop(columns=['R1', 'R2'])
output_dir = ...
path_to_mixcr_binary = ...

# 2) Create a command template for mixcr analyze. The default template is `mixcr analyze milab-human-rna-tcr-umi-multiplex -f r1 r2 output_prefix`. 
# In case of Human 7GENES DNA Multiplex MiXCR built-in preset, no template is needed.
mixcr_race_command_template = "mixcr analyze milab-human-rna-tcr-umi-race -f r1 r2 output_prefix"

# 3) Run mixcr analyze in batches
mx.mixcr4_analyze_batch(sample_df, output_dir, command_template=mixcr_race_command_template,
                        path_to_mixcr_binary, memory=32, time_estimate=1.5)
# OR 
# mx.mixcr_7genes_run_batch(sample_df, output_dir, path_to_mixcr_binary, memory=32, time_estimate=1.5)
# to check the progress, use
slurm.check_slurm_progress(os.path.join(output_dir, "mixcr_analyze_slurm_batch.log"), loop=True)

# 4) make reports (combines mixcr exportQc align, chainUsage and tags) and get report images
mx.mixcr4_reports(output_dir, mixcr_path=path_to_mixcr_binary)
slurm.check_slurm_progress(os.path.join(output_dir, "mixcr_reports_slurm_batch.log"), loop=True
mx.show_report_images(output_dir)

# 5) get a tabular report
proc_table = mx.get_processing_table(output_dir).merge(metadata)
```

??? info "Visualization"
    Properties from proc_table can be visualized in Jupyter notebook using %%R cell magic. 
    ![proc_table]("../images_docs/proc_table.png")
    ```py
    %load_ext rpy2.ipython
    %%R -i proc_table -w 900 -h 500

    params_order <- c("reads_with_umi_pc", "reads_aligned_pc","reads_per_umi","umi_in_func_clones","clones_func")

    proc_table %>%
        select(sample_id, experimental_group, subset, reads_per_umi, reads_with_umi_pc, reads_aligned_pc, clones_func, umi_in_func_clones) %>%
        pivot_longer(-c(sample_id, experimental_group, subset), names_to="parameter", values_to="value") %>%
        mutate(experimental_group=factor(experimental_group, group_order)) %>%
        mutate(parameter=factor(parameter, params_order)) %>%
        ggplot(aes(x=experimental_group, y=value, color=experimental_group)) +
            geom_boxplot(outlier.shape=NA)+
            geom_jitter()+
            facet_wrap(vars(parameter), scales="free_y")+
            scale_color_manual(values=colors_6_groups) + 
            boxplot_theme+
            theme(legend.position="none")
    ```
