# Examples

### MiXCR analyze
``` py
# Will be made executable with markdown-exec (although there should be other options)
# Not ready yet :(
import os
import pandas as pd
from IPython.display import Image, display
import seaborn as sns
import matplotlib.pyplot as plt
```

``` py
REPSEQ_PATH = '/home/epepeliaeva/soft/repseq/'

import sys

sys.path.append(REPSEQ_PATH)
from repseq import io as repseqio
from repseq import mixcr as mx
from repseq import slurm
from repseq import clonosets as cl
from repseq import stats
from repseq import clone_filter as clf
from repseq import intersections
from repseq import clustering
from repseq import logo
from repseq import vdjtools
```

``` py
MIXCR = "/projects/cdr3_software/bin/mixcr"

WORKING_DIR = "/projects/cdr3_common/repseq_demo/"
MIXCR_DIR = os.path.join(WORKING_DIR, "mixcr")

RAW_DATA_DIR = "/projects/cdr3_ngs/2023/11_room555_MiSeq_13112023/"

SAMPLE_LIST_FILENAME = os.path.join(WORKING_DIR, "sample_table.csv")
TABLE_REPORT_FILENAME = os.path.join(WORKING_DIR, "table_report.csv")

os.makedirs(WORKING_DIR, exist_ok=True)
```
``` py
WORKING_DIR = "/projects/cdr3_common/repseq_demo/"
MIXCR_DIR = os.path.join(WORKING_DIR, "mixcr")

RAW_DATA_DIR = "/projects/cdr3_ngs/2023/11_room555_MiSeq_13112023/"

SAMPLE_LIST_FILENAME = os.path.join(WORKING_DIR, "sample_table.csv")
TABLE_REPORT_FILENAME = os.path.join(WORKING_DIR, "table_report.csv")

sample_df = repseqio.read_yaml_metadata(RAW_DATA_DIR)[["sample_id", "R1", "R2"]].query('sample_id.str.contains("Rev05")')

mx.mixcr4_analyze_batch(sample_df=sample_df, 
                        output_folder = MIXCR_DIR, 
                        command_template=None,
                        mixcr_path=MIXCR, 
                        memory=32, 
                        time_estimate=1.5)

```

``` py
slurm.check_slurm_progress(os.path.join(MIXCR_DIR, "mixcr_analyze_slurm_batch.log"), loop=True)
```

```
mx.show_report_images(MIXCR_DIR)
```

``` py
proc_table = mx.get_processing_table(MIXCR_DIR)
proc_table.to_csv(TABLE_REPORT_FILENAME, index=False)
print(f"Report table saved to: {TABLE_REPORT_FILENAME}")
proc_table
```