# Page 2

## Code annotation examples


### Code blocks

Some `code` goes here

### Plain code blocks

```
def read_json_report(sample_id, folder, report_type):
    filename = os.path.join(folder, f"{sample_id}.{report_type}.report.json")
    with open(filename) as data_file:    
        for jsonObj in data_file:
            report = json.loads(jsonObj)
    return report
// some comment
```

#### Code for a specific language

Some more code with the `py` at the start

``` py
import sys
REPSEQ_PATH = '/home/mmyshkin/soft/repseq'
sys.path.append(REPSEQ_PATH)
from repseq import slurm
from repseq import io
from repseq import common_functions as cf
from repseq import clonosets as cl
from repseq import clustering
from repseq import mixcr as mx
from repseq import segment_usage as su
from repseq import stats
```

#### Code with a title

``` py title="import useful packages"
import os
import pandas as pd
from IPython.display import Image, display
import json
import re
import math
import random
import numpy as np
```

#### Add line numbers 


``` py linenums="1"
def shannon_wiener(list_of_numbers):
    list_of_numbers = list(list_of_numbers)
    total_size = sum(list_of_numbers)
    freqs = [s/total_size for s in list_of_numbers]
    diversity = len(list_of_numbers)
    sw = -sum([f*np.log(f) for f in freqs])
    sw_norm = sw/np.log(diversity)
    return sw, sw_norm, diversity
```

#### Highlighting lines

``` py hl_lines="5 6 7"
def shannon_wiener(list_of_numbers):
    list_of_numbers = list(list_of_numbers)
    total_size = sum(list_of_numbers)
    freqs = [s/total_size for s in list_of_numbers]
    diversity = len(list_of_numbers)
    sw = -sum([f*np.log(f) for f in freqs])
    sw_norm = sw/np.log(diversity)
    return sw, sw_norm, diversity
```

## Icons and Emojis

:smile:

:fontawesome-regular-face-laugh-wink:

:fontawesome-brands-twitter:{ .twitter }

:octicons-heart-fill-24:{ .heart }