import os
import pandas as pd

from .clone_filter import Filter
from .io import read_clonoset


def save_to_vdjtools(samples_df, output_folder, cl_filter=None):
    metadata_list = []
    if cl_filter is None:
        cl_filter = Filter()

    for i,r in samples_df.iterrows():
        f = r["filename"]
        sample_id = r["sample_id"]
        basename = os.path.splitext(os.path.basename(f))[0]
        new_filename = f"vdjtools.{sample_id}.txt"
        new_path = os.path.join(output_folder, new_filename)
        metadata_list.append([new_filename, sample_id])
        clonoset = read_clonoset(f)
        clonoset = cl_filter.apply(clonoset)
        clonoset.to_csv(new_path, index=False, sep="\t")
    metadata_filename = os.path.join(output_folder, "metadata.txt")
    metadata = pd.DataFrame(metadata_list, columns=["#file.name", "sample.id"])
    metadata.to_csv(metadata_filename, index=False, sep="\t")
    print(f"Saved {len(metadata)} clonosets to: {output_folder}")
    print(f"Saved sample list to: {metadata_filename}")