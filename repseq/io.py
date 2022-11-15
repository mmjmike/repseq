import yaml
import os
import pandas as pd
import dill

def read_yaml_metadata(folder, filename="metadata.yaml"):
    yaml_filename = os.path.join(folder, filename)
    with open(yaml_filename, "r") as stream:
        try:
            metadata_dict =yaml.safe_load(stream)
    #         pd.io.json.json_normalize(metadata_dict, "file", "samples", errors='ignore')
        except yaml.YAMLError as exc:
            print(exc)

    all_data = []
    for entry in metadata_dict["file"]:
        samples = entry.pop("samples")[0]
        samples.update(entry)
        all_data.append(samples)

    yaml_metadata = pd.DataFrame(all_data).rename(columns={"name": "sample_id"})
    return yaml_metadata

def save_dill_dump(obj, filename):
    with open(filename, 'wb') as f: 
        dill.dump(obj, f)

def read_dill_dump(filename):
    with open(filename, 'rb') as f:
        return dill.load(f)