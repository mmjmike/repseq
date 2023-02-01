import yaml
import os
import pandas as pd
import dill
import json

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

def read_mixcr_clonoset(filename):
    clonoset = pd.read_csv(filename, sep="\t", dtype={'cloneId': int, 'readCount': int, 'readFraction': float,
                                                          'uniqueUMICount': int, 'uniqueUMIFraction': float,
                                                          'uniqueMoleculeCount': int, 'uniqueMoleculeFraction': float,
                                                          'cloneCount': int, 'cloneFraction': float,
                                                          'targetSequences': str, 'targetQualities': str,
                                                          'allVHitsWithScore': str, 'allDHitsWithScore': str,
                                                          'allJHitsWithScore': str, 'allCHitsWithScore': str,
                                                          'allVAlignments': str, 'allDAlignments': str,
                                                          'allJAlignments': str, 'allCAlignments': str,
                                                          'nSeqFR1': str, 'minQualFR1': str,
                                                          'nSeqCDR1': str, 'minQualCDR1': str,
                                                          'nSeqFR2': str, 'minQualFR2': str,
                                                          'nSeqCDR2': str, 'minQualCDR2': str,
                                                          'nSeqFR3': str, 'minQualFR3': str,
                                                          'nSeqCDR3': str, 'minQualCDR3': str,
                                                          'nSeqFR4': str, 'minQualFR4': str,
                                                          'aaSeqFR1': str, 'aaSeqCDR1': str,
                                                          'aaSeqFR2': str, 'aaSeqCDR2': str,
                                                          'aaSeqFR3': str, 'aaSeqCDR3': str,
                                                          'aaSeqFR4': str, 'refPoints': str})
    return clonoset

def read_json_report(sample_id, folder, report_type):
    filename = os.path.join(folder, f"{sample_id}.{report_type}.report.json")
    with open(filename) as data_file:    
        for jsonObj in data_file:
            report = json.loads(jsonObj)
    return report