import yaml
import os
import pandas as pd
import dill
import json
import zipfile

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

def read_clonoset(filename):
    file_name, file_extension = os.path.splitext(filename)

    d_types_mixcr = {'cloneId': int, 'readCount': int, 'readFraction': float,
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
                    'aaSeqFR4': str, 'refPoints': str
                    }

    d_types_vdjtools = {'cdr3aa': str, 'cdr3nt': str,
                        'v': str, 'd': str, 'j': str,
                        'count': int, 'freq': float,
                        'VEnd':int, 'DStart':int, 'DEnd':int, "JStart":int
                        }
    
    d_types_bioadaptive = {'nucleotide': str, 'aminoAcid': str,
                            'count (templates/reads)': int,
                            'frequencyCount (%)': float,
                            'vGeneName': str, 'dGeneName': str,
                            'jGeneName': str, 'cdr3Length': int,
                            'n1Index': int,'dIndex': int,
                            'n2Index': int,'jIndex': int
                            }
    datatypes = {**d_types_mixcr,**d_types_vdjtools, **d_types_bioadaptive}
    if file_extension == ".zip":
        archive = zipfile.ZipFile(filename, 'r')
        inner_filename = zipfile.ZipFile.namelist(archive)[0]
        filename = archive.open(inner_filename)
    clonoset = pd.read_csv(filename, sep="\t", dtype=datatypes)
    if 'count (templates/reads)' in clonoset.columns:
        clonoset = convert_bioadaptive_clonoset(clonoset)
    return clonoset

def convert_bioadaptive_clonoset(clonoset):
    # clonoset = clonoset.loc[clonoset.sequenceStatus == "In"]
    clonoset = clonoset[["count (templates/reads)", "frequencyCount (%)", "nucleotide", "aminoAcid",
                         "vMaxResolved", "dMaxResolved", "jMaxResolved", "n1Index", "dIndex", "n2Index", "jIndex"]]
    clonoset = clonoset.rename(columns={
        "count (templates/reads)": "count",
        "frequencyCount (%)": "freq",
        "nucleotide": "cdr3nt",
        "aminoAcid": "cdr3aa",
        "vMaxResolved": "v",
        "dMaxResolved": "d",
        "jMaxResolved": "j",
        "n1Index":"VEnd",
        "dIndex":"DStart",
        "n2Index":"DEnd",
        "jIndex":"JStart"
    })
    clonoset["v"] = clonoset["v"].apply(lambda x: recode_vdj_names_bioadaptive(x))
    clonoset["d"] = clonoset["d"].apply(lambda x: recode_vdj_names_bioadaptive(x))
    clonoset["j"] = clonoset["j"].apply(lambda x: recode_vdj_names_bioadaptive(x))
    clonoset["freq"] = clonoset["count"]/clonoset["count"].sum()
    clonoset["cdr3aa"] = clonoset["cdr3aa"].fillna("")
    return clonoset

def recode_vdj_names_bioadaptive(vdj_name):
    substitution_dict = {'TRAV2-1': 'TRAV2',
                        'TRAV3-1': 'TRAV3',
                        'TRAV4-1': 'TRAV4',
                        'TRAV5-1': 'TRAV5',
                        'TRAV6-1': 'TRAV6',
                        'TRAV7-1': 'TRAV7',
                        'TRAV10-1': 'TRAV10',
                        'TRAV11-1': 'TRAV11',
                        'TRAV14-1': 'TRAV14DV4',
                        'TRAV15-1': 'TRAV15',
                        'TRAV16-1': 'TRAV16',
                        'TRAV17-1': 'TRAV17',
                        'TRAV18-1': 'TRAV18',
                        'TRAV19-1': 'TRAV19',
                        'TRAV20-1': 'TRAV20',
                        'TRAV21-1': 'TRAV21',
                        'TRAV22-1': 'TRAV22',
                        'TRAV23-1': 'TRAV23DV6',
                        'TRAV24-1': 'TRAV24',
                        'TRAV25-1': 'TRAV25',
                        'TRAV27-1': 'TRAV27',
                        'TRAV28-1': 'TRAV28',
                        'TRAV29-1': 'TRAV29DV5',
                        'TRAV30-1': 'TRAV30',
                        'TRAV31-1': 'TRAV31',
                        'TRAV32-1': 'TRAV32',
                        'TRAV33-1': 'TRAV33',
                        'TRAV34-1': 'TRAV34',
                        'TRAV35-1': 'TRAV35',
                        'TRAV36-1': 'TRAV36DV7',
                        'TRAV38-2': 'TRAV38-2DV8',
                        'TRAV39-1': 'TRAV39',
                        'TRAV40-1': 'TRAV40',
                        'TRAV41-1': 'TRAV41',
                        'TRDV1-1': 'TRDV1',
                        'TRDV2-1': 'TRDV2',
                        'TRDV3-1': 'TRDV3',
                        'TRBVA-1': 'TRBVA',
                        'TRBV1-1': 'TRBV1',
                        'TRBV2-1': 'TRBV2',
                        'TRBV9-1': 'TRBV9',
                        'TRBV13-1': 'TRBV13',
                        'TRBV14-1': 'TRBV14',
                        'TRBV15-1': 'TRBV15',
                        'TRBV16-1': 'TRBV16',
                        'TRBV17-1': 'TRBV17',
                        'TRBV18-1': 'TRBV18',
                        'TRBV19-1': 'TRBV19',
                        'TRBV20': 'TRBV20-1',
                        'TRBV27-1': 'TRBV27',
                        'TRBV28-1': 'TRBV28',
                        'TRBV30-1': 'TRBV30'}
    if not isinstance(vdj_name, str):
        return "."
    vdj_name = vdj_name.replace("TCR", "TR").split("/")[0].split("*")[0].split("-or")[0]
    name_split = vdj_name.split("-")
    if len(name_split) > 1:
        subfamily_name = str(int(name_split[1]))
        try:
            family_name = str(int(name_split[0][-2:]))
            segment = name_split[0][:-2]
        except ValueError:
            family_name = name_split[0]
            segment = ""
        vdj_name = segment + family_name + "-" + subfamily_name
    if vdj_name in substitution_dict:
        vdj_name = substitution_dict[vdj_name]
    return vdj_name

def read_json_report(sample_id, folder, report_type):
    filename = os.path.join(folder, f"{sample_id}.{report_type}.report.json")
    with open(filename) as data_file:    
        for jsonObj in data_file:
            report = json.loads(jsonObj)
    return report