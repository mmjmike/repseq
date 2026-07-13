import yaml
import os
import pandas as pd
import dill
import json
import zipfile
import requests
import io
import warnings
from .clonoset import Clonoset, standardize_to_vdjtools_columns
from .common_functions import extract_segment


def read_ngsik_metadata(folder, filename="metadata.yaml", verbose=True):
    
    """
    Reads NGSiK metadata from a given folder and converts to `pd.DataFrame`. By default 
    it searches for `metadata.yaml` file in this folder and extracts the table.
    
    Args:
        folder (str): path to NGSiK folder
        filename (str): NGSiK metadata filename
        verbose (bool): verbosity for 
    
    Returns:
        sample_df (pd.DataFrame): extracted DataFrame from metadata
    
    """
    
    
    most_important_columns = ["sample_id", "R1", "R2","libraryPerson", "projectPerson", "projectName", "species", "miNNNPattern", "SMPL", "mix_id", "preset", "startingMaterial", "libraryType"]
    yaml_filename = os.path.join(folder, filename)
    
    if os.path.isfile(yaml_filename):
        with open(yaml_filename, "r") as stream:
            try:
                metadata_dict =yaml.safe_load(stream)
        #         pd.io.json.json_normalize(metadata_dict, "file", "samples", errors='ignore')
            except yaml.YAMLError as exc:
                print(exc)

        df = pd.json_normalize(metadata_dict)
        df = df.explode("file")
        df = pd.concat([df.drop(['file'], axis=1), df['file'].apply(pd.Series)], axis=1)
        df = df.explode("samples")
        df = pd.concat([df.drop(['samples'], axis=1), df['samples'].apply(pd.Series)], axis=1)
        if 'patternGroupValues' in df.columns:
            df = pd.concat([df.drop(['patternGroupValues'], axis=1), df['patternGroupValues'].apply(pd.Series)], axis=1)
        df["R1"] = df["R1"].apply(lambda x: os.path.join(folder, x))
        df["R2"] = df["R2"].apply(lambda x: os.path.join(folder, x))
        df = df.rename(columns={"name": "sample_id"})
        
        for col_name in most_important_columns[::-1]:
            if col_name in df.columns:
                first_column = df.pop(col_name) 
                df.insert(0, col_name, first_column)
        
        return df.reset_index(drop=True)
    else:
        if verbose:
            print(f"Metadata file '{yaml_filename}' not found. Nothing to return")
        return pd.DataFrame()


def read_yaml_metadata(folder, filename="metadata.yaml", verbose=True):
    """
    Deprecated alias for `read_ngsik_metadata`.

    Use `read_ngsik_metadata` instead. This alias will be removed in a future
    version.
    """
    warnings.warn(
        "`read_yaml_metadata` is deprecated; use `read_ngsik_metadata` instead. "
        "`read_yaml_metadata` will be removed in a future version.",
        DeprecationWarning,
        stacklevel=2,
    )
    return read_ngsik_metadata(folder, filename=filename, verbose=verbose)


def read_clonoset(
    filename,
    *,
    as_clonoset=False,
    standardize=False,
    chain=None,
    sample_id=None,
    metadata=None,
    by_umi=False,
    validate=True,
):
    """
    Reads generic clonoset files. 
    Easyly reads `csv`, `tsv`, `txt` or `gz` files.
    Reads first found file inside `zip` files.
    Clonosets should be in tab-separated format: MiXCR (v3 or v4), vdjtools, Bioadaptive
    
    Args:
        filename (str): path to clonoset file
        as_clonoset (bool): If `True`, return a `Clonoset` object instead of
            a plain pandas DataFrame. Defaults to `False` for backward
            compatibility.
        standardize (bool): If `True`, return a table with canonical
            VDJtools-like column names (`freq`, `count`, `cdr3nt`, `cdr3aa`,
            `v`, `d`, `j`). This is applied automatically when
            `as_clonoset=True`.
        chain (str, optional): Receptor chain metadata for `Clonoset`.
        sample_id (str, optional): Sample identifier metadata for `Clonoset`.
        metadata (dict, optional): Additional metadata for `Clonoset`.
        by_umi (bool): If `True`, use UMI or molecule count/fraction columns
            for canonical `count` and `freq` when standardizing.
        validate (bool): If `True`, validate the minimal `Clonoset` schema
            when `as_clonoset=True`.

    Returns:
        clonoset (pd.DataFrame or Clonoset): DataFrame representation of the
            clonoset, or a `Clonoset` object when `as_clonoset=True`.
            Bioadaptive clonosets are converted to vdjtools-like format.
    """
    
    
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
                        'CDR3aa': str, 'CDR3nt': str,
                        'V': str, 'D': str, 'J': str,
                        'C': str, "frequency": float#,
                        #'count': int, 'freq': float#,
                        #'VEnd':int, 'DStart':int, 'DEnd':int, "JStart":int
                        }
    
    d_types_bioadaptive = {'nucleotide': str, 'aminoAcid': str,
                            'count (templates/reads)': int,
                            'frequencyCount (%)': float,
                            'count': int,
                            'frequencyCount': float,
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
    if "nucleotide" in clonoset.columns and "aminoAcid" in clonoset.columns:
        clonoset = convert_bioadaptive_clonoset(clonoset)
    if as_clonoset:
        return Clonoset(
            clonoset,
            chain=chain,
            sample_id=sample_id,
            metadata=metadata or {},
            standardize=True,
            by_umi=by_umi,
            validate=validate,
        )
    if standardize:
        clonoset = standardize_to_vdjtools_columns(
            clonoset,
            by_umi=by_umi,
            copy=False,
        )
    return clonoset


def _vdjtools_filename(row, has_chain):
    sample_id = row["sample_id"]
    if has_chain:
        return f"vdjtools.{sample_id}.{row['chain']}.txt"
    return f"vdjtools.{sample_id}.txt"


def save_to_vdjtools(samples_df, output_folder, cl_filter=None, force_overwrite=False):
    """
    Save clonosets in VDJtools-like format and write VDJtools metadata.

    Args:
        samples_df (pd.DataFrame): table with `sample_id` and `filename`
            columns. If a `chain` column is present, output filenames are
            written as `vdjtools.{sample_id}.{chain}.txt`; otherwise they are
            written as `vdjtools.{sample_id}.txt`.
        output_folder (str): folder in which to save converted clonosets and
            `metadata.txt`.
        cl_filter (Filter, optional): filter used to convert and filter each
            clonoset. Defaults to `Filter()`.
        force_overwrite (bool): if `False`, check all target filenames before
            writing. If any target already exists, print a warning and do not
            write any files. Set to `True` to overwrite existing files.

    Returns:
        pd.DataFrame or None: VDJtools metadata table when files are written,
            otherwise `None` if conflicting files were found.
    """
    from .clone_filter import Filter

    if "sample_id" not in samples_df.columns:
        raise ValueError("samples_df must contain 'sample_id' column")
    if "filename" not in samples_df.columns:
        raise ValueError("samples_df must contain 'filename' column")
    if cl_filter is None:
        cl_filter = Filter()

    os.makedirs(output_folder, exist_ok=True)
    has_chain = "chain" in samples_df.columns
    metadata_list = []
    targets = []

    for _, row in samples_df.iterrows():
        new_filename = _vdjtools_filename(row, has_chain)
        new_path = os.path.join(output_folder, new_filename)
        targets.append(new_path)
        metadata_list.append([new_filename, row["sample_id"]])

    metadata_filename = os.path.join(output_folder, "metadata.txt")
    targets.append(metadata_filename)
    conflicts = [filename for filename in targets if os.path.exists(filename)]
    if conflicts and not force_overwrite:
        print("WARNING! The following output files already exist:")
        for filename in conflicts:
            print(filename)
        print("No files were written. Set force_overwrite=True to overwrite existing files.")
        return None

    for (_, row), new_path in zip(samples_df.iterrows(), targets[:-1]):
        clonoset = read_clonoset(row["filename"])
        clonoset = cl_filter.apply(clonoset)
        clonoset.to_csv(new_path, index=False, sep="\t")

    metadata = pd.DataFrame(metadata_list, columns=["#file.name", "sample.id"])
    metadata.to_csv(metadata_filename, index=False, sep="\t")
    print(f"Saved {len(metadata)} clonosets to: {output_folder}")
    print(f"Saved sample list to: {metadata_filename}")
    return metadata


def read_json_report(sample_id, folder, report_type):
    """
    Reads MiXCR4 json reports into a Python mixed data structure.
    This function takes the last json record, if for example MiXCR adds up several records 
    to json file (it happens, when the program is rerun several times on the same data).
    Program also includes cases when Sample-barcodes are used.

    Args:
        sample_id (str): sample_id used when running the MiXCR program
        folder (str): folder in which the MiXCR output is stored
        report_type (str): align, refine, assemble

    Returns:
        report (dict): mixed dict/list python structure, representing the json report
    """


    filename = os.path.join(folder, f"{sample_id}.{report_type}.report.json")
    if "." in sample_id:
        sample_id2 = ".".join(sample_id.split(".")[:-1])
        filename2 = os.path.join(folder, f"{sample_id2}.{report_type}.report.json")
        try:
            report = open_json_report(filename)
        except FileNotFoundError:
            report = open_json_report(filename2)
    else:
        report = open_json_report(filename)
    return report

def open_json_report(filename):
    """
    Supporting function for `read_json_report`. Reads the last record from json file.
    """

    with open(filename) as data_file:
        contents = data_file.read()

    decoder = json.JSONDecoder()
    report = None
    index = 0
    while index < len(contents):
        while index < len(contents) and contents[index].isspace():
            index += 1
        if index >= len(contents):
            break

        try:
            report, index = decoder.raw_decode(contents, index)
        except json.JSONDecodeError:
            next_object = contents.find("{", index + 1)
            next_array = contents.find("[", index + 1)
            next_positions = [pos for pos in [next_object, next_array] if pos != -1]
            if len(next_positions) == 0:
                break
            index = min(next_positions)

    if report is None:
        raise json.JSONDecodeError("Expecting JSON value", contents, 0)
    return report


def convert_bioadaptive_clonoset(clonoset):
    # clonoset = clonoset.loc[clonoset.sequenceStatus == "In"]
    putative_colnames = ["count (templates/reads)", "count", "frequencyCount", "frequencyCount (%)", "nucleotide", "aminoAcid",
                         "vMaxResolved", "dMaxResolved", "jMaxResolved", "n1Index", "dIndex", "n2Index", "jIndex"]

    clonoset = clonoset[[c for c in putative_colnames if c in clonoset.columns]]
    clonoset = clonoset.rename(columns={
        "count (templates/reads)": "count",
        "frequencyCount (%)": "freq",
        "frequencyCount": "freq",
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
                        'TRBV1-1': 'TRBV1',
                        'TRBV2-1': 'TRBV2',
                        'TRBV03' : 'TRBV3',
                        'TRBV04' : 'TRBV4',
                        'TRBV05' : 'TRBV5',
                        'TRBV06' : 'TRBV6',
                        'TRBV07' : 'TRBV7',
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
        subfamily_name = str(int(name_split[1][:2])) + name_split[1][2:]
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

def save_dill_dump(obj, filename):
    with open(filename, 'wb') as f: 
        dill.dump(obj, f)

def read_dill_dump(filename):
    with open(filename, 'rb') as f:
        return dill.load(f)

# def read_mixcr_clonoset(filename):
#     # DEPRECATED
#     clonoset = pd.read_csv(filename, sep="\t", dtype={'cloneId': int, 'readCount': int, 'readFraction': float,
#                                                           'uniqueUMICount': int, 'uniqueUMIFraction': float,
#                                                           'uniqueMoleculeCount': int, 'uniqueMoleculeFraction': float,
#                                                           'cloneCount': int, 'cloneFraction': float,
#                                                           'targetSequences': str, 'targetQualities': str,
#                                                           'allVHitsWithScore': str, 'allDHitsWithScore': str,
#                                                           'allJHitsWithScore': str, 'allCHitsWithScore': str,
#                                                           'allVAlignments': str, 'allDAlignments': str,
#                                                           'allJAlignments': str, 'allCAlignments': str,
#                                                           'nSeqFR1': str, 'minQualFR1': str,
#                                                           'nSeqCDR1': str, 'minQualCDR1': str,
#                                                           'nSeqFR2': str, 'minQualFR2': str,
#                                                           'nSeqCDR2': str, 'minQualCDR2': str,
#                                                           'nSeqFR3': str, 'minQualFR3': str,
#                                                           'nSeqCDR3': str, 'minQualCDR3': str,
#                                                           'nSeqFR4': str, 'minQualFR4': str,
#                                                           'aaSeqFR1': str, 'aaSeqCDR1': str,
#                                                           'aaSeqFR2': str, 'aaSeqCDR2': str,
#                                                           'aaSeqFR3': str, 'aaSeqCDR3': str,
#                                                           'aaSeqFR4': str, 'refPoints': str})
#     return clonoset
