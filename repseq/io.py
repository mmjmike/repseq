import yaml
import os
import pandas as pd
import dill
import json
import csv
import bz2
import gzip
import lzma
import zipfile
import io
import warnings
from contextlib import contextmanager
from .clonoset import Clonoset, standardize_to_vdjtools_columns
from .common_functions import extract_segment


_STRING_DTYPE = "string"
_INTEGER_DTYPE = "Int64"
_FLOAT_DTYPE = "float64"
_BOOLEAN_DTYPE = "boolean"

MIXCR_DTYPES = {
    "cloneId": _INTEGER_DTYPE,
    "readCount": _INTEGER_DTYPE,
    "readFraction": _FLOAT_DTYPE,
    "uniqueUMICount": _INTEGER_DTYPE,
    "uniqueUMIFraction": _FLOAT_DTYPE,
    "uniqueMoleculeCount": _INTEGER_DTYPE,
    "uniqueMoleculeFraction": _FLOAT_DTYPE,
    "cloneCount": _INTEGER_DTYPE,
    "cloneFraction": _FLOAT_DTYPE,
    "targetSequences": _STRING_DTYPE,
    "targetQualities": _STRING_DTYPE,
    "allVHitsWithScore": _STRING_DTYPE,
    "allDHitsWithScore": _STRING_DTYPE,
    "allJHitsWithScore": _STRING_DTYPE,
    "allCHitsWithScore": _STRING_DTYPE,
    "allVAlignments": _STRING_DTYPE,
    "allDAlignments": _STRING_DTYPE,
    "allJAlignments": _STRING_DTYPE,
    "allCAlignments": _STRING_DTYPE,
    "nSeqFR1": _STRING_DTYPE,
    "minQualFR1": _STRING_DTYPE,
    "nSeqCDR1": _STRING_DTYPE,
    "minQualCDR1": _STRING_DTYPE,
    "nSeqFR2": _STRING_DTYPE,
    "minQualFR2": _STRING_DTYPE,
    "nSeqCDR2": _STRING_DTYPE,
    "minQualCDR2": _STRING_DTYPE,
    "nSeqFR3": _STRING_DTYPE,
    "minQualFR3": _STRING_DTYPE,
    "nSeqCDR3": _STRING_DTYPE,
    "minQualCDR3": _STRING_DTYPE,
    "nSeqFR4": _STRING_DTYPE,
    "minQualFR4": _STRING_DTYPE,
    "aaSeqFR1": _STRING_DTYPE,
    "aaSeqCDR1": _STRING_DTYPE,
    "aaSeqFR2": _STRING_DTYPE,
    "aaSeqCDR2": _STRING_DTYPE,
    "aaSeqFR3": _STRING_DTYPE,
    "aaSeqCDR3": _STRING_DTYPE,
    "aaSeqFR4": _STRING_DTYPE,
    "refPoints": _STRING_DTYPE,
}

VDJTOOLS_DTYPES = {
    "count": _INTEGER_DTYPE,
    "freq": _FLOAT_DTYPE,
    "frequency": _FLOAT_DTYPE,
    "cdr3aa": _STRING_DTYPE,
    "cdr3nt": _STRING_DTYPE,
    "v": _STRING_DTYPE,
    "d": _STRING_DTYPE,
    "j": _STRING_DTYPE,
    "VEnd": _INTEGER_DTYPE,
    "DStart": _INTEGER_DTYPE,
    "DEnd": _INTEGER_DTYPE,
    "JStart": _INTEGER_DTYPE,
}

TRUST4_DTYPES = {
    "#count": _INTEGER_DTYPE,
    "count": _INTEGER_DTYPE,
    "frequency": _FLOAT_DTYPE,
    "CDR3nt": _STRING_DTYPE,
    "CDR3aa": _STRING_DTYPE,
    "CDR3_dna": _STRING_DTYPE,
    "CDR3_amino_acids": _STRING_DTYPE,
    "V": _STRING_DTYPE,
    "D": _STRING_DTYPE,
    "J": _STRING_DTYPE,
    "C": _STRING_DTYPE,
}

AIRR_STRING_COLUMNS = (
    "sequence_id", "sequence", "quality", "sequence_aa", "sequence_alignment",
    "quality_alignment", "sequence_alignment_aa", "germline_alignment",
    "germline_alignment_aa", "v_call", "d_call", "d2_call", "j_call", "c_call",
    "junction", "junction_aa", "np1", "np1_aa", "np2", "np2_aa", "np3", "np3_aa",
    "fwr1", "fwr1_aa", "cdr1", "cdr1_aa", "fwr2", "fwr2_aa", "cdr2",
    "cdr2_aa", "fwr3", "fwr3_aa", "cdr3", "cdr3_aa", "fwr4", "fwr4_aa",
    "v_cigar", "d_cigar", "d2_cigar", "j_cigar", "c_cigar",
    "v_sequence_alignment", "v_sequence_alignment_aa", "v_germline_alignment",
    "v_germline_alignment_aa", "d_sequence_alignment", "d_sequence_alignment_aa",
    "d_germline_alignment", "d_germline_alignment_aa", "d2_sequence_alignment",
    "d2_sequence_alignment_aa", "d2_germline_alignment", "d2_germline_alignment_aa",
    "j_sequence_alignment", "j_sequence_alignment_aa", "j_germline_alignment",
    "j_germline_alignment_aa", "c_sequence_alignment", "c_sequence_alignment_aa",
    "c_germline_alignment", "c_germline_alignment_aa", "locus", "locus_species",
    "cell_id", "clone_id", "repertoire_id", "reactivity_id", "reactivity_ref",
    "sample_processing_id", "data_processing_id", "rearrangement_type",
    "rearrangement_id", "rearrangement_set_id", "germline_database",
)
AIRR_INTEGER_COLUMNS = (
    "v_sequence_start", "v_sequence_end", "v_germline_start", "v_germline_end",
    "v_alignment_start", "v_alignment_end", "d_sequence_start", "d_sequence_end",
    "d_germline_start", "d_germline_end", "d_alignment_start", "d_alignment_end",
    "d2_sequence_start", "d2_sequence_end", "d2_germline_start", "d2_germline_end",
    "d2_alignment_start", "d2_alignment_end", "j_sequence_start", "j_sequence_end",
    "j_germline_start", "j_germline_end", "j_alignment_start", "j_alignment_end",
    "c_sequence_start", "c_sequence_end", "c_germline_start", "c_germline_end",
    "c_alignment_start", "c_alignment_end", "cdr1_start", "cdr1_end", "cdr2_start",
    "cdr2_end", "cdr3_start", "cdr3_end", "fwr1_start", "fwr1_end", "fwr2_start",
    "fwr2_end", "fwr3_start", "fwr3_end", "fwr4_start", "fwr4_end",
    "junction_length", "junction_aa_length", "np1_length", "np2_length", "np3_length",
    "n1_length", "n2_length", "n3_length", "p3v_length", "p5d_length", "p3d_length",
    "p5d2_length", "p3d2_length", "p5j_length", "duplicate_count", "consensus_count",
    "umi_count", "d_frame", "d2_frame",
)
AIRR_FLOAT_COLUMNS = (
    "v_score", "v_identity", "v_support", "d_score", "d_identity", "d_support",
    "d2_score", "d2_identity", "d2_support", "j_score", "j_identity", "j_support",
    "c_score", "c_identity", "c_support",
)
AIRR_BOOLEAN_COLUMNS = (
    "rev_comp", "productive", "vj_in_frame", "stop_codon", "complete_vdj",
    "v_frameshift", "j_frameshift",
)
AIRR_DTYPES = {
    **{column: _STRING_DTYPE for column in AIRR_STRING_COLUMNS},
    **{column: _INTEGER_DTYPE for column in AIRR_INTEGER_COLUMNS},
    **{column: _FLOAT_DTYPE for column in AIRR_FLOAT_COLUMNS},
    **{column: _BOOLEAN_DTYPE for column in AIRR_BOOLEAN_COLUMNS},
}

IGBLAST_DTYPES = {
    **{column.upper(): datatype for column, datatype in AIRR_DTYPES.items()},
    "SEQUENCE_INPUT": _STRING_DTYPE,
    "FUNCTIONAL": _STRING_DTYPE,
}

CLONOSET_FORMAT_DTYPES = {
    "MiXCR3": MIXCR_DTYPES,
    "MiXCR4": MIXCR_DTYPES,
    "VDJtools": VDJTOOLS_DTYPES,
    "AIRR": AIRR_DTYPES,
    "IgBLAST": IGBLAST_DTYPES,
    "TRUST4": TRUST4_DTYPES,
}


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


@contextmanager
def _open_clonoset_text(filename):
    path = os.fspath(filename)
    lower_path = path.lower()
    if lower_path.endswith(".zip"):
        with zipfile.ZipFile(path) as archive:
            members = [name for name in archive.namelist() if not name.endswith("/")]
            if not members:
                raise ValueError(f"Zip archive contains no files: {path}")
            with archive.open(members[0]) as binary_file:
                with io.TextIOWrapper(binary_file, encoding="utf-8-sig", newline="") as text_file:
                    yield text_file
        return

    opener = open
    if lower_path.endswith(".gz"):
        opener = gzip.open
    elif lower_path.endswith(".bz2"):
        opener = bz2.open
    elif lower_path.endswith((".xz", ".lzma")):
        opener = lzma.open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as text_file:
        yield text_file


def _clonoset_separator(filename, line_count=10):
    with _open_clonoset_text(filename) as text_file:
        sample = "".join(text_file.readline() for _ in range(line_count))
    if not sample:
        raise ValueError(f"Clonoset file is empty: {filename}")
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,").delimiter
    except csv.Error:
        header = sample.splitlines()[0]
        if "\t" in header:
            return "\t"
        if "," in header:
            return ","
        raise ValueError(f"Could not detect tab or comma separator in: {filename}")


def _clonoset_columns(filename, separator):
    with _open_clonoset_text(filename) as text_file:
        return list(pd.read_csv(text_file, sep=separator, nrows=0).columns)


def _detect_clonoset_format(columns):
    column_set = set(columns)
    lower_columns = {str(column).lower() for column in columns}

    if {"#count", "cdr3nt", "cdr3aa"}.issubset(lower_columns) or {
        "cdr3_dna", "cdr3_amino_acids"
    }.issubset(lower_columns):
        return "TRUST4"
    if "readCount" in column_set or "readFraction" in column_set:
        return "MiXCR4"
    if "cloneCount" in column_set or "cloneFraction" in column_set:
        return "MiXCR3"
    if {"count", "freq"}.issubset(lower_columns) and (
        "cdr3nt" in lower_columns or "cdr3aa" in lower_columns
    ):
        return "VDJtools"
    if "SEQUENCE_ID" in column_set or "SEQUENCE_INPUT" in column_set:
        return "IgBLAST"
    if {"sequence_id", "v_call", "j_call"}.issubset(lower_columns):
        return "AIRR"
    return None


def _read_typed_clonoset(filename, separator, columns, clonoset_format):
    dtype_registry = CLONOSET_FORMAT_DTYPES.get(clonoset_format, {})
    datatypes = {column: dtype_registry[column] for column in columns if column in dtype_registry}
    read_options = {
        "sep": separator,
        "dtype": datatypes,
        "true_values": ["T", "TRUE", "True", "true"],
        "false_values": ["F", "FALSE", "False", "false"],
    }
    with _open_clonoset_text(filename) as text_file:
        return pd.read_csv(text_file, **read_options)


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
    Reads `csv`, `tsv`, `txt`, compressed files, and the first file in a zip.
    The separator and clonoset format are detected from the first few lines,
    then format-specific pandas dtypes are used for the full read.
    Supported typed formats are MiXCR3, MiXCR4, VDJtools, AIRR, IgBLAST,
    and TRUST4. Unknown column layouts are still read with inferred dtypes.
    
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
    
    
    separator = _clonoset_separator(filename)
    columns = _clonoset_columns(filename, separator)
    clonoset_format = _detect_clonoset_format(columns)
    clonoset = _read_typed_clonoset(filename, separator, columns, clonoset_format)
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
            written as `vdjtools.{sample_id}.txt`. All additional columns are
            preserved in `metadata.txt`.
        output_folder (str): folder in which to save converted clonosets and
            `metadata.txt`.
        cl_filter (Filter, optional): filter used to convert and filter each
            clonoset. Defaults to `Filter()`.
        force_overwrite (bool): if `False`, check all target filenames before
            writing. If any target already exists, print a warning and do not
            write any files. Set to `True` to overwrite existing clonoset files
            and merge new rows into existing `metadata.txt`, preserving rows
            for output files that still exist.

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
    output_filenames = []
    targets = []

    for _, row in samples_df.iterrows():
        new_filename = _vdjtools_filename(row, has_chain)
        new_path = os.path.join(output_folder, new_filename)
        targets.append(new_path)
        output_filenames.append(new_filename)

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

    metadata = _vdjtools_metadata(
        samples_df,
        output_filenames,
        output_folder,
        metadata_filename,
        force_overwrite=force_overwrite,
    )
    metadata.to_csv(metadata_filename, index=False, sep="\t")
    print(f"Saved {len(output_filenames)} clonosets to: {output_folder}")
    print(f"Saved sample list to: {metadata_filename}")
    return metadata


def _new_vdjtools_metadata(samples_df, output_filenames):
    metadata = samples_df.copy().reset_index(drop=True)
    if "#file.name" in metadata.columns:
        metadata = metadata.drop(columns=["#file.name"])
    if "filename" in metadata.columns:
        metadata = metadata.rename(columns={"filename": "original_filename"})
    metadata.insert(0, "#file.name", output_filenames)
    if "sample.id" not in metadata.columns:
        metadata.insert(1, "sample.id", metadata["sample_id"])
    return metadata


def _vdjtools_metadata(samples_df, output_filenames, output_folder, metadata_filename, force_overwrite=False):
    new_metadata = _new_vdjtools_metadata(samples_df, output_filenames)
    if not force_overwrite or not os.path.exists(metadata_filename):
        return new_metadata

    existing_metadata = pd.read_csv(metadata_filename, sep="\t")
    if "filename" in existing_metadata.columns and "original_filename" not in existing_metadata.columns:
        existing_metadata = existing_metadata.rename(columns={"filename": "original_filename"})
    if "#file.name" not in existing_metadata.columns:
        return new_metadata

    new_filenames = set(output_filenames)
    keep_existing = existing_metadata["#file.name"].apply(
        lambda filename: (
            filename not in new_filenames
            and os.path.exists(os.path.join(output_folder, str(filename)))
        )
    )
    existing_metadata = existing_metadata.loc[keep_existing]
    return pd.concat([existing_metadata, new_metadata], ignore_index=True, sort=False)


AIRR_SOURCE_ALIASES = {
    "sequence_id": ("sequence_id", "SEQUENCE_ID", "cloneId", "clone_id", "cid"),
    "sequence": ("sequence", "SEQUENCE_INPUT", "targetSequences"),
    "quality": ("quality", "targetQualities"),
    "sequence_aa": ("sequence_aa",),
    "rev_comp": ("rev_comp",),
    "productive": ("productive", "FUNCTIONAL"),
    "vj_in_frame": ("vj_in_frame",),
    "stop_codon": ("stop_codon",),
    "complete_vdj": ("complete_vdj",),
    "v_call": ("v_call", "V_CALL", "allVHitsWithScore", "v", "V", "bestVGene"),
    "d_call": ("d_call", "D_CALL", "allDHitsWithScore", "d", "D", "bestDGene"),
    "j_call": ("j_call", "J_CALL", "allJHitsWithScore", "j", "J", "bestJGene"),
    "c_call": ("c_call", "C_CALL", "allCHitsWithScore", "c", "C", "bestCGene"),
    "junction": (
        "junction", "JUNCTION", "nSeqCDR3", "cdr3nt", "CDR3nt", "CDR3_dna",
        "cdr3.nucleotide.sequence",
    ),
    "junction_aa": (
        "junction_aa", "JUNCTION_AA", "aaSeqCDR3", "cdr3aa", "CDR3aa",
        "CDR3_amino_acids", "cdr3.amino.acid.sequence",
    ),
    "fwr1": ("fwr1", "nSeqFR1"),
    "fwr1_aa": ("fwr1_aa", "aaSeqFR1"),
    "cdr1": ("cdr1", "nSeqCDR1"),
    "cdr1_aa": ("cdr1_aa", "aaSeqCDR1"),
    "fwr2": ("fwr2", "nSeqFR2"),
    "fwr2_aa": ("fwr2_aa", "aaSeqFR2"),
    "cdr2": ("cdr2", "nSeqCDR2"),
    "cdr2_aa": ("cdr2_aa", "aaSeqCDR2"),
    "fwr3": ("fwr3", "nSeqFR3"),
    "fwr3_aa": ("fwr3_aa", "aaSeqFR3"),
    "fwr4": ("fwr4", "nSeqFR4"),
    "fwr4_aa": ("fwr4_aa", "aaSeqFR4"),
    "duplicate_count": (
        "duplicate_count", "readCount", "cloneCount", "count", "#count",
        "count (templates/reads)",
    ),
    "umi_count": ("umi_count", "uniqueUMICount", "uniqueMoleculeCount"),
}


def _series_has_data(series):
    populated = series.notna()
    if pd.api.types.is_string_dtype(series.dtype) or series.dtype == object:
        populated &= series.astype("string").str.strip().ne("")
    return bool(populated.any())


def _first_populated_column(clonoset, aliases):
    for column in aliases:
        if column in clonoset.columns and _series_has_data(clonoset[column]):
            return column
    return None


def _clean_airr_gene_calls(series):
    return series.astype("string").str.replace(r"\([^,()]*\)", "", regex=True)


def _to_airr_boolean(series):
    if pd.api.types.is_bool_dtype(series.dtype):
        return series.astype(_BOOLEAN_DTYPE)
    normalized = series.astype("string").str.strip().str.lower()
    return normalized.map({
        "t": True,
        "true": True,
        "1": True,
        "yes": True,
        "functional": True,
        "productive": True,
        "f": False,
        "false": False,
        "0": False,
        "no": False,
        "non-functional": False,
        "nonfunctional": False,
        "unproductive": False,
    }).astype(_BOOLEAN_DTYPE)


def _infer_productive(junction_aa):
    amino_acids = junction_aa.astype("string")
    return (
        amino_acids.notna()
        & amino_acids.str.strip().ne("")
        & ~amino_acids.str.contains(r"\*|_", na=True)
    ).astype(_BOOLEAN_DTYPE)


def _infer_locus(airr_clonoset):
    calls = None
    for column in ("v_call", "j_call", "c_call"):
        if column in airr_clonoset.columns:
            calls = airr_clonoset[column] if calls is None else calls.fillna(airr_clonoset[column])
    if calls is None:
        return None
    return calls.astype("string").str.extract(r"\b(IG[HKL]|TR[ABDG])", expand=False)


def _to_airr_clonoset(clonoset):
    airr_columns = {}
    for target, aliases in AIRR_SOURCE_ALIASES.items():
        source = _first_populated_column(clonoset, aliases)
        if source is None:
            continue
        values = clonoset[source].copy()
        if target.endswith("_call") and source.startswith("all") and source.endswith("HitsWithScore"):
            values = _clean_airr_gene_calls(values)
        airr_columns[target] = values

    for column in AIRR_DTYPES:
        if column not in airr_columns and column in clonoset.columns and _series_has_data(clonoset[column]):
            airr_columns[column] = clonoset[column].copy()

    airr_clonoset = pd.DataFrame(airr_columns, index=clonoset.index)
    if "productive" not in airr_clonoset.columns and "junction_aa" in airr_clonoset.columns:
        airr_clonoset["productive"] = _infer_productive(airr_clonoset["junction_aa"])
    if "junction_length" not in airr_clonoset.columns and "junction" in airr_clonoset.columns:
        airr_clonoset["junction_length"] = airr_clonoset["junction"].astype("string").str.len().astype(_INTEGER_DTYPE)
    if "junction_aa_length" not in airr_clonoset.columns and "junction_aa" in airr_clonoset.columns:
        airr_clonoset["junction_aa_length"] = airr_clonoset["junction_aa"].astype("string").str.len().astype(_INTEGER_DTYPE)
    if "locus" not in airr_clonoset.columns:
        locus = _infer_locus(airr_clonoset)
        if locus is not None and _series_has_data(locus):
            airr_clonoset["locus"] = locus

    for column in list(airr_clonoset.columns):
        datatype = AIRR_DTYPES[column]
        if datatype == _BOOLEAN_DTYPE:
            airr_clonoset[column] = _to_airr_boolean(airr_clonoset[column])
        else:
            airr_clonoset[column] = airr_clonoset[column].astype(datatype)

    ordered_columns = [column for column in AIRR_DTYPES if column in airr_clonoset.columns]
    return airr_clonoset.loc[:, ordered_columns].reset_index(drop=True)


def _airr_filename(source_filename):
    return f"{os.path.basename(os.fspath(source_filename))}.airr.tsv"


def _airr_metadata(samples_df, output_filenames, output_folder, metadata_filename, force_overwrite=False):
    new_metadata = samples_df.copy().reset_index(drop=True)
    new_metadata["filename"] = output_filenames
    if not force_overwrite or not os.path.exists(metadata_filename):
        return new_metadata

    existing_metadata = pd.read_csv(metadata_filename)
    if "filename" not in existing_metadata.columns:
        return new_metadata
    new_filenames = set(output_filenames)
    keep_existing = existing_metadata["filename"].apply(
        lambda filename: (
            filename not in new_filenames
            and os.path.exists(os.path.join(output_folder, str(filename)))
        )
    )
    existing_metadata = existing_metadata.loc[keep_existing]
    return pd.concat([existing_metadata, new_metadata], ignore_index=True, sort=False)


def save_to_airr(samples_df, output_folder, cl_filter=None, force_overwrite=False):
    """
    Save clonosets as AIRR rearrangement TSV files and write `metadata.csv`.

    Output filenames are the basename of each source filename with an
    additional `.airr.tsv` suffix. Metadata preserves the original
    `samples_df` columns and replaces only `filename` with that output basename.
    AIRR columns are copied or derived when supported by source data; columns
    with no data are omitted.

    Args:
        samples_df (pd.DataFrame): table containing a `filename` column.
        output_folder (str): destination for AIRR TSV files and `metadata.csv`.
        cl_filter (Filter, optional): filter applied before AIRR conversion.
            The default preserves source columns while applying no filtering.
        force_overwrite (bool): use the same all-or-nothing conflict behavior
            as `save_to_vdjtools`; when true, overwrite and merge metadata rows
            for output files that still exist.

    Returns:
        pd.DataFrame or None: AIRR metadata, or `None` when conflicts prevent
            writing.
    """
    from .clone_filter import Filter

    if "filename" not in samples_df.columns:
        raise ValueError("samples_df must contain 'filename' column")
    if cl_filter is None:
        cl_filter = Filter(convert=False)

    os.makedirs(output_folder, exist_ok=True)
    output_filenames = [_airr_filename(filename) for filename in samples_df["filename"]]
    output_paths = [os.path.join(output_folder, filename) for filename in output_filenames]
    metadata_filename = os.path.join(output_folder, "metadata.csv")
    targets = output_paths + [metadata_filename]
    conflicts = [filename for filename in targets if os.path.exists(filename)]
    if conflicts and not force_overwrite:
        print("WARNING! The following output files already exist:")
        for filename in conflicts:
            print(filename)
        print("No files were written. Set force_overwrite=True to overwrite existing files.")
        return None

    for (_, row), output_path in zip(samples_df.iterrows(), output_paths):
        clonoset = read_clonoset(row["filename"])
        clonoset = cl_filter.apply(clonoset)
        airr_clonoset = _to_airr_clonoset(clonoset)
        airr_clonoset.to_csv(output_path, index=False, sep="\t")

    metadata = _airr_metadata(
        samples_df,
        output_filenames,
        output_folder,
        metadata_filename,
        force_overwrite=force_overwrite,
    )
    metadata.to_csv(metadata_filename, index=False)
    print(f"Saved {len(output_filenames)} clonosets to: {output_folder}")
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
