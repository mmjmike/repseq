import pandas as pd

from repseq.clone_filter import Filter


def test_filter_by_umi_falls_back_to_clone_count_when_umi_columns_absent():
    clonoset = pd.DataFrame(
        [
            {
                "cloneId": 0,
                "cloneCount": 10,
                "cloneFraction": 0.8,
                "nSeqCDR3": "TGTGCC",
                "aaSeqCDR3": "CASSLG",
                "allVHitsWithScore": "TRBV1*01(100)",
                "allDHitsWithScore": "TRBD1*01(10)",
                "allJHitsWithScore": "TRBJ1*01(80)",
            },
            {
                "cloneId": 1,
                "cloneCount": 2,
                "cloneFraction": 0.2,
                "nSeqCDR3": "TGTGCT",
                "aaSeqCDR3": "CAS*LG",
                "allVHitsWithScore": "TRBV2*01(90)",
                "allDHitsWithScore": "TRBD1*01(10)",
                "allJHitsWithScore": "TRBJ2*01(70)",
            },
        ]
    )

    result = Filter(functionality="f", by_umi=True).apply(clonoset)

    assert result["count"].tolist() == [10]
    assert result["freq"].tolist() == [1.0]
    assert result["cdr3aa"].tolist() == ["CASSLG"]
    assert result["v"].tolist() == ["TRBV1"]


def test_filter_treats_region_not_covered_as_nonfunctional_only_in_cdr3():
    clonoset = pd.DataFrame(
        [
            {
                "cloneId": 0,
                "cloneCount": 10,
                "cloneFraction": 0.5,
                "nSeqFR1": "region_not_covered",
                "aaSeqCDR1": "region_not_covered",
                "nSeqCDR3": "TGTGCC",
                "aaSeqCDR3": "CASSLG",
                "allVHitsWithScore": "TRBV1*01(100)",
                "allJHitsWithScore": "TRBJ1*01(80)",
            },
            {
                "cloneId": 1,
                "cloneCount": 6,
                "cloneFraction": 0.3,
                "nSeqCDR3": "region_not_covered",
                "aaSeqCDR3": "CASSQG",
                "allVHitsWithScore": "TRBV2*01(90)",
                "allJHitsWithScore": "TRBJ2*01(70)",
            },
            {
                "cloneId": 2,
                "cloneCount": 4,
                "cloneFraction": 0.2,
                "nSeqCDR3": "TGTGCT",
                "aaSeqCDR3": "region_not_covered",
                "allVHitsWithScore": "TRBV3*01(85)",
                "allJHitsWithScore": "TRBJ2*01(65)",
            },
        ]
    )

    functional = Filter(functionality="f").apply(clonoset)
    nonfunctional = Filter(functionality="n").apply(clonoset)

    assert functional["cdr3aa"].tolist() == ["CASSLG"]
    assert nonfunctional["cdr3aa"].tolist() == ["CASSQG", "region_not_covered"]


def test_filter_can_retain_alleles_for_all_segment_types():
    clonoset = pd.DataFrame(
        [
            {
                "cloneCount": 10,
                "cloneFraction": 1.0,
                "nSeqCDR3": "TGTGCC",
                "aaSeqCDR3": "CASSLG",
                "allVHitsWithScore": "TRBV6-2*000(742.2),TRBV6-3*001(742.2)",
                "allDHitsWithScore": "TRBD1*LONG(20),TRBD2*02(10)",
                "allJHitsWithScore": "TRBJ2-7*A1(80),TRBJ2-3*01(70)",
                "allCHitsWithScore": "TRBC2*XYZ(50),TRBC1*01(40)",
            }
        ]
    )

    result = Filter(retain_alleles=True).apply(clonoset)

    assert result.loc[0, ["v", "d", "j", "c"]].tolist() == [
        "TRBV6-2*000",
        "TRBD1*LONG",
        "TRBJ2-7*A1",
        "TRBC2*XYZ",
    ]


def test_filter_removes_alleles_by_default_and_uses_first_hit():
    clonoset = pd.DataFrame(
        [
            {
                "cloneCount": 10,
                "cloneFraction": 1.0,
                "nSeqCDR3": "TGTGCC",
                "aaSeqCDR3": "CASSLG",
                "allVHitsWithScore": "TRBV6-2*000(742.2),TRBV6-3*001(742.2)",
                "allJHitsWithScore": "TRBJ2-7*A1(80),TRBJ2-3*01(70)",
            }
        ]
    )

    result = Filter().apply(clonoset)

    assert result.loc[0, ["v", "j"]].tolist() == ["TRBV6-2", "TRBJ2-7"]


def test_filter_spawn_preserves_retain_alleles():
    assert Filter(retain_alleles=True).spawn().retain_alleles is True


def _filter_test_clonoset():
    return pd.DataFrame(
        [
            {
                "count": 10,
                "freq": 0.5,
                "cdr3nt": "TGTGCCAAA",
                "cdr3aa": "CASSLGQYF",
                "v": "TRBV7-8",
                "d": ".",
                "j": "TRBJ2-1",
            },
            {
                "count": 6,
                "freq": 0.3,
                "cdr3nt": "TGTGCCTTT",
                "cdr3aa": "CASSPGQYF",
                "v": "TRBV7-3",
                "d": ".",
                "j": "TRBJ2-7",
            },
            {
                "count": 4,
                "freq": 0.2,
                "cdr3nt": "TGTGGGAAA",
                "cdr3aa": "CATDAGNTIYF",
                "v": "TRBV20-1",
                "d": ".",
                "j": "TRBJ1-1",
            },
        ]
    )


def test_filter_white_list_exact_v_segment():
    result = Filter(white_list=[{"v": "TRBV7-8"}]).apply(_filter_test_clonoset())

    assert result["v"].tolist() == ["TRBV7-8"]


def test_filter_white_list_v_family_contains():
    result = Filter(white_list=[{"v": {"contains": "TRBV7"}}]).apply(_filter_test_clonoset())

    assert result["v"].tolist() == ["TRBV7-8", "TRBV7-3"]


def test_filter_white_list_wildcard_cdr3aa():
    result = Filter(white_list=[{"cdr3aa": "CASS*QYF"}]).apply(_filter_test_clonoset())

    assert result["cdr3aa"].tolist() == ["CASSLGQYF", "CASSPGQYF"]


def test_filter_white_list_regex_and_pattern_cdr3aa():
    regex_result = Filter(
        white_list=[{"cdr3aa": {"regex": r"CASS[LPS]GQYF"}}]
    ).apply(_filter_test_clonoset())
    pattern_result = Filter(
        white_list=[{"cdr3aa": {"pattern": r"CASS[LP]GQYF"}}]
    ).apply(_filter_test_clonoset())

    assert regex_result["cdr3aa"].tolist() == ["CASSLGQYF", "CASSPGQYF"]
    assert pattern_result["cdr3aa"].tolist() == ["CASSLGQYF", "CASSPGQYF"]


def test_filter_white_list_combined_cdr3nt_and_v():
    result = Filter(
        white_list=[{"cdr3nt": "TGTGCCAAA", "v": "TRBV7-8"}]
    ).apply(_filter_test_clonoset())

    assert result["cdr3nt"].tolist() == ["TGTGCCAAA"]


def test_filter_black_list_flexible_rule():
    result = Filter(black_list=[{"v": {"contains": "TRBV7"}}]).apply(_filter_test_clonoset())

    assert result["v"].tolist() == ["TRBV20-1"]


def test_filter_legacy_tuple_rule_still_works():
    result = Filter(white_list=[("CASSLGQYF", "TRBV7-8")]).apply(_filter_test_clonoset())

    assert result["cdr3aa"].tolist() == ["CASSLGQYF"]


def test_filter_tuple_rule_can_ignore_cdr3aa_and_select_v_only():
    result = Filter(white_list=[(None, "TRBV7-8")]).apply(_filter_test_clonoset())

    assert result["cdr3aa"].tolist() == ["CASSLGQYF"]


def test_filter_tuple_rule_can_ignore_cdr3aa_and_select_vj():
    result = Filter(white_list=[(None, "TRBV7-3", "TRBJ2-7")]).apply(_filter_test_clonoset())

    assert result["cdr3aa"].tolist() == ["CASSPGQYF"]


def test_filter_tuple_rule_can_ignore_v_and_select_cdr3aa_j():
    result = Filter(white_list=[("CATDAGNTIYF", None, "TRBJ1-1")]).apply(_filter_test_clonoset())

    assert result["v"].tolist() == ["TRBV20-1"]


def test_filter_tuple_rule_can_use_empty_string_as_ignored_position():
    result = Filter(white_list=[("", "TRBV7-3", "")]).apply(_filter_test_clonoset())

    assert result["cdr3aa"].tolist() == ["CASSPGQYF"]
