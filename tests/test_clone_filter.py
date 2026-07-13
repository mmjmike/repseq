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
