from repseq import common_functions
from repseq.common_functions import center_5, diversity_metrics, run_parallel_calculation
import pytest

def test_center_5():
    assert center_5("ABCD") == "ABCD"
    assert center_5("ABCDEF") == "ABCDE"
    assert center_5("ABCDEFGHK") == "CDEFG"


def divide(a,b):
    if b == 0:
        raise ValueError("Cannot divide by zero")
    return a/b


def test_some_function():
    with pytest.raises(ValueError, match="Cannot divide by zero"):
        divide(10,0)


def identity_task(value):
    return value


def test_run_parallel_calculation_prints_worker_message_before_progress(monkeypatch, capsys):
    class FakeExecutor:
        def __init__(self, max_workers=None):
            self.max_workers = max_workers

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, traceback):
            return False

        def map(self, function, tasks):
            return [function(task) for task in tasks]

    monkeypatch.setattr(common_functions.concurrent.futures, "ProcessPoolExecutor", FakeExecutor)

    result = run_parallel_calculation(
        identity_task,
        [1, 2],
        "TestProgram",
        object_name="item(s)",
        verbose=True,
        cpu=None,
    )

    captured = capsys.readouterr().out
    assert result == [1, 2]
    assert "Using None cores" not in captured
    assert captured.startswith("Using default number of worker processes\nTestProgram")


def test_diversity_metrics_include_expected_metrics_and_order():
    result = diversity_metrics([5, 3, 2])

    assert list(result)[:5] == [
        "diversity",
        "norm_shannon_wiener",
        "clonality",
        "shannon_wiener",
        "chao1",
    ]
    assert result["diversity"] == 3
    assert result["richness"] == 3
    assert result["chao1"] == 3
    assert result["ace"] == 3
    assert result["goods_coverage"] == 1
    assert result["d50"] == pytest.approx(1 / 3)
    assert result["simpson"] == pytest.approx(0.38)
    assert result["inverse_simpson"] == pytest.approx(1 / 0.38)
    assert result["gini_simpson"] == pytest.approx(0.62)
    assert result["berger_parker"] == pytest.approx(0.5)
    assert result["gini_coefficient"] == pytest.approx(0.2)
    assert result["clonality"] == pytest.approx(1 - result["norm_shannon_wiener"])
