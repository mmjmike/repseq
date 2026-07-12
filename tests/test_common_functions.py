from repseq import common_functions
from repseq.common_functions import center_5, run_parallel_calculation
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

