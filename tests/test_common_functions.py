from repseq.common_functions import center_5
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


