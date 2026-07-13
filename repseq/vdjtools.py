import warnings

from .io import save_to_vdjtools as _save_to_vdjtools


def save_to_vdjtools(*args, **kwargs):
    """
    Deprecated alias for `repseq.io.save_to_vdjtools`.

    Use `repseq.io.save_to_vdjtools` instead. This alias will be removed in a
    future version.
    """
    warnings.warn(
        "`repseq.vdjtools.save_to_vdjtools` is deprecated; use "
        "`repseq.io.save_to_vdjtools` instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _save_to_vdjtools(*args, **kwargs)
