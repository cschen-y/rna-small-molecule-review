






"""Common utility functions for various Bio submodules."""

import os
from typing import Any
from collections.abc import Callable
from typing import cast
from typing import Optional
from typing import Protocol
from typing import TypeVar



F = TypeVar("F", bound=Callable[..., object])


class _FunctionWithPrevious(Protocol[F]):
    previous: int | None
    __call__: F


def function_with_previous(func: F) -> _FunctionWithPrevious[F]:
    """Decorate a function as having an attribute named 'previous'."""
    function_with_previous = cast(_FunctionWithPrevious[F], func)
    
    function_with_previous.previous = None
    return function_with_previous


def find_test_dir(start_dir: str | None = None) -> str:
    """Find the absolute path of Biopython's Tests directory.

    Arguments:
    start_dir -- Initial directory to begin lookup (default to current dir)

    If the directory is not found up the filesystem's root directory, an
    exception will be raised.

    """
    if not start_dir:
        
        
        
        start_dir = "."

    target = os.path.abspath(start_dir)
    while True:
        if os.path.isdir(os.path.join(target, "Bio")) and os.path.isdir(
            os.path.join(target, "Tests")
        ):
            
            return os.path.abspath(os.path.join(target, "Tests"))
        
        
        new, tmp = os.path.split(target)
        if target == new:
            
            break
        target = new
    raise ValueError(
        f"Not within Biopython source tree: {os.path.abspath(start_dir)!r}"
    )


def run_doctest(target_dir: str | None = None, *args: Any, **kwargs: Any) -> None:
    """Run doctest for the importing module."""
    import doctest

    
    doctest_attributes: dict[str, Any] = {"optionflags": doctest.ELLIPSIS}
    doctest_attributes.update(kwargs)

    cur_dir = os.path.abspath(os.curdir)

    print("Running doctests...")
    try:
        os.chdir(find_test_dir(target_dir))
        doctest.testmod(*args, **doctest_attributes)
    finally:
        
        os.chdir(cur_dir)
    print("Done")


if __name__ == "__main__":
    run_doctest()
