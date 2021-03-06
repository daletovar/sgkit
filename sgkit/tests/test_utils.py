from typing import Any, List

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from sgkit.typing import ArrayLike
from sgkit.utils import check_array_like, encode_array, split_array_chunks


def test_check_array_like():
    with pytest.raises(TypeError):
        check_array_like("foo")
    a = np.arange(100, dtype="i4")
    with pytest.raises(TypeError):
        check_array_like(a, dtype="i8")
    with pytest.raises(TypeError):
        check_array_like(a, dtype={"i1", "i2"})
    with pytest.raises(TypeError):
        check_array_like(a, kind="f")
    with pytest.raises(ValueError):
        check_array_like(a, ndim=2)
    with pytest.raises(ValueError):
        check_array_like(a, ndim={2, 3})


def test_encode_array():
    def check(x: ArrayLike, values: ArrayLike, names: List[Any]) -> None:
        v, n = encode_array(x)
        np.testing.assert_equal(v, values)
        np.testing.assert_equal(n, names)

    check([], [], [])
    check(["a"], [0], ["a"])
    check(["a", "b"], [0, 1], ["a", "b"])
    check(["b", "a"], [0, 1], ["b", "a"])
    check(["a", "b", "b"], [0, 1, 1], ["a", "b"])
    check(["b", "b", "a"], [0, 0, 1], ["b", "a"])
    check(["b", "b", "a", "a"], [0, 0, 1, 1], ["b", "a"])
    check(["c", "a", "a", "b"], [0, 1, 1, 2], ["c", "a", "b"])
    check(["b", "b", "c", "c", "c", "a", "a"], [0, 0, 1, 1, 1, 2, 2], ["b", "c", "a"])
    check(["b", "c", "b", "c", "a"], [0, 1, 0, 1, 2], ["b", "c", "a"])
    check([2, 2, 1, 3, 1, 5, 5, 1], [0, 0, 1, 2, 1, 3, 3, 1], [2.0, 1.0, 3.0, 5.0])
    check(
        [2.0, 2.0, 1.0, 3.0, 1.0, 5.0, 5.0, 1.0],
        [0, 0, 1, 2, 1, 3, 3, 1],
        [2.0, 1.0, 3.0, 5.0],
    )


@pytest.mark.parametrize(  # type: ignore[misc]
    "n,blocks,expected_chunks",  # type: ignore[no-untyped-def]
    [
        (1, 1, [1]),
        (2, 1, [2]),
        (2, 2, [1] * 2),
        (3, 1, [3]),
        (3, 3, [1] * 3),
        (3, 2, [2, 1]),
        (7, 2, [4, 3]),
        (7, 3, [3, 2, 2]),
        (7, 7, [1] * 7),
    ],
)
def test_split_array_chunks__precomputed(
    n: int, blocks: int, expected_chunks: List[int]
):
    assert split_array_chunks(n, blocks) == tuple(expected_chunks)


@given(st.integers(1, 50), st.integers(0, 50))  # type: ignore[misc]
@settings(max_examples=50)  # type: ignore[misc]
def test_split_array_chunks__size(a: int, b: int) -> None:
    res = split_array_chunks(a + b, a)
    assert sum(res) == a + b
    assert len(res) == a


def test_split_array_chunks__raise_on_blocks_gt_n():
    with pytest.raises(
        ValueError,
        match=r"Number of blocks .* cannot be greater than number of elements",
    ):
        split_array_chunks(3, 10)


def test_split_array_chunks__raise_on_blocks_lte_0():
    with pytest.raises(ValueError, match=r"Number of blocks .* must be >= 0"):
        split_array_chunks(3, 0)


def test_split_array_chunks__raise_on_n_lte_0():
    with pytest.raises(ValueError, match=r"Number of elements .* must be >= 0"):
        split_array_chunks(0, 0)
