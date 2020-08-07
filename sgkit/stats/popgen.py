from typing import Hashable

from xarray import DataArray, Dataset

from .aggregation import count_alleles


def diversity(
    ds: Dataset, allele_counts: Hashable = "variant_allele_counts",
) -> DataArray:
    if allele_counts not in ds:
        ds[allele_counts] = count_alleles(ds)
    ac = ds[allele_counts]
    an = ac.sum(axis=1)
    n_pairs = an * (an - 1) / 2
    n_same = (ac * (ac - 1) / 2).sum(axis=1)
    n_diff = n_pairs - n_same
    # Let's ignore missing data and division by zero for now.
    pi = n_diff / n_pairs
    # Because we're not providing any arguments on windowing, etc,
    # we return the total over the whole region. Maybe this isn't
    # the behaviour we want, but it's a starting point. Note that
    # this is different to the tskit default behaviour where we
    # normalise by the size of windows so that results
    # in different windows are comparable. However, we don't have
    # any information about the overall length of the sequence here
    # so we can't normalise by it.
    return pi.sum()  # type: ignore[no-any-return]


def divergence(
    ds1: Dataset, ds2: Dataset, allele_counts: Hashable = "variant_allele_counts",
) -> DataArray:
    if allele_counts not in ds1:
        ds1[allele_counts] = count_alleles(ds1)
    ac1 = ds1[allele_counts]
    if allele_counts not in ds2:
        ds2[allele_counts] = count_alleles(ds2)
    ac2 = ds2[allele_counts]
    an1 = ds1[allele_counts].sum(axis=1)
    an2 = ds2[allele_counts].sum(axis=1)

    n_pairs = an1 * an2
    n_same = (ac1 * ac2).sum(axis=1)
    n_diff = n_pairs - n_same
    # Ignore missing data and division by zero for now.
    div = n_diff / n_pairs
    return div.sum()  # type: ignore[no-any-return]


def test(ac1: DataArray) -> DataArray:
    return ac1
