import pytest
import msprime
import numpy as np

from sgkit import ts_to_dataset, diversity, divergence

# @pytest.mark.parametrize(
#    "ts",
#    [
#        msprime.simulate(...),
#        msprime.simulate(...),
#    ]
# )
# @pytest.mark.parametrize(
#    "windows",
#    [
#
#
#    ]
# )
def test_diversity(
    ts,
    # windows,
    # span_normalize
):
    ds = ts_to_dataset(ts)
    div = diversity(ds).compute()
    ts_div = ts.diversity(span_normalise=False)
    assert np.allclose(div, ts_div)


# @pytest.mark.parametrize(
#    "ts",
#    [
#        msprime.simulate(...),
#        msprime.simulate(...),
#    ]
# )
# @pytest.mark.parametrize("samples")
# @pytest.mark.parametrize(
#    "windows",
#    [
#
#
#    ]
# )
def test_divergence(ts):
    subset_1 = ts.samples()[: ts.num_samples // 2]
    subset_2 = ts.samples()[ts.num_samples // 2 :]
    ds1 = ts_to_dataset(ts, subset_1)
    ds2 = ts_to_dataset(ts, subset_2)
    ac1 = count_alleles(ds1).compute()
    ac2 = count_alleles(ds2).compute()

    div = divergence(ac1, ac2)
    ts_div = ts.divergence([subset_1, subset_2], span_normalise=False)
    assert np.allclose(div, ts_div)
