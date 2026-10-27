import numpy as np
import pytest

import search


def test_model_peak_holds_its_events():
    """With no background, the Bs peak spreads n_bs events over the bins."""
    x, _ = search.histogram([])
    counts = search.model(x, n_bkg=0, slope=1, n_b0=0, n_bs=100, sigma=0.022)
    assert counts.sum() == pytest.approx(100, rel=0.01)

def test_fit_returns_expected_shapes():
    """The fit returns parameters, uncertainties and a chi-square."""
    rng = np.random.default_rng(42)
    m = np.concatenate([
        rng.uniform(4.9, 5.9, 1000),
        rng.normal(search.M_B0, 0.022, 100),
        rng.normal(search.M_BS, 0.022, 50),
    ])

    params, errors, chi2 = search.fit(m, sigma=0.022)

    assert params.shape == (4,)
    assert errors.shape == (4,)
    assert np.all(np.isfinite(params))
    assert np.all(np.isfinite(errors))
    assert np.isfinite(chi2)
    assert np.all(errors >= 0)


def test_fit_without_bs_returns_three_parameters():
    """When Bs is excluded, the fit returns three parameters."""
    rng = np.random.default_rng(42)
    m = rng.uniform(4.9, 5.9, 1000)

    params, errors, chi2 = search.fit(m, sigma=0.022, with_bs=False)

    assert params.shape == (3,)
    assert errors.shape == (3,)
    assert np.all(np.isfinite(params))
    assert np.isfinite(chi2)


def test_pvalue_returns_valid_results():
    """The p-value and significance are finite and in valid ranges."""
    rng = np.random.default_rng(42)
    m = np.concatenate([
        rng.uniform(4.9, 5.9, 1000),
        rng.normal(search.M_B0, 0.022, 100),
        rng.normal(search.M_BS, 0.022, 50),
    ])

    p, significance = search.pvalue(m, sigma=0.022)

    assert 0.0 <= p <= 1.0
    assert significance >= 0.0
    assert np.isfinite(p)
    assert np.isfinite(significance)
