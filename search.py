"""Search for Bs -> mu mu: fit of the mass spectrum and p-value."""

import numpy as np
from scipy import optimize, stats

M_B0 = 5.27966          # GeV, PDG
M_BS = 5.36692
BINS = np.linspace(4.9, 5.9, 51)        # 20 MeV bins


def histogram(m, bins=BINS):
    """Bin centres and counts of the masses m."""
    counts, edges = np.histogram(m, bins)
    return 0.5 * (edges[1:] + edges[:-1]), counts


def model(x, n_bkg, slope, n_b0, n_bs, sigma):
    """Expected counts per bin at masses x: exponential background, B0 and Bs peaks.

    n_b0 and n_bs are numbers of events; n_bkg is the background per bin at 5 GeV.
    """
    width = BINS[1] - BINS[0]
    peak = lambda mass: width * stats.norm.pdf(x, mass, sigma)
    return n_bkg * np.exp(-slope * (x - 5)) + n_b0 * peak(M_B0) + n_bs * peak(M_BS)

def fit(
        m: np.ndarray,
        sigma: float,
        with_bs: bool = True,
    ) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Fit the invariant mass distribution to the signal and background histogram model.
    
    If with_bs is True, the model includes the background, the B0 peak and
    the Bs peak. Otherwise, the Bs contribution is fixed to zero.

    Inputs
    ----------
    m : np.ndarray (the masses in GeV)
    sigma : float (Constant mass sigma)
    with_bs : bool, default=True
        Whether to include the Bs signal in the fit.
        If False, n_bs is fixed to zero.

    Returns
    -------
    params : np.ndarray
        Best-fit parameters. If with_bs is True, the order is
        [n_bkg, slope, n_b0, n_bs]. Otherwise, the order is
        [n_bkg, slope, n_b0].
    param_errors : np.ndarray
        One-standard-deviation uncertainties of the fitted parameters,
        estimated from the covariance matrix.
    chi2 : float
        Chi-square statistic of the fit.

    May raise RuntimeError if the nonlinear fit fails to converge.
    """
    x, counts = histogram(m)
    errors = np.sqrt(np.maximum(counts, 1))

    if with_bs:
        p0 = [max(float(np.mean(counts)), 1.0), 1.0, 10.0, 10.0]
        bounds = (
            [0.0, -np.inf, 0.0, 0.0],
            [np.inf, np.inf, np.inf, np.inf],
        )
        fit_function = lambda x, n_bkg, slope, n_b0, n_bs: model(
            x, n_bkg, slope, n_b0, n_bs, sigma
        )
    else:
        p0 = [max(float(np.mean(counts)), 1.0), 1.0, 10.0]
        bounds = (
            [0.0, -np.inf, 0.0],
            [np.inf, np.inf, np.inf],
        )
        fit_function = lambda x, n_bkg, slope, n_b0: model(
            x, n_bkg, slope, n_b0, 0.0, sigma
        )

    params, covariance = optimize.curve_fit(
        fit_function,
        x,
        counts,
        p0=p0,
        sigma=errors,
        absolute_sigma=True,
        bounds=bounds,
        maxfev=10000,
    )

    param_errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))

    expected = fit_function(x, *params)
    chi2 = float(np.sum(((counts - expected) / errors) ** 2))

    return params, param_errors, chi2


def pvalue(
        m: np.ndarray,
        sigma: float,
    ) -> tuple[float, float]:
    """
    Calculate the p-value and significance of the Bs signal.

    Parameters
    ----------
    m : np.ndarray
        One-dimensional array of invariant masses, in GeV.
    sigma : float
        Fixed mass resolution of the signal peaks, in GeV. Must be positive.

    Returns
    -------
    p : float (p-value)
    significance : float
        Significance expressed in standard deviations (sigma), calculated
        from the two-sided standard normal tail probability p.

    """
    _, _, chi2_0 = fit(m, sigma, with_bs=False)
    _, _, chi2_1 = fit(m, sigma, with_bs=True)

    delta_chi2 = max(chi2_0 - chi2_1, 0.0)
    p = float(stats.chi2.sf(delta_chi2, df=1))
    significance = float(np.sqrt(delta_chi2))

    return p, significance
