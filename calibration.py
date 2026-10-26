"""Momentum calibration of the muons, from the J/psi peak."""

import json

import numpy as np
from scipy import optimize

import dimuon

M_JPSI = 3.0969         # GeV, PDG
ETA_SPLIT = 3.5         # boundary between the 'low' and 'high' eta regions


def fit_peak(m, lo=2.85, hi=3.25, bins=80):
    """Fit a Gaussian plus a straight line to the histogram of m in (lo, hi).

    Returns the mean, its error and the width of the Gaussian.
    """
    #raise NotImplementedError       # TODO (B)
    counts,edges=np.histogram(m,bins=bins,range=(lo,hi))
    x=(edges[:-1]+edges[1:])/2
    x0=(lo+hi)/2
    def model(x,amplitude,mean,width,background,slope):
        peak=amplitude*np.exp((-(x-mean)**2)/(2*width**2))
        return peak+background+slope*(x-x0)
    background0=np.median(counts)
    amplitude0=max(counts.max()-background0,1)
    mean0=x[np.argmax(counts)]
    p0=[amplitude0,mean0,0.015,background0,0]
    count_errors=np.sqrt(np.maximum(counts,1))
    params,covariance=optimize.curve_fit(model,x,counts,p0=p0,sigma=count_errors,absolute_sigma=True,
                                         bounds=(
                                             [0, lo, 1e-6, -np.inf, -np.inf],
                                             [np.inf, hi, hi - lo, np.inf, np.inf],
                                         ))
    mean=params[1]
    error=np.sqrt(covariance[1,1])
    width=params[2]
    return mean,error,width
class Calibration:
    """Momentum scale of the muons in two eta regions: measured pT = scale * true pT."""

    def __init__(self, low=1.0, high=1.0):
        self.low, self.high = low, high

    def __repr__(self):
        return f'Calibration(low={self.low:.4f}, high={self.high:.4f})'

    def save(self, path='calibration.json'):
        with open(path, 'w') as f:
            json.dump({'low': self.low, 'high': self.high}, f)

    @classmethod
    def load(cls, path='calibration.json'):
        with open(path) as f:
            return cls(**json.load(f))

    def apply(self, df):
        """Return a copy of df with mu1_pt and mu2_pt corrected (divided by their scale).

        The interface agreed by the two teams. For now a dummy: it returns the data unchanged.
        """
        return df.copy()            # TODO (B): divide each pT by the scale of its eta region

    def fit(self, control):
        """Fit the J/psi peak in events with both muons in the same region; set the scales.

        control: selected events of the control sample. Returns self.
        """
        raise NotImplementedError   # TODO (B)
