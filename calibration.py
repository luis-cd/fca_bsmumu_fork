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
        res = df.copy() #creamos una copia exacta del DataFrame
    
        #primero creamos una condición para saber que escala le toca a cada muon  (máscara)
    
        low1 = res['mu1_eta'] < ETA_SPLIT  #comprobamos la primera columna para ver fila por fila si es menor que el límite establecido
        low2 = res['mu2_eta'] < ETA_SPLIT  #lo mismo pero para el segundo muon
    
        res['mu1_pt'] = res['mu1_pt']/np.where(low1, self.low, self.high) #dividimos el momento trasversal, dividiéndolo por su factor de escala
        #si está por debajo, se cumple la condicion low1 y aplicamos la escala self.low
        #si está por encima o es igual, la condición no se cumple y utilizamos self.high como escala
    
        res['mu2_pt'] = res['mu2_pt']/np.where(low2, self.low, self.high) #lo mismo para el momento transversal del segundo muon
        return res

    def fit(self, control):
        """Fit the J/psi peak in events with both muons in the same region; set the scales.

        control: selected events of the control sample. Returns self.
        """

        #REGIÓN LOW (eta < 3.5):
        low_events = control[(control.mu1_eta < ETA_SPLIT) & (control.mu2_eta < ETA_SPLIT)]

        #REGIÓN HIGH (eta > 3.5):
        high_events = control[(control.mu1_eta > ETA_SPLIT) & (control.mu2_eta > ETA_SPLIT)]

        #Calculo masas invariantes:
        mass_low = dimuon.invariant_mass(low_events)
        mass_high = dimuon.invariant_mass(high_events)   

        #Ajustamos el pico de J/psi en cada región:
        mean_low, error_low, sigma_low = fit_peak(mass_low)
        mean_high, error_high, sigma_high = fit_peak(mass_high)

        #Calculamos los factores de calibración:
        self.low = mean_low / M_JPSI
        self.high = mean_high / M_JPSI

        return self

