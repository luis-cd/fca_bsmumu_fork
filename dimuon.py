"""Muon pairs: invariant mass and selection."""

from functools import reduce

import numpy as np

M_MU = 0.1056584        # GeV


def invariant_mass(df):

    # Muón 1: componentes del momento
    px1 = df["mu1_pt"] * np.cos(df["mu1_phi"])
    py1 = df["mu1_pt"] * np.sin(df["mu1_phi"])
    pz1 = df["mu1_pt"] * np.sinh(df["mu1_eta"])

    # Muón 2: componentes del momento
    px2 = df["mu2_pt"] * np.cos(df["mu2_phi"])
    py2 = df["mu2_pt"] * np.sin(df["mu2_phi"])
    pz2 = df["mu2_pt"] * np.sinh(df["mu2_eta"])

    # Energía de cada muón
    E1 = np.sqrt(px1**2 + py1**2 + pz1**2 + M_MU**2)
    E2 = np.sqrt(px2**2 + py2**2 + pz2**2 + M_MU**2)

    # Masa invariante del par
    m2 = (
        (E1 + E2)**2
        - (px1 + px2)**2
        - (py1 + py2)**2
        - (pz1 + pz2)**2
    )

    # Evitar pequeños negativos por errores de redondeo
    return np.sqrt(np.maximum(m2, 0))


# A cut is a function: DataFrame -> boolean Series (True = keep the row)
# the cuts (opposite charges, muon identification, vertex)
def has_trigger(df):
    return df['trigger']

def opposite_charges(df):
    return df["mu1_q"] * df["mu2_q"] < 0

#Parametro de la señal que se ve en la distribution
c_pid = 0.8
c_vtx = 5.0

def mu_pid(df):
    return(
        (df["mu1_pid"] > c_pid) 
        & (df["mu2_pid"] > c_pid)
    )

def vtx_chi(df):
    return df["vtx_chi2"] < c_vtx


# TODO (C): more cuts: opposite charges (mu1_q, mu2_q), and on mu1_pid, mu2_pid and vtx_chi2

CUTS = [has_trigger, opposite_charges, mu_pid, vtx_chi]


def select(df, cuts=CUTS):
    """Rows of df without missing values that pass all the cuts.
    Combine the cuts with functools.reduce.
    """
    df = df.dropna()

    mask = reduce(
        lambda a, b: a & b,
        (cut(df) for cut in cuts)
    )
    
    return df.loc[mask]       #TODO(C)
