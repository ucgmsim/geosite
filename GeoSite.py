"""
GeoSite - Version 1.2
Authors: Felipe Kuncar (felipe.kuncar@pg.canterbury.ac.nz | fkw21)
         Christopher de la Torre (christopher.delatorre@canterbury.ac.nz | cde62)
Created on 16-09-2021
Last modification on 04-02-2022
Python Version: 3.8
"""

# =====================================================================================================================
# IMPORT LIBRARIES AND FUNCTIONS
# =====================================================================================================================

import numpy as np
from scipy.optimize import fsolve
import os.path
import pandas as pd
from math import sin, pi
import math
import os
import scipy.integrate as integrate
import matplotlib.pyplot as plt


# =====================================================================================================================
# PATHS
# =====================================================================================================================

# Define the path to the root directory
rootDir = 'C:/Users/Felipe Kuncar/Desktop/PhD/7. RESEARCH/6. Objectives/Obj1/3. Site Characterization/GeoSite/version 1/version 1.2'

# =====================================================================================================================
# DEFINE CONSTANTS
# =====================================================================================================================

# Atmospheric pressure (kPa)
pa = 101.325

# Unit weight of water (kN/m3)
gamma_w = 9.80665

# Gravity (m/s2)
g = 9.80665

# Poisson's ratio of soil
nu = 0.25

# =====================================================================================================================
# CLASS CPTu
# =====================================================================================================================

class CPTu(object):
    """
    A cone penetration test with pore pressure measurement (CPTu) object

    Parameters
    ----------
    depth (array_like):
        Depth - in (m)
    qc (array_like):
        Cone resistance - in (kPa)
    fs (array_like):
        Sleeve friction resistance - in (kPA)
    u2 (array_like):
        Pore pressure, measured just behind the cone - in (kPa)
    GWL (float (default = 0.00)):
        Groundwater level - in (m)
    a (float (default = 0.80)):
        Net area ratio for cone - dimensionless
    """

    def __init__(self, depth, qc, fs, u2, GWL=0.00, a=0.80):
        self.depth = np.array(depth)
        self.qc_orig = np.array(qc)
        self.qc = np.where(np.array(qc) > 0, np.array(qc), 0.01)
        self.fs_orig = np.array(fs)
        self.fs = np.where(np.array(fs) > 0, np.array(fs), 0.01)
        self.u2 = np.array(u2)
        self.GWL = GWL if GWL > 0 else -GWL
        self.a = a
        self._sigma_v0 = None
        self._u0 = None

    @property
    def qt(self):
        """ Corrected cone resistance (corrected for pore water effects) - in (kPa) """
        qt = self.qc + self.u2 * (1 - self.a)
        return qt

    @property
    def Rf(self):
        """ Friction ratio - in (%) """
        Rf = (self.fs / self.qt) * 100
        return Rf

    @property
    def u0(self):
        """ In-situ equilibrium pore pressure - in (%) """
        if self._u0 is None:
            u0 = np.zeros(len(self.depth))
            mask = self.depth < self.GWL
            u0[mask] = 0.0
            u0[~mask] = gamma_w * (self.depth[~mask] - self.GWL)
            self._u0 = u0
        return self._u0

    @property
    def gamma(self):
        """
        It estimates the soil total unit weight - in (kN/m3)

        According to Robertson & Cabal (2010)

        References
        ----------
        Robertson P.K., Cabal K.L. (2010). Estimating soil unit weight from CPT. 2nd International Symposium on Cone Penetration Test, CPT'10, Huntington Beach, CA, USA.
        """
        gamma = (0.27 * np.log10(self.Rf) + 0.36 * np.log10(self.qt / pa) + 1.236) * gamma_w
        # If the values of qc or fs are zero, negative or non-existent, then enforce gamma = 18 (kN/m3)
        gamma = np.where((self.qc == 0.01) | (self.fs == 0.01), 18, gamma)
        return gamma

    @property
    def sigma_v0(self):
        """ In-situ vertical total stress - in (kPa) """
        if self._sigma_v0 is None:
            sigma_v0 = np.zeros(len(self.depth))
            sigma_v0[0] = sigma_v0[0] + self.gamma[0] * self.depth[0]
            for i in range(1, len(self.depth)):
                sigma_v0[i] = sigma_v0[i - 1] + self.gamma[i] * (self.depth[i] - self.depth[i - 1])
            self._sigma_v0 = sigma_v0
        return self._sigma_v0

    @property
    def sigma_v0_eff(self):
        """ In-situ vertical effective stress - in (kPa) """
        sigma_v0_eff = self.sigma_v0 - self.u0
        return sigma_v0_eff

    @property
    def Qt(self):
        """ Normalized cone penetration resistance - dimensionless """
        Qt = (self.qt - self.sigma_v0) / self.sigma_v0_eff
        return Qt

    @property
    def Fr(self):
        """ Normalized friction ratio - in (%) """
        Fr = (self.fs / (self.qt - self.sigma_v0)) * 100
        return Fr

    def Ic(self, version=1):
        """
        Soil behaviour type index - dimensionless

        Parameters
        ----------
        version: str, optional (default = 1)
            If version = 1, it computes the normalized soil behaviour type index,
            according to Robertson (2009)
            If version = 2, it computes the non-normalized soil behaviour type index,
            according to Robertson (2010)

        """
        if version == 1:
            Ic = np.zeros(len(self.depth))
            Qtn = np.zeros(len(self.depth))
            for i in range(0, len(self.depth)):
                Fr = self.Fr[i]
                qt = self.qt[i]
                sigma_v0 = self.sigma_v0[i]
                sigma_v0_eff = self.sigma_v0_eff[i]
                def equations(p):
                    Ic_, Qtn_, n_ = p
                    eq1 = Ic_ - (((3.47 - np.log10(Qtn_)) ** 2 + (np.log10(Fr) + 1.22) ** 2) ** 0.5)
                    eq2 = Qtn_ - (((qt - sigma_v0) / pa) * ((pa / sigma_v0_eff) ** n_))
                    eq3 = n_ - min(0.381 * Ic_ + 0.05 * (sigma_v0_eff / pa) - 0.15, 1)
                    return (eq1, eq2, eq3)
                Ic_, Qtn_, n_ = fsolve(equations, [1.80, 100.00, 0.50])
                Ic[i] = Ic_
                Qtn[i] = Qtn_

        if version == 2:
            Ic = ((3.47 - np.log10(self.qt / pa)) ** 2 + (np.log10(self.Rf) + 1.22) ** 2) ** 0.5
            Qtn = None

        # Only consider Ic if qc and fs are realistic values (existent values greater than 0)
        Ic = np.where((self.qc_orig > 0) & (self.fs_orig > 0), Ic, float('NaN'))

        return Ic, Qtn

    def Vs(self, Ic, correlation=1):
        """
        It estimates the shear-wave velocity of the soil, Vs - in (m/s)

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If version = 1, it computes Vs according to McGann et al. (2015)
            If version = 2, it computes Vs according to CPT Guide 6th Edition (2015)

        References
        ----------
        McGann, C. R., Bradley, B. A., Taylor, M. L., Wotherspoon, L. M., & Cubrinovski, M. (2015). Development of an
        empirical correlation for predicting shear wave velocity of Christchurch soils from cone penetration test data.
        Soil Dynamics and Earthquake Engineering, 75, 66–75.

        Robertson, P.K., Cabal K,L. (2015). Guide to cone penetration testing for geotechnical engineering, 6th Edition,
        GREGG.
        """
        if correlation == 1:
            Vs = 18.4 * (self.qt ** 0.144) * (self.fs ** 0.0823) * (self.depth ** 0.278)

        elif correlation == 2:
            alpha = 10 ** (0.55 * Ic + 1.68)
            Vs = (alpha * (self.qt - self.sigma_v0) / pa) ** 0.5

        Vs = np.where((self.qc_orig < 0) | (self.fs_orig < 0), float('NaN'), Vs)

        return Vs

    def relative_density(self, Ic, Qtn, correlation=1):
        """
        It estimates the relative density of the soil, Dr - in (%)

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If correlation = 1, it computes Dr according to Idriss & Boulanger (2008)
            If correlation = 2, it computes Dr according to Kulhawy & Mayne (1990) - simplified expression from CPT
            Guide 6th Edition (2015)

        References
        ----------
        Idriss, I. M., and Boulanger, R. W. (2008) "Soil Liquefaction during Earthquakes." MNO-12, Earthquake
        Engineering Research Institute, Oakland, CA.

        Kulhawy, F.H., Mayne, P.H., 1990. Manual on estimating soil properties for foundation design, Report EL-6800
        Power Research Institute, EPRI, August 1990.

        Robertson, P.K., Cabal K,L. (2015). Guide to cone penetration testing for geotechnical engineering, 6th Edition,
        GREGG.
        """
        if correlation == 1:
            relative_density = np.zeros(len(self.depth))
            for i in range(0, len(self.depth)):
                qt = self.qt[i]
                sigma_v0_eff = self.sigma_v0_eff[i]
                def equations(p):
                    Dr, Cn, m = p
                    eq1 = Dr - (0.478 * ((Cn * qt / pa) ** 0.264) - 1.063)
                    eq2 = Cn - (min((pa / sigma_v0_eff) ** m, 1.70))
                    eq3 = m - (0.784 - (0.521 * Dr))
                    return (eq1, eq2, eq3)
                Dr, Cn, m_ = fsolve(equations, [0.50, 1.50, 0.52])
                if Ic[i] <= 2.6:
                    relative_density[i] = min(Dr * 100, 100)
                else:
                    relative_density[i] = None

        if correlation == 2:
            relative_density = np.where(Ic <= 2.6, np.where(np.sqrt(Qtn / 350) <= 1.00, (np.sqrt(Qtn / 350)), 1.00) * 100, float('NaN'))

        # Only consider Dr if qc and fs are realistic values (existent values greater than 0)
        relative_density = np.where((self.qc_orig > 0) & (self.fs_orig > 0), relative_density, float('NaN'))

        return relative_density

    def friction_angle(self, Ic, Qtn, correlation=1):
        """
        It estimates the friction angle of the soil, Dr - in (°)

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If correlation = 1, it computes the peak friction angle according to Robertson (2010)
            If correlation = 2, it computes the peak friction angle according to Kulhawy & Mayne (1990)
            If correlation = 3, it computes the peak friction angle according to Robertson & Campanella (1983)

        References
        ----------
        Kulhawy, F.H., Mayne, P.H. (1990). Manual on estimating soil properties for foundation design, Report EL-6800
        Power Research Institute, EPRI, August 1990.

        Robertson, P.K., Campanella, R.G. (1983). Interpretation of cone penetration tests – Part I (sand). Canadian
        Geotechnical Journal, 20(4):718-733.

        """
        if correlation == 1:
            # correction factor
            Kc = np.where(Ic <= 1.64, 1.00, 5.581 * (Ic ** 3) - 0.403 * (Ic ** 4) - 21.63 * (Ic ** 2) + 33.75 * Ic - 17.88)
            # equivalent clean sand normalized cone resistance
            Qtncs = Kc * Qtn
            # constant volume (or critical state) friction angle
            # typically about 33 (°) for quartz sands but can be as high as 40 (°) for felspathic sand
            phi_cv = 33 * np.ones(len(self.depth))
            # peak friction angle
            friction_angle = np.where(Ic < 2.60, phi_cv + 15.84 * np.log10(Qtncs) - 26.88 * np.ones(len(self.depth)), float('NaN'))

        elif correlation == 2:
            friction_angle = np.where(Ic < 2.60, 17.6 + 11 * np.log10(Qtn), float('NaN'))

        elif correlation == 3:
            friction_angle = np.where(Ic < 2.6, np.arctan((1/2.68) * (np.log10(self.qt / self.sigma_v0_eff) + 0.29)) * 360 / (2 * np.pi), float('NaN'))

        # Only consider friction_angle if qc and fs are realistic values (existent values greater than 0)
        friction_angle = np.where((self.qc_orig > 0) & (self.fs_orig > 0), friction_angle, float('NaN'))

        return friction_angle

# =====================================================================================================================
# CLASS Vs
# =====================================================================================================================

class Vs(object):
    """
    A shear wave velocity (Vs) profile object

    Parameters
    ----------
    depth (array_like):
        Depth - in (m)
    Vs_values (array_like):
        Cone resistance - in (kPa)
    GWL (float (default = 0.00)):
        Groundwater level - in (m)
    """

    def __init__(self, depth, values, GWL=0.00):
        self.depth = np.array(depth)
        self.values = np.array(values)
        self.GWL = GWL if GWL > 0 else -GWL
        self._sigma_v0 = None
        self._u0 = None

    @property
    def gamma(self):
        """
        It estimates the soil total unit weight - in (kN/m3)

        According to NCHRP (2019)

        References
        ----------
        NCHRP (2019). Manual on subsurface investigations. Web-Only Document 258, NCHRP Project 21-10.

        """
        gamma = 7.83 * np.log10(self.values) - 0.125 * np.ones(len(self.depth))

        return gamma

    def gamma_corr(self, correlation=1):
        """
        Estimate the soil total unit weight - in (kN/m3)

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If correlation = 1, it computes the soil total unit weight according to NCHRP (2019).
            If correlation = 2, it computes the soil total unit weight according to Boore (2015).

        References
        ----------
        NCHRP (2019). Manual on subsurface investigations. Web-Only Document 258, NCHRP Project 21-10.

        Boore D. (2015). Notes on relating density to velocity for use in site amplification calculations.

        """
        if correlation == 1:
            gamma_corr = 7.83 * np.log10(self.values) - 0.125 * np.ones(len(self.depth))

        elif correlation == 2:
            # conversion of Vs from (m/s) to (km/s)
            Vs = self.values / 1000
            # density in (g/cm3) - equation valid for Vs < 3550 (m/s)
            rho = np.where(Vs < 0.30, 1 + ((1.53 * (Vs ** 0.85))/(0.35 + 1.889 * (Vs ** 1.7))), 1.74 * (0.9409 + 2.0947 * Vs - 0.8206 * (Vs ** 2) + 0.2683 * (Vs ** 3) - 0.0251 * (Vs ** 4)) ** 0.25)
            # soil total unit weight
            gamma_corr = rho * 9.81

        return gamma_corr

    @property
    def u0(self):
        """ In-situ equilibrium pore pressure - in (%) """
        if self._u0 is None:
            u0 = np.zeros(len(self.depth))
            mask = self.depth < self.GWL
            u0[mask] = 0.0
            u0[~mask] = gamma_w * (self.depth[~mask] - self.GWL)
            self._u0 = u0
        return self._u0

    @property
    def sigma_v0(self):
        """ In-situ vertical total stress - in (kPa) """
        if self._sigma_v0 is None:
            sigma_v0 = np.zeros(len(self.depth))
            sigma_v0[0] = sigma_v0[0] + self.gamma[0] * self.depth[0]
            for i in range(1, len(self.depth)):
                sigma_v0[i] = sigma_v0[i - 1] + self.gamma[i] * (self.depth[i] - self.depth[i - 1])
            self._sigma_v0 = sigma_v0
        return self._sigma_v0

    @property
    def sigma_v0_eff(self):
        """ In-situ vertical effective stress - in (kPa) """
        sigma_v0_eff = self.sigma_v0 - self.u0
        return sigma_v0_eff

    def N1_60(self, correlation=1):
        """
        Corrected SPT N-Values (N1)60 - dimensionless

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If correlation = 1, it estimates the N60 according to NCHRP (2019).

        References
        ----------
        NCHRP (2019). Manual on subsurface investigations. Web-Only Document 258, NCHRP Project 21-10.
        """

        if correlation == 1:
            N_60 = (self.values/97) ** (1 / 0.314)
            C_N = np.where((pa/self.sigma_v0_eff) ** 0.5 <= 2, (pa/self.sigma_v0_eff) ** 0.5, 2)
            N1_60 = C_N * N_60

        return N1_60

    def relative_density(self, N1_60, correlation=1):
        """
        It estimates the relative density of the soil, Dr - in (%)

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If correlation = 1, it computes Dr according to Idriss & Boulanger (2008)
            If correlation = 2, it computes Dr according to Kulhawy & Mayne (1990)

        References
        ----------
        Idriss, I. M., and Boulanger, R. W. (2008) "Soil Liquefaction during Earthquakes." MNO-12, Earthquake
        Engineering Research Institute, Oakland, CA.

        Kulhawy, F.H., Mayne, P.H., 1990. Manual on estimating soil properties for foundation design, Report EL-6800
        Power Research Institute, EPRI, August 1990.

        """
        if correlation == 1:
            relative_density = np.sqrt(N1_60 / 46) * 100

        return relative_density

    def friction_angle(self, N1_60, correlation=1):
        """
        It estimates the friction angle of the soil - in (°)

        Parameters
        ----------
        correlation: str, optional (default = 1)
            If correlation = 1, it estimates the friction angle according to NCHRP (2019).

        References
        ----------
        NCHRP (2019). Manual on subsurface investigations. Web-Only Document 258, NCHRP Project 21-10.
        """

        if correlation == 1:
            friction_angle = 20 + np.sqrt(15.4 * N1_60)

        return friction_angle


# =====================================================================================================================
# FUNCTIONS FOR GENERATING TRIAL (OR INITIAL) MODELS
# =====================================================================================================================

def trialModelsFromCPTu(layersThickness, CPTu):

    # Obtain Ic
    Ic, Qtn = CPTu.Ic(version=1)

    # Obtain layersDepth from layersThickness
    layersDepth = np.zeros(1 + len(layersThickness))
    for i in range(1, len(layersDepth)):
        layersDepth[i] = layersDepth[i - 1] + layersThickness[i - 1]
    print(layersDepth)

    # Obtain counter that determine data within each layer
    j = 1
    counterLim = [0]
    for i in range(0, len(CPTu.depth)):
        if CPTu.depth[i] > layersDepth[j]:
            counterLim.append(i)
            j = j + 1

    layersThickness_trial = layersThickness[0:len(counterLim) - 1]

    # Initialize arrays
    Ic_trial = []
    Vs_trial1 = []
    Vs1 = CPTu.Vs(Ic, correlation=1)
    Vs_trial2 = []
    Vs2 = CPTu.Vs(Ic, correlation=2)
    gamma_trial = []
    gamma = CPTu.gamma
    Dr_trial1 = []
    Dr1 = CPTu.relative_density(Ic, Qtn, correlation=1)
    Dr_trial2 = []
    Dr2 = CPTu.relative_density(Ic, Qtn, correlation=2)
    phi_trial1 = []
    phi1 = CPTu.friction_angle(Ic, Qtn, correlation=1)
    phi_trial2 = []
    phi2 = CPTu.friction_angle(Ic, Qtn, correlation=2)
    phi_trial3 = []
    phi3 = CPTu.friction_angle(Ic, Qtn, correlation=3)

    for i in range(1, len(layersThickness_trial) + 1):
        Ic_trial.append(np.nanmean(Ic[counterLim[i-1]:counterLim[i]]))
        Vs_trial1.append(np.nanmean(Vs1[counterLim[i-1]:counterLim[i]]))
        Vs_trial2.append(np.nanmean(Vs2[counterLim[i - 1]:counterLim[i]]))
        gamma_trial.append(np.nanmean(gamma[counterLim[i - 1]:counterLim[i]]))
        Dr_trial1.append(np.nanmean(Dr1[counterLim[i - 1]:counterLim[i]]))
        Dr_trial2.append(np.nanmean(Dr2[counterLim[i - 1]:counterLim[i]]))
        phi_trial1.append(np.nanmean(phi1[counterLim[i - 1]:counterLim[i]]))
        phi_trial2.append(np.nanmean(phi2[counterLim[i - 1]:counterLim[i]]))
        phi_trial3.append(np.nanmean(phi3[counterLim[i - 1]:counterLim[i]]))


    return layersThickness_trial, Ic_trial, Vs_trial1, Vs_trial2, gamma_trial, Dr_trial1, Dr_trial2, phi_trial1, phi_trial2, phi_trial3


'''
def trialModelsFromVs(layersThickness, CPTu):
'''

# =====================================================================================================================
# CLASS PDMY02
# =====================================================================================================================

class PDMY02(object):
    """
    A PDMY02 object that allows the estimation of model-specific parameters from the relative density values, according
    to the calibration of Karimi & Dashti (2015) for Nevada sand

    Parameters
    ----------
    depth (array_like):
        Depth - in (m)
    Dr (array_like):
        Relative density - in (%)

    References
    ----------
    Karimi, Z., & Dashti, S. (2016). Numerical and Centrifuge Modeling of Seismic Soil–Foundation–Structure Interaction
    on Liquefiable Ground. In Journal of Geotechnical and Geoenvironmental Engineering (Vol. 142, Issue 1, p. 04015061).
    American Society of Civil Engineers (ASCE).
    """

    def __init__(self, depth, Dr):
        self.depth = np.array(depth)
        self.Dr = np.array(Dr)

    def PTAng(self):
        """
        Estimate the phase transformation angle - in (°)

        """
        PTAng = np.where(self.Dr < 30.0,  31.0,
                np.where(self.Dr < 40.0,  31.0 + -0.1000 * (self.Dr - 30.0),
                np.where(self.Dr < 50.0,  30.0 + -0.4500 * (self.Dr - 40.0),
                np.where(self.Dr < 63.0,  25.5 +  0.0769 * (self.Dr - 50.0),
                np.where(self.Dr < 68.0,  26.5 + -0.1000 * (self.Dr - 63.0),
                np.where(self.Dr < 90.0,  26.0 +  0.0227 * (self.Dr - 68.0), 26.5))))))

        return PTAng

    def e(self):
        """
        It estimates the void ratio - dimensionless

        """
        e = np.where(self.Dr < 30.0,  0.76,
            np.where(self.Dr < 40.0,  0.76 + -0.0030 * (self.Dr - 30.0),
            np.where(self.Dr < 50.0,  0.73 + -0.0030 * (self.Dr - 40.0),
            np.where(self.Dr < 63.0,  0.70 + -0.0031 * (self.Dr - 50.0),
            np.where(self.Dr < 68.0,  0.66 + -0.0020 * (self.Dr - 63.0),
            np.where(self.Dr < 90.0,  0.65 + -0.0032 * (self.Dr - 68.0), 0.58))))))

        return e

    def contrac1(self):
        """
        It estimates the parameter contrac1

        """
        contrac1 = np.where(self.Dr < 30.0,  0.087,
                   np.where(self.Dr < 40.0,  0.087 + -0.0020 * (self.Dr - 30.0),
                   np.where(self.Dr < 50.0,  0.067 + -0.0017 * (self.Dr - 40.0),
                   np.where(self.Dr < 63.0,  0.050 + -0.0008 * (self.Dr - 50.0),
                   np.where(self.Dr < 68.0,  0.040 + -0.0040 * (self.Dr - 63.0),
                   np.where(self.Dr < 90.0,  0.020 + -0.0002 * (self.Dr - 68.0), 0.016))))))

        return contrac1

    def contrac2(self):
        """
        It estimates the parameter contrac2

        """
        contrac2 = np.where(self.Dr < 30.0,  5.00,
                   np.where(self.Dr < 40.0,  5.00 + -0.0500 * (self.Dr - 30.0),
                   np.where(self.Dr < 50.0,  4.50 + -0.0500 * (self.Dr - 40.0),
                   np.where(self.Dr < 63.0,  4.50 + -0.1154 * (self.Dr - 50.0),
                   np.where(self.Dr < 68.0,  2.50 + -0.2000 * (self.Dr - 63.0),
                   np.where(self.Dr < 90.0,  1.50 + -0.0023 * (self.Dr - 68.0), 1.40))))))

        return contrac2

    def contrac3(self):
        """
        It estimates the parameter contrac3

        """
        contrac3 = np.where(self.Dr < 30.0,  0.30,
                   np.where(self.Dr < 40.0,  0.30 + -0.0030 * (self.Dr - 30.0),
                   np.where(self.Dr < 50.0,  0.27 + -0.0020 * (self.Dr - 40.0),
                   np.where(self.Dr < 63.0,  0.25 + -0.0038 * (self.Dr - 50.0),
                   np.where(self.Dr < 68.0,  0.20 + -0.0100 * (self.Dr - 63.0),
                   np.where(self.Dr < 90.0,  0.15 + -0.0005 * (self.Dr - 68.0), 0.14))))))

        return contrac3

    def dilat1(self):
        """
        It estimates the parameter dilat1

        """
        dilat1 = np.where(self.Dr < 30.0,  0.01,
                 np.where(self.Dr < 40.0,  0.01 + 0.0010 * (self.Dr - 30.0),
                 np.where(self.Dr < 50.0,  0.02 + 0.0040 * (self.Dr - 40.0),
                 np.where(self.Dr < 63.0,  0.06 + 0.0008 * (self.Dr - 50.0),
                 np.where(self.Dr < 68.0,  0.07 + 0.0160 * (self.Dr - 63.0),
                 np.where(self.Dr < 90.0,  0.15 + 0.0045 * (self.Dr - 68.0), 0.25))))))

        return dilat1

    def dilat2(self):
        """
        It estimates the parameter dilat2

        """
        dilat2 = 3 * np.ones(len(self.Dr))

        return dilat2

    def dilat3(self):
        """
        It estimates the parameter dilat3

        """
        dilat3 = np.zeros(len(self.Dr))

        return dilat3

# =====================================================================================================================
# FUNCTION create_inputOpenSees
# =====================================================================================================================

def create_inputOpenSees(siteID, PathOutputs, layerThick, GWL, gamma, Vs, phi, Dr, refDepth, pressCoeff):

    GWL = GWL if GWL > 0 else -GWL

    # reverse order so that idx = 0 is deepest layer (i.e., model base)
    layerThick = layerThick[::-1]
    rho = gamma[::-1] / 9.81
    Vs = Vs[::-1]
    phi = phi[::-1]
    Dr = Dr[::-1]
    refDepth = refDepth[::-1]
    pressCoeff = pressCoeff[::-1]

    # define PDMY02 object
    PDMY02model = PDMY02(None, Dr)

    # define some variables to be used
    numLayers = len(layerThick) - 1
    soilThick = np.cumsum(layerThick)[-1]

    # open the .tcl file and write
    inputFile = open(os.path.join(PathOutputs, '%s_input.tcl' % siteID), 'w')

    inputFile.write('# define gravity and pi \n')
    inputFile.write('set g 9.80665 \n')
    inputFile.write('set pi [expr atan(1)*4] \n')

    inputFile.write('\n')
    inputFile.write('#---SOIL GEOMETRY \n')
    inputFile.write('\n')

    inputFile.write('# thicknesses of soil profile (m) \n')
    inputFile.write('set soilThick %s \n' % soilThick)
    inputFile.write('\n')

    inputFile.write('# number of soil layers (not counting layer 0 which is bedrock) \n')
    inputFile.write('set numLayers %s \n' % numLayers)
    inputFile.write('\n')

    inputFile.write('# layer thicknesses \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set layerThick(%s) %s \n' % (ele, layerThick[ele]))
    inputFile.write('\n')

    inputFile.write('# depth of water table (create a new layer at WT) \n')
    inputFile.write('# if water not present set waterTable anywhere below depth of model \n')
    inputFile.write('set waterTable %s \n' % GWL)
    inputFile.write('\n')

    inputFile.write('#---BASIC MATERIAL PROPERTIES \n')
    inputFile.write('\n')

    inputFile.write('# soil mass density (Mg/m^3) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set rho(%s) %s \n' % (ele, rho[ele]))
    inputFile.write('set rhoWater 1.0 \n')

    inputFile.write('\n')
    inputFile.write('# soil shear wave velocity for each layer (m/s) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set Vs(%s) %s \n' % (ele, Vs[ele]))
    inputFile.write('\n')

    inputFile.write('\n')
    inputFile.write('# Poisson ratio of soil \n')
    inputFile.write('set nu 0.25 \n')
    inputFile.write('\n')

    inputFile.write('# rock elastic properties) \n')
    inputFile.write('# bedrock shear wave velocity (m/s) \n')
    inputFile.write('set rockVS $Vs(0) \n')
    inputFile.write('# bedrock mass density (Mg/m^3) \n')
    inputFile.write('set rockDen $rho(0) \n')
    inputFile.write('\n')

    inputFile.write('\n')
    inputFile.write('#---OTHER PDMY02 MODEL PARAMETERS \n')
    inputFile.write('\n')

    inputFile.write('# soil friction angle (°) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set phi(%s) %s \n' % (ele, phi[ele]))
    inputFile.write('\n')

    inputFile.write('# peak shear strain \n')
    inputFile.write('set gammaPeak 0.1 \n')
    inputFile.write('\n')

    inputFile.write('# relative density (%) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('# Layer %s: Dr = %s \n' % (ele, Dr[ele]))
    inputFile.write('\n')

    inputFile.write('# reference pressure is computed as mean confining pressure \n')
    inputFile.write('# at refDepth for each layer (0 is ToL, 1 is BoL) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set refDepth(%s) %s \n' % (ele, refDepth[ele]))
    inputFile.write('\n')

    inputFile.write('# pressure dependency coefficient \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set pressCoeff(%s) %s \n' % (ele, pressCoeff[ele]))
    inputFile.write('\n')

    inputFile.write('# phase transformation angle (not for layer 0) (°) \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set phaseAng(%s) %s \n' % (ele, PDMY02model.PTAng()[ele]))
    inputFile.write('\n')

    inputFile.write('# contraction (not for layer 0) \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set contract1(%s) %s \n' % (ele, PDMY02model.contrac1()[ele]))
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set contract2(%s) %s \n' % (ele, PDMY02model.contrac2()[ele]))
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set contract3(%s) %s \n' % (ele, PDMY02model.contrac3()[ele]))
    inputFile.write('\n')

    inputFile.write('# dilation coefficients (not for layer 0) \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set dilate1(%s) %s \n' % (ele, PDMY02model.dilat1()[ele]))
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set dilate2(%s) %s \n' % (ele, PDMY02model.dilat2()[ele]))
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set dilate3(%s) %s \n' % (ele, PDMY02model.dilat3()[ele]))
    inputFile.write('\n')

    inputFile.write('# void ratio (need it for layer 0 for element definition) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set voidR(%s) %s \n' % (ele, PDMY02model.e()[ele]))
    inputFile.write('\n')

    return

# =====================================================================================================================
# FUNCTION create_modelOpenSees
# =====================================================================================================================

def create_modelOpenSees(siteID, PathOutputs, layerThick, GWL, gamma, LR, Vs, constModelFlag, phi, Dr, refDepth, pressCoeff):

    GWL = GWL if GWL > 0 else -GWL

    # Reverse order so that idx = 0 is deepest layer (i.e., model base)
    layerThick = layerThick[::-1]
    rho = gamma[::-1] / 9.81
    Vs = Vs[::-1]
    phi = phi[::-1]
    Dr = Dr[::-1]
    refDepth = refDepth[::-1]
    pressCoeff = pressCoeff[::-1]

    # Define PDMY02 object
    PDMY02model = PDMY02(None, Dr)

    # Define some variables to be used
    numLayers = len(layerThick) - 1
    soilThick = np.cumsum(layerThick)[-1]

    # Define reference Vs
    Vs_ref = Vs[0]

    # Open the .tcl file and write
    inputFile = open(os.path.join(PathOutputs, '%s_TotalStress_LR%.1f_%.0f.tcl' % (siteID, LR, Vs_ref)), 'w')

    inputFile.write('########################################################### \n')
    inputFile.write('#                                                         # \n')
    inputFile.write('# Total Stress Site Response Analysis                     # \n')
    inputFile.write('# Constitutive Models: PDMY02 (sands) & PIMY (clays)      # \n')
    inputFile.write('# Pore Pressure Generation Allowed: No                    # \n')
    inputFile.write('#                                                         # \n')
    inputFile.write('########################################################### \n')
    inputFile.write('\n')

    inputFile.write('# Extract exterior inputs and give them a variable name \n')
    inputFile.write('set site [lindex $argv 0] \n')
    inputFile.write('set modelID [lindex $argv 1] \n')
    inputFile.write('set gMotionPath [lindex $argv 2] \n')
    inputFile.write('set gMotionName [lindex $argv 3] \n')
    inputFile.write('set npts [lindex $argv 4] \n')
    inputFile.write('set dt [lindex $argv 5] \n')
    inputFile.write('set saveDir [lindex $argv 6] \n')
    inputFile.write('\n')

    inputFile.write('set dash "_" \n')
    inputFile.write('set analID $site$dash$modelID$dash$gMotionName\n')
    inputFile.write('set siteModelID $site$dash$modelID\n')
    inputFile.write('\n')

    inputFile.write('# Define gravity and pi \n')
    inputFile.write('set g 9.80665 \n')
    inputFile.write('set pi [expr atan(1)*4] \n')
    inputFile.write('\n')

    inputFile.write('# ----------------------------------------------------------------------------------------- \n')
    inputFile.write('#  1. DEFINE SOIL GEOMETERY AND MATERIAL PARAMETERS \n')
    inputFile.write('# ----------------------------------------------------------------------------------------- \n')

    inputFile.write('\n')
    inputFile.write('#---SOIL GEOMETRY \n')
    inputFile.write('\n')

    inputFile.write('# Thicknesses of soil profile (m) \n')
    inputFile.write('set soilThick %s \n' % soilThick)
    inputFile.write('\n')

    inputFile.write('# Number of soil layers (not counting layer 0 which is bedrock) \n')
    inputFile.write('set numLayers %s \n' % numLayers)
    inputFile.write('\n')

    inputFile.write('# Layer thicknesses \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set layerThick(%s) %s \n' % (ele, layerThick[ele]))
    inputFile.write('\n')

    inputFile.write('# Depth of water table (create a new layer at WT) \n')
    inputFile.write('# If water not present set waterTable anywhere below depth of model \n')
    inputFile.write('set waterTable %s \n' % GWL)
    inputFile.write('\n')

    inputFile.write('#---BASIC MATERIAL PROPERTIES \n')
    inputFile.write('\n')

    inputFile.write('# Soil mass density (Mg/m^3) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set rho(%s) %s \n' % (ele, rho[ele]))
    inputFile.write('set rhoWater 1.0 \n')

    inputFile.write('\n')
    inputFile.write('# Soil shear wave velocity for each layer (m/s) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set Vs(%s) %s \n' % (ele, Vs[ele]))
    inputFile.write('\n')

    inputFile.write('# Poisson ratio of soil \n')
    inputFile.write('set nu 0.25 \n')
    inputFile.write('\n')

    inputFile.write('# Rock elastic properties) \n')
    inputFile.write('# Bedrock shear wave velocity (m/s) \n')
    inputFile.write('set rockVS $Vs(0) \n')
    inputFile.write('# Bedrock mass density (Mg/m^3) \n')
    inputFile.write('set rockDen $rho(0) \n')
    inputFile.write('\n')

    inputFile.write('#--- MODEL PARAMETERS \n')
    inputFile.write('\n')

    inputFile.write('# Consitutive model \n')
    inputFile.write('# constModelFlag: 1 for PDMY02 (sands), 0 for PIDMY (clays) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set constModelFlag(%s) %s \n' % (ele, constModelFlag[ele]))
    inputFile.write('\n')

    inputFile.write('# Soil friction angle (°) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set phi(%s) %s \n' % (ele, phi[ele]))
    inputFile.write('\n')

    inputFile.write('# Peak shear strain \n')
    inputFile.write('set gammaPeak 0.1 \n')
    inputFile.write('\n')

    inputFile.write('# Relative density (%) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('# Layer %s: Dr = %s \n' % (ele, Dr[ele]))
    inputFile.write('\n')

    inputFile.write('# Reference pressure is computed as mean confining pressure \n')
    inputFile.write('# at refDepth for each layer (0 is ToL, 1 is BoL, 0.5 is in the middle) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set refDepth(%s) %s \n' % (ele, refDepth[ele]))
    inputFile.write('\n')

    inputFile.write('# Pressure dependency coefficient \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set pressCoeff(%s) %s \n' % (ele, pressCoeff[ele]))
    inputFile.write('\n')

    inputFile.write('# Phase transformation angle (not for layer 0) (°) \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set phaseAng(%s) %s \n' % (ele, PDMY02model.PTAng()[ele]))
    inputFile.write('\n')

    inputFile.write('# Contraction coefficients (not for layer 0) \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set contract1(%s) %s \n' % (ele, PDMY02model.contrac1()[ele]))
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set contract2(%s) %s \n' % (ele, PDMY02model.contrac2()[ele]))
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set contract3(%s) %s \n' % (ele, PDMY02model.contrac3()[ele]))
    inputFile.write('\n')

    inputFile.write('# Dilation coefficients (not for layer 0) \n')
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set dilate1(%s) %s \n' % (ele, PDMY02model.dilat1()[ele]))
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set dilate2(%s) %s \n' % (ele, PDMY02model.dilat2()[ele]))
    for ele in reversed(np.arange(1, len(layerThick))):
        inputFile.write('set dilate3(%s) %s \n' % (ele, PDMY02model.dilat3()[ele]))
    inputFile.write('\n')

    inputFile.write('# Void ratio (need it for layer 0 for element definition) \n')
    for ele in reversed(np.arange(0, len(layerThick))):
        inputFile.write('set voidR(%s) %s \n' % (ele, PDMY02model.e()[ele]))
    inputFile.write('\n')

    # Close input file
    inputFile.close()

    # Open and copy the lines written in the OpenSees model template
    OpenSeesTemplate = open(os.path.join(rootDir, 'OpenSeesTemplates', 'totalStressModel.tcl'))
    addLines = OpenSeesTemplate.read()
    OpenSeesTemplate.close()

    # Path to the input file
    inputFilePath = os.path.join(PathOutputs, '%s_TotalStress_LR%.1f_%.0f.tcl' % (siteID, LR, Vs_ref))

    # Add the copied lines to the input file
    with open(inputFilePath, 'a') as finalFile:
        finalFile.write(addLines)

    # Close final file
    finalFile.close()

    return

# =====================================================================================================================
# FUNCTION readVsData
# =====================================================================================================================

def readVsData(PathVsSW):

    Vs = pd.read_csv(PathVsSW, skiprows=0)
    LayerThickness = Vs.iloc[:, 1]
    VsValue = Vs.iloc[:, 2]

    depth = np.zeros(1 + len(LayerThickness))
    for i in range(1, len(depth)):
        depth[i] = depth[i-1] + LayerThickness[i-1]
    depth = np.repeat(depth, 2)
    depth = np.delete(depth, [0, len(depth) - 1])
    Vs = np.repeat(VsValue, 2)
    Vs = Vs.to_numpy()

    return depth, Vs

# =====================================================================================================================
# FUNCTION PressDepVs
# =====================================================================================================================

# Created by Chris de la Torre
def PressDepVs(PathOutputs, thick, Vs, gamma, phi, siteID, pressDependCoe, refDepth, VsInversionTop, waterDepth):

    numLayers = len(thick)

    waterDepth = waterDepth if waterDepth > 0 else -waterDepth

    # ------------------------------------------------------------------------------------------------------------------
    # Define the size of the arrays to be generated containing the reference parameters
    # ------------------------------------------------------------------------------------------------------------------

    # Define the size of the arrays to be generated
    # Size = number of layers
    vertStress = np.zeros(numLayers)  # 1) Effective vertical stress (kN/m2 = kPa)
    kNot = np.zeros(numLayers)        # 2) Coefficient of lateral earth pressure at rest
    refStress = np.zeros(numLayers)   # 3) Reference mean effective confining pressure (kN/m2 = kPa)
    Gref = np.zeros(numLayers)        # 4) Reference low-strain shear modulus (kN/m2 = kPa)

    # ------------------------------------------------------------------------------------------------------------------
    # Compute the reference parameters for the top layer
    # ------------------------------------------------------------------------------------------------------------------

    # 1) Effective vertical stress, vertStress (kN/m2)
    if thick[0] > waterDepth:  # just in the case where waterDepth = 0
        vertStress[0] = gamma[0] * thick[0] * refDepth[0] - gamma_w * (thick[0] * refDepth[0] - waterDepth)
    else:
        vertStress[0] = gamma[0] * thick[0] * refDepth[0]

    # 2) Coefficient of lateral earth pressure at rest, kNot
    kNot[0] = 1.0 - sin(phi[0] * 2.0 * pi / 360.0)

    # 3) Reference mean effective confining pressure, refStress (kN/m2)
    # p' = (sigma'_v + 2*sigma'_h)/3 = (sigma'_v + 2*KNot*sigma'_v)/3 = sigma'_v*(1 + 2*kNot)/3
    refStress[0] = vertStress[0] * (1.0 + 2.0 * kNot[0]) / 3.0

    # 4) Reference low-strain shear modulus, Gref (kN/m2 = kPa)
    # Vs_not (m/s) is the constant required for an exponential function of Vs to have equal travel time to a layer of
    # constant Vs
    Vs_not = Vs[0] / ((1.0 - pressDependCoe[0] / 2.0) * thick[0] ** (pressDependCoe[0] / 2))
    if VsInversionTop == 'No':
        Gref[0] = (gamma[0] / g) * (Vs_not ** 2.0) * (gamma[0] * (1.0 + 2.0 * kNot[0]) / (3.0 * refStress[0])) ** (-pressDependCoe[0])
        # Gref[0] = rho[0] * (Vs_not ** 2.0) * ((thick[0] * refDepth[0]) ** pressDependCoe[0])
    else:
        Gref[0] = (gamma[0] / g) * Vs[0] ** 2.0

    # ------------------------------------------------------------------------------------------------------------------
    # Compute the reference parameters for all other layers
    # ------------------------------------------------------------------------------------------------------------------

    # Compute depth to bottom and midpoint of each layer (m)
    bottomDepth = np.zeros(numLayers)
    midDepth = np.zeros(numLayers)
    bottomDepth[0] = thick[0]
    midDepth[0] = thick[0] / 2.0
    for i in range(1, numLayers):
        bottomDepth[i] = bottomDepth[i-1] + thick[i]
        midDepth[i] = (bottomDepth[i] + bottomDepth[i-1]) / 2.0

    # For each layer compute:
    # 1) Effective vertical stress, vertStress (kN/m2)
    # 2) Coefficient of lateral earth pressure at rest, kNot
    # 3) Reference mean effective confining pressure, refStress (kN/m2)
    for i in range(1, numLayers):
        if bottomDepth[i] > waterDepth and bottomDepth[i-1] > waterDepth:
            vertStress[i] = vertStress[i - 1] + ((gamma[i - 1] - gamma_w) * thick[i - 1] * (1 - refDepth[i - 1])) + ((gamma[i] - gamma_w) * thick[i] * refDepth[i])
        elif bottomDepth[i] > waterDepth:
            vertStress[i] = vertStress[i - 1] + (gamma[i - 1] * thick[i - 1] * (1 - refDepth[i - 1])) + ((gamma[i] - gamma_w) * thick[i] * refDepth[i])
        else:
            vertStress[i] = vertStress[i - 1] + (gamma[i - 1] * thick[i - 1] * (1 - refDepth[i - 1])) + (gamma[i] * thick[i] * refDepth[i])
        kNot[i] = 1.0 - sin(phi[i] * 2.0 * pi / 360.0)
        refStress[i] = vertStress[i] * (1.0 + 2.0 * kNot[i]) / 3.0

    # 4) Reference low-strain shear modulus, Gref (kN/m2)
    # For top layer function of Vs_not and refPress to match pressure independent (PI) travel time
    # For all other layers, directly computed from Vs (of constant profile)
    for i in range(1, numLayers):
        if i == 0:
            Gref[i] = (gamma[i] / g) * (Vs_not ** 2.0) * (gamma[0] * (1.0 + 2.0 * kNot[0]) / (3.0 * refStress[0])) ** (-pressDependCoe[0])
        else:
            Gref[i] = (gamma[i] / g) * Vs[i] ** 2.0

    # ------------------------------------------------------------------------------------------------------------------
    # Compute the pressure dependent shear modulus
    # ------------------------------------------------------------------------------------------------------------------

    # Define the depths to consider for evaluating and plotting
    totalDepth = np.sum(thick)  # Total depth of the deposit (m)
    depthIncr = 0.001  # Increment of depth (m)
    depths = np.arange(0, totalDepth + depthIncr, depthIncr) # Return evenly spaced values within the interval 'depthIcr' (m)
                                                             # It doesn't include the value 'totalDepth + depthIncr'

    # Effective vertical stress at the bottom of each layer, bottomVertStress (kN/m2)
    bottomVertStress = np.zeros(len(bottomDepth))  # Define the size of the array
    # Top layer
    if thick[0] > waterDepth:  # Just in the case where waterDepth = 0
        bottomVertStress[0] = (gamma[0] - gamma_w) * thick[0]
    else:
        bottomVertStress[0] = gamma[0] * thick[0]
    # All other layers
    for i in range(1, numLayers):
        if bottomDepth[i] > waterDepth:
            bottomVertStress[i] = bottomVertStress[i - 1] + (gamma[i] - gamma_w) * thick[i]
        else:
            bottomVertStress[i] = bottomVertStress[i - 1] + gamma[i] * thick[i]

    # Generate two arrays to be used in subsequent computations
    # np.hstack: stack arrays in sequence
    interfaceDepth = np.hstack(([0], bottomDepth)) # Array with depths of interfaces between layers, including 0 (m)
    interfaceVertStress = np.hstack(([0], bottomVertStress)) # Array with vertical effective stresses at interfaces between layers, including 0 (kPa) at surface (kN/m2 = kPa)

    # Generate an array containing the vertical effective stresses evaluated in 'depths' (kN/m2)
    # np.interp: one-dimensional linear interpolation for monotonically increasing sample points
    # np.iterp(x,xp,fp)
    # x: The x-coordinates at which to evaluate the interpolated values
    # xp: The x-coordinates of the data points. Must be increasing if argument period is not specified
    # fp: The y-coordinates of the data points. Same length as xp
    vertStress = np.interp(depths, interfaceDepth, interfaceVertStress) #lin interp between stresses at layer bottoms

    #### compute pressure dependent shear modulus
    # Define the size of the arrays to be generated
    G_independ = np.zeros(len(depths))           # Reference (pressure-independent) low-strain shear modulus (kN/m2)
    G_depend = np.zeros(len(depths))             # Pressure-dependent low-strain shear modulus (kN/m2)
    Vs_depend  = np.zeros(len(depths))           # Pressure-dependent shear wave velocity profile (m/s)
    VsArray = np.zeros(len(depths))              # Pressure-independent shear wave velocity profile (m/s)
    pressDependCoeArray = np.zeros(len(depths))  # Pressure dependent coefficient
    rhoArray = np.zeros(len(depths))             # Density (Mg/m3)
    refStressArray = np.zeros(len(depths))       # Reference mean effective confining pressure, p' (kN/m2)
    meanStress = np.zeros(len(depths))           # Mean effective confining pressure, p'r (kN/m2)

    limit = np.zeros(len(depths))

    # Fill some of the arrays with the information already known
    # Loop over layers
    for i in range(0, numLayers):
        # Loop over 'depths'
        for j in range(0, len(depths)):
            # Check if the depth depths[i] is in layer i
            if depths[j] <= bottomDepth[i] and G_independ[j] == 0.:
                G_independ[j] = Gref[i]
                pressDependCoeArray[j] = pressDependCoe[i]
                refStressArray[j] = refStress[i]
                rhoArray[j] = gamma[i] / g
                VsArray[j] = Vs[i]
                meanStress[j] = vertStress[j] * (1.0 + 2.0 *kNot[i]) / 3.0
                limit[j] = 20

    # Compute the pressure-dependent low-strain shear modulus and shear wave velocity for each depth in 'depths'
    for i in range(0, len(depths)):
        G_depend[i] = G_independ[i] * (meanStress[i] / refStressArray[i]) ** pressDependCoeArray[i]
        Vs_depend[i] = np.sqrt(G_depend[i] / rhoArray[i])

    # Compute the travel time of shear wave through both pressure-independent and -dependent profiles

    # Travel time (s) of the pressure-independent Vs profile
    travelTimeLayerIndep = np.zeros(numLayers)
    for i in range(0, numLayers):
        travelTimeLayerIndep[i] = thick[i]/Vs[i]
    travelTimeIndep = np.sum(travelTimeLayerIndep)

    # Travel time (s) of the pressure-dependent Vs profile
    travelTimeLayerDepend = np.zeros(numLayers)
    # integrate.cumtrapz(y,x): Cumulatively integrate y(x) using the composite trapezoidal rule
    travelTimeDependIncr = integrate.cumtrapz(1 / Vs_depend[1:], depths[1:])
    # integrate.trapz(y,x): Integrate y(x) using the composite trapezoidal rule.
    travelTimeDependTotal = integrate.trapz(1 / Vs_depend[1:], depths[1:])
    # Travel time (s) of the pressure-dependent Vs profile, by layer
    for i in range(0, numLayers):
        if i == 0:
            travelTimeLayerDepend[i] = travelTimeDependIncr[int(bottomDepth[i] / depthIncr - 2.0)] # Need to subtract 1 b/c leaving out first increment from 0 to depthIncr...
                                                                                              # and subtract 1 b/c 1st index is 0 index
        else:
            travelTimeLayerDepend[i] = travelTimeDependIncr[int(bottomDepth[i] / depthIncr - 2.0)] - sum(travelTimeLayerDepend[0:i])

    # ------------------------------------------------------------------------------------------------------------------
    # PLOTTING
    # ------------------------------------------------------------------------------------------------------------------

    plt.figure(figsize=(12, 14))

    plt.plot(VsArray[1:], depths[1:], 'k-')
    plt.plot(Vs_depend[1:], depths[1:], 'b-')
    plt.legend(['Pressure Independent (PI)', 'Pressure Dependent (PD)'], prop={'size':18}, fontsize=30, loc='lower left')
    plt.xlabel('Vs (m/s)', size=18)
    plt.ylabel('Depth (m)', size=18)
    #plt.title('%s' %siteID, y = 1.07, size=26)
    plt.xlim([0, None])
    plt.ylim([totalDepth, 0])

    # Plot the pressure dependent coefficients
    plt.text(Vs[0] - 50.0, midDepth[0], 'd=%.2f' %pressDependCoe[0], ha='right', va='top', fontsize=14)
    for i in range(1, numLayers):
        plt.text(Vs[i] - 50.0, midDepth[i], 'd=%.2f' %pressDependCoe[i], ha='right', va='center', fontsize=14)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    plt.tight_layout()

    plt.savefig(os.path.join(PathOutputs, '%s_PD_Vs.pdf' %siteID), format='pdf')

    print(siteID)
    print('Travel time PI Vs:', travelTimeIndep)
    print('Travel time PD Vs:', travelTimeDependTotal)

    return depths, VsArray, Vs_depend