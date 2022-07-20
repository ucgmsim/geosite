"""
GeoSite - Version 1.2
Authors: Felipe Kuncar (felipe.kuncar@pg.canterbury.ac.nz | fkw21)
         Christopher de la Torre (christopher.delatorre@canterbury.ac.nz | cde62)
Created on 16-09-2021
Last modification on 04-02-2022
Python Version: 3.8
"""

# =====================================================================================================================
# IMPORT LIBRARIES, FUNCTIONS, CLASSES
# =====================================================================================================================

import numpy as np
import pandas as pd
import os.path
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import time

from GeoSite import CPTu
from GeoSite import readVsData
from GeoSite import Vs
from GeoSite import trialModelsFromCPTu
from GeoSite import PDMY02
from GeoSite import create_inputOpenSees
from GeoSite import create_modelOpenSees
from GeoSite import PressDepVs

start = time.time()

# =====================================================================================================================
# SITE INFORMATION
# =====================================================================================================================

# Site ID
SiteID = 'LRSS'

# Water table (positive value in (m))
# (always at the interface of two model layers)
GWL = 3.50

# =====================================================================================================================
# PATHS
# =====================================================================================================================

# Define the path to the root directory
rootDir = 'C:/Users/Felipe Kuncar/Desktop/PhD/7. RESEARCH/6. Objectives/Obj1/3. Site Characterization/GeoSite/version 1/version 1.2'

# Define the path to the directory that contains the site characterization data to read
inputPath = os.path.join(rootDir, 'Inputs')
# Define the path to the directory that will save the results
outputPath = os.path.join(rootDir, 'Outputs')
if not os.path.exists(os.path.join(rootDir, 'Outputs')):
    os.makedirs(os.path.join(rootDir, 'Outputs'))

# =====================================================================================================================
# INITIALIZE LISTS THAT WILL SAVE DATA
# =====================================================================================================================

CPTu_Path = []
CPTu_Data = []
SCPT_Path = []
SCPT_Data = []
SW_Path = []
SW_Data = []

# =====================================================================================================================
# FILES TO BE READ
# =====================================================================================================================

# CPTu
# Format (Excel): depth (m) | qc (kPa) | fs (kPa) | u2 (kPa)
CPTu_Path += [os.path.join(inputPath, 'LRSS_CPTu1.xlsx')]

# SCPT
# Format (Excel): depth (m) | Vs (m/s)
SCPT_Path += [os.path.join(inputPath, 'LRSS_SCPT1.xlsx')]

# SW testing
# Format (csv file): data_point_number,depth (m),Vs(m/s)
SW_Path += [os.path.join(inputPath, 'LRSS_LR15_vs_profiles.csv')]

# =====================================================================================================================
# READ DATA AND DEFINE OBJECTS
# =====================================================================================================================

### Define CPTu objects

# CPTu 1
# Read data
CPTu_read = pd.read_excel(CPTu_Path[0])
# Assign class CPTu
CPTu_depth = CPTu_read.iloc[:, 0]
qc = CPTu_read.iloc[:, 1]
fs = CPTu_read.iloc[:, 2]
u2 = CPTu_read.iloc[:, 3]
CPTu_Data.append(CPTu(CPTu_depth, qc, fs, u2, GWL))

### Define Vs objects

# From SCPT

# SCPT 1
# Read data
Vs_read = pd.read_excel(SCPT_Path[0])
# Assign class Vs
Vs_depth = Vs_read.iloc[:, 0]
Vs_values = Vs_read.iloc[:, 1]
SCPT_Data.append(Vs(Vs_depth, Vs_values, GWL))

# From SW testing

# SW 1
# Read data
Vs_depth, Vs_values = readVsData(SW_Path[0])
# Assign class Vs
SW_Data.append(Vs(Vs_depth, Vs_values, GWL))
# Layering ratio used
LR = 1.5

# =====================================================================================================================
# ADJUST SITE-RESPONSE MODEL PARAMETERS
# =====================================================================================================================

# Model layers
layersThickness_OpenSees = np.array([1.40, 2.10, 1.50, 1.30, 1.90, 1.275, 1.275, 5.02, 3.44, 1.00])

# Model parameters
# Manually adjusted by the user
# Constitutive model flag: 1 for PDMY02 (sands) | 2 for PIMY (clays)
# TODO: add model flag 0 for lineal-elastic model
constModelFlag = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
Ic_OpenSees = np.array([1.57, 1.20, 1.40, 1.22, 1.84, 2.39, 2.58, None, None, None])
Vs_OpenSees = np.array([110.0, 195.0, 195.0, 195.0, 175.0, 165.0, 165.0, 197.0, 273.0, 500.0])
gamma_OpenSees = np.array([18.87, 19.92, 18.98, 20.19, 19.12, 18.45, 18.08, 18.00, 18.80, 20.00])
Dr_OpenSees = np.array([71.0, 99.5, 89.6, 99.4, 59.2, 29.7, 19.85, 40.0, 65.0, 100.0])
phi_OpenSees = np.array([38.0, 40.0, 37.0, 40.0, 40.0, 35.0, 35.0, 35.0, 37.0, 40.0])
cohesion_OpenSees = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# Depth-dependent parameters
pressDependCoe = np.array([0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.50, 0.50, 0.50, 0.00])
refDepth = np.array([0.50, 1.00, 0.50, 0.50, 0.50, 0.50, 0.00, 0.50, 0.50, 0.50])

# Plot the pressure dependent Vs profile to manually adjust pressDependCoe and refDepth
VsInversionTop = 'No'  # 'Yes' if the Vs inversion is immediately below top layer, 'No' otherwise or if no inversions present

# =====================================================================================================================
# DEFINE MODEL OBJECTS
# =====================================================================================================================

# Modify the arrays to plot
depth_model = np.zeros(1 + len(layersThickness_OpenSees))
for i in range(1, len(depth_model)):
    depth_model[i] = depth_model[i-1] + layersThickness_OpenSees[i-1]
depth_model = np.repeat(depth_model, 2)
depth_model = np.delete(depth_model, [0, len(depth_model) - 1])
Ic_model = np.repeat(Ic_OpenSees, 2)
gamma_model = np.repeat(gamma_OpenSees, 2)
Vs_model = np.repeat(Vs_OpenSees, 2)
phi_model = np.repeat(phi_OpenSees, 2)
Dr_model = np.repeat(Dr_OpenSees, 2)

### Define the model Vs object
# Assign class Vs
Vs_m = Vs(depth_model, Vs_model, GWL)

### Define PDMY02 object
# Assign class PDMY02
PDMY02_Model = PDMY02(depth_model, Dr_model)

# =====================================================================================================================
# CREATE OUTPUT FILES
# =====================================================================================================================

### CPTu
Ic, Qtn = CPTu_Data[0].Ic(version=1)
I_SBT, none = CPTu_Data[0].Ic(version=2)

basic_results = np.array([CPTu_Data[0].depth, (CPTu_Data[0].qc / 1000), CPTu_Data[0].fs, CPTu_Data[0].u2, CPTu_Data[0].u0, (CPTu_Data[0].qt / 1000), CPTu_Data[0].Rf, CPTu_Data[0].gamma, CPTu_Data[0].sigma_v0, I_SBT, Qtn, CPTu_Data[0].Fr, Ic])
basic_results = pd.DataFrame(np.transpose(basic_results), columns=['Depth (m)', 'qc (Mpa)', 'fs (kPa)', 'u2 (kPa)', 'u0 (kPa)', 'qt (MPa)', 'Rf (%)', 'Unit weight (kN/m3)', 'sigma_v0 [kPa]', 'I SBT', 'Qtn', 'Fr (%)', 'Ic'])
basic_results.to_excel(os.path.join(outputPath, '%s_CPTu_basic_results.xlsx' % SiteID))

estimated_parameters = np.array([CPTu_Data[0].depth, CPTu_Data[0].Vs(Ic, correlation=1), CPTu_Data[0].gamma, CPTu_Data[0].relative_density(Ic, Qtn, correlation=1), CPTu_Data[0].friction_angle(Ic, Qtn, correlation=1)])
estimated_parameters = pd.DataFrame(np.transpose(estimated_parameters), columns=['Depth (m)', 'Vs (m/s)', 'Unit weight (kN/m3)', 'Dr (%)', 'Friction angle (°)'])
estimated_parameters.to_excel(os.path.join(outputPath, '%s_CPTu_estimated_parameters.xlsx' % SiteID))

# Suggested or trial models obtained from CPTu-based correlations
layersThickness_trial, Ic_trial, Vs_trial1, Vs_trial2, gamma_trial, Dr_trial1, Dr_trial2, phi_trial1, phi_trial2, phi_trial3 = trialModelsFromCPTu(layersThickness_OpenSees, CPTu_Data[0])
trial_model_CPTu = np.array([layersThickness_trial, Ic_trial, Vs_trial1, Vs_trial2, gamma_trial,  Dr_trial1, Dr_trial2, phi_trial1, phi_trial2, phi_trial3])
trial_model_CPTu = pd.DataFrame(np.transpose(trial_model_CPTu), columns=['Thickness (m)', 'Ic', 'Vs (m/s) | McGann et al. (2015)', 'Vs (m/s) | CPT Guide (2015)', 'Unit weight (kN/m3) | Robertson & Cabal (2010)', 'Dr (%) | Idriss & Boulanger (2008)', 'Dr (%) | Kulhawy & Mayne (1990)', 'Friction angle (°) | Robertson (2010)', 'Friction angle (°) | Kulhawy & Mayne (1990)', 'Friction angle (°) | Robertson & Campanella (1983)'])
trial_model_CPTu.to_excel(os.path.join(outputPath, '%s_CPTu_trial_models.xlsx' % SiteID))

### Vs
N1_60 = Vs_m.N1_60(correlation=1)

### Write tcl file for OpenSees
# Only input, for running analyses in laptop
create_inputOpenSees(SiteID, outputPath, layersThickness_OpenSees, GWL, gamma_OpenSees, Vs_OpenSees, phi_OpenSees, Dr_OpenSees, refDepth, pressDependCoe)
# Complete model, for running analyses in workflow implemented in Mahuika
create_modelOpenSees(SiteID, outputPath, layersThickness_OpenSees, GWL, gamma_OpenSees, LR, Vs_OpenSees, constModelFlag, phi_OpenSees, Dr_OpenSees, refDepth, pressDependCoe)

# =====================================================================================================================
# CREATE FIGURES
# =====================================================================================================================

# x limit
x_lim_Vs = max(Vs_model) + 30
# y limit
y_lim_cpt = max(CPTu_Data[0].depth)
y_lim_model = max(depth_model) + 5

#----------------------------------------------------------------------------------------------------------------------
# PRESSURE-DEPENDENT VS PROFILE
# ---------------------------------------------------------------------------------------------------------------------

PressDepVs(outputPath, layersThickness_OpenSees, Vs_OpenSees, gamma_OpenSees, phi_OpenSees, SiteID, pressDependCoe, refDepth,
           VsInversionTop, GWL)

# ---------------------------------------------------------------------------------------------------------------------
# BASIC CPTu FIGURE
# ---------------------------------------------------------------------------------------------------------------------

# Create a figure
fig1 = plt.figure(figsize=(15, 8))

# Title
fig1.suptitle("CPTu Information - CPTu1 - Site %s" % SiteID, fontsize=16)

# Create the first gridspec for ploting the model biases
gs1 = fig1.add_gridspec(nrows=1, ncols=5, left=0.05, bottom=0.08, right=0.97, top=0.90, wspace=0.35)

# Cone resistance (MPa)
ax1 = fig1.add_subplot(gs1[0, 0])
ax1.plot(CPTu_Data[0].qc / 1000, CPTu_Data[0].depth, color='0.6', linewidth=1.0, label='$q_c$')
ax1.plot(CPTu_Data[0].qt / 1000, CPTu_Data[0].depth, 'k-', linewidth=1.0, label='$q_t$')
ax1.legend(loc=4)
ax1.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax1.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax1.set_title('Cone resistance', size=14)
ax1.set_xlabel('Tip resistance (MPa)', size=12)
ax1.set_ylabel('Depth (m)', size=12)
ax1.set_xlim([min(CPTu_Data[0].qt) / 1000, max(CPTu_Data[0].qt) / 1000])
ax1.set_xlim([0, None])
ax1.set_ylim([0, y_lim_cpt])
plt.gca().invert_yaxis()

# Sleeve friction (kPa)
ax2 = fig1.add_subplot(gs1[0, 1])
ax2.plot(CPTu_Data[0].fs, CPTu_Data[0].depth, 'k-', linewidth=1.0)
ax2.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax2.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax2.set_title('Sleeve friction', size=14)
ax2.set_xlabel('$f_s$ (kPa)', size=12)
ax2.set_ylabel('Depth (m)', size=12)
ax2.set_xlim([min(CPTu_Data[0].fs), max(CPTu_Data[0].fs)])
ax2.set_xlim([0, None])
ax2.set_ylim([0, y_lim_cpt])
plt.gca().invert_yaxis()

# Friction ratio (%)
ax3 = fig1.add_subplot(gs1[0, 2])
ax3.plot(CPTu_Data[0].Rf, CPTu_Data[0].depth, 'k-', linewidth=1.0)
ax3.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax3.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax3.set_title('Friction ratio', size=14)
ax3.set_xlabel('$R_s$ (%)', size=12)
ax3.set_ylabel('Depth (m)', size=12)
ax3.set_xlim([min(CPTu_Data[0].Rf), max(CPTu_Data[0].Rf)])
ax3.set_xlim([0, 10])
ax3.set_ylim([0, y_lim_cpt])
plt.gca().invert_yaxis()

# Pore pressure (kPa)
ax4 = fig1.add_subplot(gs1[0, 3])
ax4.plot([np.minimum(np.min(CPTu_Data[0].u2), np.min(CPTu_Data[0].u0)), np.maximum(np.max(CPTu_Data[0].u2), np.max(CPTu_Data[0].u0))], [CPTu_Data[0].GWL, CPTu_Data[0].GWL], color='cyan', linestyle=(0, (5, 5)), linewidth=1.5)
ax4.plot(CPTu_Data[0].u0, CPTu_Data[0].depth, 'b-', linewidth=1.0, label='$u_0$')
ax4.plot(CPTu_Data[0].u2, CPTu_Data[0].depth, 'k-', linewidth=1.0, label='$u_2$')
ax4.legend(loc=4)
ax4.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax4.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax4.set_title('Pore pressure (kPa)', size=14)
ax4.set_xlabel('Pressure (kPa)', size=12)
ax4.set_ylabel('Depth (m)', size=12)
ax4.set_xlim([min(CPTu_Data[0].u2), max(CPTu_Data[0].u2)])
ax4.set_ylim([0, y_lim_cpt])
plt.gca().invert_yaxis()

# Soil behaviour type index
ax5 = fig1.add_subplot(gs1[0, 4])
ax5.plot(Ic, CPTu_Data[0].depth, color='yellow', linewidth=1.5)
ax5.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax5.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax5.set_title('Soil Behaviour Type Index', size=14)
ax5.set_xlabel('$I_{c}$', size=12)
ax5.set_ylabel('Depth (m)', size=12)
ax5.set_xlim([1, 4])
ax5.set_ylim([0, y_lim_cpt])
ax5.add_patch(patches.Rectangle((0.00, 0.00), 1.31, y_lim_cpt, facecolor='whitesmoke'))
ax5.add_patch(patches.Rectangle((1.31, 0.00), 0.74, y_lim_cpt, facecolor='gainsboro'))
ax5.add_patch(patches.Rectangle((2.05, 0.00), 0.55, y_lim_cpt, facecolor='silver'))
ax5.add_patch(patches.Rectangle((2.60, 0.00), 0.35, y_lim_cpt, facecolor='darkgray'))
ax5.add_patch(patches.Rectangle((2.95, 0.00), 0.65, y_lim_cpt, facecolor='gray'))
ax5.add_patch(patches.Rectangle((3.60, 0.00), 0.40, y_lim_cpt, facecolor='dimgray'))
ax5.text(1.00+(1.31-1.00)/2, y_lim_cpt/2, "gravelly sand - dense sand", ha='center', va='center', rotation=90, size=10, color="black")
ax5.text(1.31+0.74/2, y_lim_cpt/2, " sands: clean sand - silty sand", ha='center', va='center', rotation=90, size=10, color="black")
ax5.text(2.05+0.55/2, y_lim_cpt/2, "sand mixtures: silty sand - sandy silt", ha='center', va='center', rotation=90, size=10, color="black")
ax5.text(2.60+0.35/2, y_lim_cpt/2, "silt mixtures: clayey silt - silty clay", ha='center', va='center', rotation=90, size=10, color="black")
ax5.text(2.95+0.65/2, y_lim_cpt/2, "clays: silty clay - clay", ha='center', va='center', rotation=90, size=10, color="black")
ax5.text(3.60+(4.00-3.60)/2, y_lim_cpt/2, "organic soils: peats", ha='center', va='center', rotation=90, size=10, color="black")
plt.gca().invert_yaxis()

# Save figure
plt.savefig(os.path.join(outputPath, '%s_CPTu_basic.pdf' % SiteID), dpi=300)

# ---------------------------------------------------------------------------------------------------------------------
# ESTIMATED PROPERTIES
# ---------------------------------------------------------------------------------------------------------------------

# Create a figure
fig2 = plt.figure(figsize=(15, 8))

# Title
fig2.suptitle("Parameters Estimate - Site %s" % (SiteID), fontsize=16)

# Create the first gridspec for ploting the model biases
gs2 = fig2.add_gridspec(nrows=1, ncols=5, left=0.05, bottom=0.08, right=0.97, top=0.90, wspace=0.35)

# Soil behaviour type index
ax1 = fig2.add_subplot(gs2[0, 0])
ax1.plot(Ic, CPTu_Data[0].depth, color='yellow', linewidth=1.2)
ax1.plot(Ic_model, depth_model, 'r', linewidth=1.2, linestyle=(0, (5, 3)))
for i in range (1, len(depth_model)):
    ax1.plot([1, 4], [depth_model[i], depth_model[i]], color='lightcoral', linewidth=0.2)
ax1.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax1.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax1.set_title('Soil Behaviour Type Index', size=14)
ax1.set_xlabel('$I_{c}$', size=12)
ax1.set_ylabel('Depth (m)', size=12)
ax1.set_xlim([1, 4])
ax1.set_ylim([0, y_lim_model])
ax1.add_patch(patches.Rectangle((0.00, 0.00), 1.31, y_lim_model, facecolor='whitesmoke'))
ax1.add_patch(patches.Rectangle((1.31, 0.00), 0.74, y_lim_model, facecolor='gainsboro'))
ax1.add_patch(patches.Rectangle((2.05, 0.00), 0.55, y_lim_model, facecolor='silver'))
ax1.add_patch(patches.Rectangle((2.60, 0.00), 0.35, y_lim_model, facecolor='darkgray'))
ax1.add_patch(patches.Rectangle((2.95, 0.00), 0.65, y_lim_model, facecolor='gray'))
ax1.add_patch(patches.Rectangle((3.60, 0.00), 0.40, y_lim_model, facecolor='dimgray'))
ax1.text(1.00+(1.31-1.00)/2, y_lim_model/2, "gravelly sand - dense sand", ha='center', va='center', rotation=90, size=10, color="black")
ax1.text(1.31+0.74/2, y_lim_model/2, " sands: clean sand - silty sand", ha='center', va='center', rotation=90, size=10, color="black")
ax1.text(2.05+0.55/2, y_lim_model/2, "sand mixtures: silty sand - sandy silt", ha='center', va='center', rotation=90, size=10, color="black")
ax1.text(2.60+0.35/2, y_lim_model/2, "silt mixtures: clayey silt - silty clay", ha='center', va='center', rotation=90, size=10, color="black")
ax1.text(2.95+0.65/2, y_lim_model/2, "clays: silty clay - clay", ha='center', va='center', rotation=90, size=10, color="black")
ax1.text(3.60+(4.00-3.60)/2, y_lim_model/2, "organic soils: peats", ha='center', va='center', rotation=90, size=10, color="black")
plt.gca().invert_yaxis()

# Shear-wave velocity (m/s)
ax2 = fig2.add_subplot(gs2[0, 1])
# Direct measurements:
ax2.plot(SCPT_Data[0].values, SCPT_Data[0].depth, 'k-', linewidth=1.2, label='SCPT test')
ax2.plot(SW_Data[0].values, SW_Data[0].depth, color='grey', linestyle=(0, (1, 1)), linewidth=1.1, label='SW testing')
# CPTu-based correlations:
ax2.plot(CPTu_Data[0].Vs(Ic, correlation=1), CPTu_Data[0].depth, 'b', linewidth=1.0, label='CPTu - McGann et al. (2015)')
ax2.plot(CPTu_Data[0].Vs(Ic, correlation=2), CPTu_Data[0].depth, 'c', linewidth=1.0, label='CPTu - CPT Guide (2015)')
# Site-response model:
ax2.plot(Vs_model, depth_model, 'r', linewidth=1.4, linestyle=(0, (5, 3)), label='Site-response model')
ax2.legend(loc=3, fontsize=7)
ax2.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax2.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax2.set_title('Shear-wave velocity', size=14)
ax2.set_xlabel('$V_s$ (m/s)', size=12)
ax2.set_ylabel('Depth (m)', size=12)
ax2.set_xlim([0, x_lim_Vs])
ax2.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# Soil total unit weight (kPa)
ax3 = fig2.add_subplot(gs2[0, 2])
# CPTu-based correlations:
ax3.plot(CPTu_Data[0].gamma, CPTu_Data[0].depth, 'b-', linewidth=1.0, label='CPTu - Robertson & Cabal (2010)')
# Vs-based correlations:
ax3.plot(Vs_m.gamma_corr(correlation=1), Vs_m.depth, 'green', linewidth=1.1, label='Vs - NCHRP (2019)')
ax3.plot(Vs_m.gamma_corr(correlation=2), Vs_m.depth, color='limegreen', linewidth=1.1, label='Vs - Boore (2015)')
# Site-response model:
ax3.plot(gamma_model, depth_model, 'r', linewidth=1.4, linestyle=(0, (5, 3)), label='Site-response model')
ax3.legend(loc=3, fontsize=7)
ax3.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax3.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax3.set_title('Soil total unit weight', size=14)
ax3.set_xlabel('$\gamma$ (kN/m3)', size=12)
ax3.set_ylabel('Depth (m)', size=12)
ax3.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# Relative density (%)
ax4 = fig2.add_subplot(gs2[0, 3])
# CPTu-based correlations:
ax4.plot(CPTu_Data[0].relative_density(Ic, Qtn, correlation=1), CPTu_Data[0].depth, 'b-', linewidth=1.0, label='CPTu - Idriss & Boulanger (2008)')
ax4.plot(CPTu_Data[0].relative_density(Ic, Qtn, correlation=2), CPTu_Data[0].depth, 'c--', linewidth=1.0, label='CPTu - Kulhawy & Mayne (1990)')
# Vs-based correlations:
ax4.plot(Vs_m.relative_density(N1_60, correlation=1), Vs_m.depth, 'g-', linewidth=1.1, label='Vs - Idriss & Boulanger (2008)')
# Site-response model:
ax4.plot(Dr_model, depth_model, 'r--', linewidth=1.4, label='Site-response model', linestyle=(0, (5, 3)))
ax4.legend(loc=3, fontsize=7)
ax4.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax4.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax4.set_title('Relative density', size=14)
ax4.set_xlabel('$D_r$ (%)', size=12)
ax4.set_ylabel('Depth (m)', size=12)
ax4.set_xlim([0, 101])
ax4.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# Friction angle (°)
ax5 = fig2.add_subplot(gs2[0, 4])
# CPTu-based correlations:
ax5.plot(CPTu_Data[0].friction_angle(Ic, Qtn, correlation=1), CPTu_Data[0].depth, 'b', linewidth=1.0, label='CPTu - Robertson (2010)')
ax5.plot(CPTu_Data[0].friction_angle(Ic, Qtn, correlation=2), CPTu_Data[0].depth, 'darkviolet', linewidth=1.0, label='CPTu - Kulhawy & Mayne (1990)')
ax5.plot(CPTu_Data[0].friction_angle(Ic, Qtn, correlation=3), CPTu_Data[0].depth, 'c', linewidth=1.0, label='CPTu - Robertson & Campanella (1983)')
# Vs-based correlations:
ax5.plot(Vs_m.friction_angle(N1_60, correlation=1), Vs_m.depth, 'g-', linewidth=1.1, label='Vs - NCHRP (2019)')
# Site-response model:
ax5.plot(phi_model, depth_model, 'r', linewidth=1.4, label='Site-response model', linestyle=(0, (5, 3)))
ax5.legend(loc=3, fontsize=7)
ax5.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax5.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax5.set_title('Friction angle', size=14)
ax5.set_xlabel('$\phi$ (°)', size=12)
ax5.set_ylabel('Depth (m)', size=12)
ax5.set_xlim([20, 50])
ax5.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# Save figure
plt.savefig(os.path.join(outputPath, '%s_parameters_estimate.pdf' % SiteID), dpi=300)


# ---------------------------------------------------------------------------------------------------------------------
# OTHERS PDMY02 MODEL PARAMETERS
# ---------------------------------------------------------------------------------------------------------------------

# create a figure
fig3 = plt.figure(figsize=(15, 8))

# title
fig3.suptitle("PDMY02 Parameters Estimate - Site %s" % (SiteID), fontsize=16)

# create the first gridspec for ploting the model biases
gs3 = fig3.add_gridspec(nrows=1, ncols=8, left=0.05, bottom=0.08, right=0.97, top=0.90, wspace=0.35)

# PTAng
ax1 = fig3.add_subplot(gs3[0, 0])
# site-response model:
ax1.plot(PDMY02_Model.PTAng(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax1.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax1.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax1.set_title('$PTAng (°)', size=14)
ax1.set_ylabel('Depth (m)', size=12)
ax1.set_xlim([20, 35])
ax1.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# e
ax2 = fig3.add_subplot(gs3[0, 1])
# site-response model:
ax2.plot(PDMY02_Model.e(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax2.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax2.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax2.set_title('$e', size=14)
ax2.set_xlim([0.50, 0.90])
ax2.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# contrac1
ax3 = fig3.add_subplot(gs3[0, 2])
# site-response model:
ax3.plot(PDMY02_Model.contrac1(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax3.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax3.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax3.set_title('$contrac1', size=14)
ax3.set_xlim([0.00, 0.10])
ax3.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# contrac2
ax4 = fig3.add_subplot(gs3[0, 3])
# site-response model:
ax4.plot(PDMY02_Model.contrac2(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax4.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax4.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax4.set_title('$contrac2', size=14)
ax4.set_xlim([0.00, 6.00])
ax4.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# contrac3
ax5 = fig3.add_subplot(gs3[0, 4])
# site-response model:
ax5.plot(PDMY02_Model.contrac3(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax5.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax5.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax5.set_title('$contrac3', size=14)
ax5.set_xlim([0.00, 0.40])
ax5.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# dilat1
ax6 = fig3.add_subplot(gs3[0, 5])
# site-response model:
ax6.plot(PDMY02_Model.dilat1(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax6.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax6.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax6.set_title('$dilat1', size=14)
ax6.set_xlim([0.00, 0.40])
ax6.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# dilat2
ax7 = fig3.add_subplot(gs3[0, 6])
# site-response model:
ax7.plot(PDMY02_Model.dilat2(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax7.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax7.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax7.set_title('$dilat2', size=14)
ax7.set_xlim([0.00, 4.00])
ax7.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# dilat3
ax8 = fig3.add_subplot(gs3[0, 7])
# site-response model:
ax8.plot(PDMY02_Model.dilat3(), PDMY02_Model.depth, 'r', linewidth=1.4, linestyle=(0, (5, 3)))
ax8.grid(color='0.6', linestyle=(0, (5, 10)), which='major', linewidth=0.3)
ax8.grid(color='0.6', linestyle=(0, (5, 10)), which='minor', linewidth=0.3)
ax8.set_title('$dilat3', size=14)
ax8.set_xlim([-0.01, 0.50])
ax8.set_ylim([0, y_lim_model])
plt.gca().invert_yaxis()

# save figure
plt.savefig(os.path.join(outputPath, '%s_PDMY02_parameters_estimate.pdf' % SiteID), dpi=300)

end = time.time()
print("Done in {:10.1f} secs".format(end-start))
