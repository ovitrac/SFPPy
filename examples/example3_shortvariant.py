#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Variant of example3 including potential release (PR)

Created on Mon Sep  1 17:35:17 2025

@author: olivi
"""
# =============================================================================
# SFPPy Path Bootstrap (allows running from examples/ without installation)
# =============================================================================
import sys
from pathlib import Path
_sfppy_root = Path(__file__).resolve().parent.parent
if str(_sfppy_root) not in sys.path:
    sys.path.insert(0, str(_sfppy_root))
# =============================================================================

from patankar.loadpubchem import migrant
from patankar.layer import gPET,wPET,PP
from patankar.geometry import Packaging3D
from patankar.food import setoff,ambient,hotfilled,realfood,liquid,fat
from patankar.migration import CFSimulationContainer as store

# --- 1) Geometry ----------------------------------------------------------
box = Packaging3D("box_container", length=(19,"cm"),
                  width=(10,"cm"), height=(8,"cm"))
Vint, Acontact = box.get_volume_and_area()

# --- 2) Migrant and layers ------------------------------------------------
m = migrant("limonene")
A1 = wPET(l=(30,"um"), migrant=m, C0=0)
B  = PP(l=(0.5,"mm"), migrant=m, CP0=200) # CP0 or C0 are aliases
A2 = gPET(l=(30,"um"), migrant=m, C0=0)
ABA = A1 + B + A2  # trilayer structure

# --- 3) Contact conditions ------------------------------------------------
class Step1(setoff,ambient): contacttime=(4,"months"); contacttemperature=(20,"C")
class Step2(hotfilled,realfood,liquid,fat): pass
class Step3(ambient,realfood,liquid,fat): contacttime=(6,"months"); contacttemperature=(30,"C")
F1,F2,F3 = Step1(), Step2(), Step3()
box >> F1 >> F2 >> F3  # propagate geometry

# --- 4) Chained simulation (>> operator) ----------------------------------
F1 >> ABA >> F1 >> F2 >> F3
sol_ref = F1.lastsimulation + F2.lastsimulation + F3.lastsimulation

# --- 5) Variants ----------------------------------------------------------
# operator % is used to inject substances
m2 = migrant("toluene")
m2 % F1 >> ABA >> F1 >> F2 >> F3
sol_v1 = F1.lastsimulation + F2.lastsimulation + F3.lastsimulation
m % F1 >> ABA.copy(l=[10e-6,0.5e-3,10e-6],migrant=m) >> F1>>F2>>F3
sol_v2 = F1.lastsimulation + F2.lastsimulation + F3.lastsimulation
m2 % F1 @ ABA.copy(l=[10e-6,0.5e-3,10e-6],migrant=m2) >> F1>>F2>>F3
sol_v3 = F1.lastsimulation + F2.lastsimulation + F3.lastsimulation


# --- 6) Compare all cases -------------------------------------------------
comp = store(name="ABA study")
comp.add(sol_ref, "Limonene with 20 µm thick FB (ref)", "Teal", linewidth=3, linestyle='--')
comp.add(sol_v1, "Toluene with 20 µm thick FB (v1)", "Crimson", linewidth=3)
comp.add(sol_v2, "Limonene with 10 µm thick FB (v2)", "deepskyblue", linewidth=2, linestyle='--')
comp.add(sol_v3, "Toluene-FB with 10 µm thick FB (v3)", "tomato", linewidth=2)
hfig = comp.plotCF()
hfig.print('cmp_pltCF_container_with_FB',destinationfolder="tmp/",overwrite=True)



# %% Example 3 - Extension (1/2) - Potential Release (PR)
# ---------------------------------------------------------
"""
Analysis at the scale of the full packaging (the initial distribution of contaminants is preserved)
    PR = PRE * PRT

PR : ndarray of shape (n_steps,)
    Array of **intrinsic potential release values** for each step in the
    current sequence of transfer (e.g. storage, hot fill, long-term storage).

"""
import numpy as np
# effective potential release of steps 1+2+3, 2+3, 1+2
# Fj @ ABA or Fj >> ABA = thermalize
# ABA >> Fj >> Fj+1 propgates migration from ABA to Fj and then Fj+1
sim123 = m2 % F1 >> ABA >> F1 >> F2 >> F3
PRall123 = np.array([m.lastsimulation.PR.PRtarget_effective for m in [F1,F2,F3]])
F2 @ ABA >> F2 >> F3 # step 1 omitted
PRall23 = np.array([m.lastsimulation.PR.PRtarget_effective for m in [F2,F3]])
F1 @ ABA >> F1 >> F3 # step 2 omitted
PRall13 = np.array([m.lastsimulation.PR.PRtarget_effective for m in [F1,F3]])

# intrinsic potential release for steps, 1=1+2+3\2+3,2=1+2+3\1+3,3=1+2+3\1+2
PR = np.zeros((len(PRall123),),dtype=float)
PRN = PRall123[-1,0] # potential release with all steps
PR[0] = 1.0 - (1.0-PRN)/(1.0-PRall23[1,0])
PR[1] = 1.0 - (1.0-PRN)/(1.0-PRall13[1,0])
PR[2] = 1.0 - (1.0-PRN)/(1.0-PRall123[1,0])
print("PR =", PR)
print("PRT =", PR/F3.lastsimulation.PR.PRE_effective)
print("score = PR/PRN =", PR/PRN)
# compute the error on intrinsic PR
PRNcontrol = (1-np.prod(1-PR))
print(f"PR dependency ~ {(PRNcontrol/PRN - 1).item()*100:0.3f} %")


# %% Example 3 - Extension (2/2) - layer-based Potential Release (PR)
# --------------------------------------------------------------------
"""
Analysis at the scale of each layer (j=1,2,3).

PR_layers : ndarray of shape (n_steps, n_layers)
    Array of **intrinsic potential release values per step and per layer**
    in a multilayer system, reconstructed by comparison of full and
    step-omitted sequences.

"""

R123 = F1.potentialRelease(ABA) >> F2 >> F3
R23 = F2.potentialRelease(ABA) >> F3
R13 = F1.potentialRelease(ABA) >> F3

# intrinsic layer-based potential release
PR_layers = np.zeros((len(R123),len(ABA)),dtype=float) #nsteps x nlayers
PRN_layers = R123.PR[2,:]
PR_layers[0,:] = 1.0 - (1.0-PRN_layers)/(1.0-R23.PR[1,:])
PR_layers[1,:] = 1.0 - (1.0-PRN_layers)/(1.0-R13.PR[1,:])
PR_layers[2,:] = 1.0 - (1.0-PRN_layers)/(1.0-R123.PR[1,:])
PRN_layers_control = 1 - np.prod(1-PR_layers,axis=0)
