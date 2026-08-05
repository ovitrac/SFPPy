"""
===============================================================================
SFPPy Example: cosmetic exposure and the tiered decision
===============================================================================

# Example 6: A Recycled PET Bottle for a Body Lotion
---------------------------------------------------

## Overview
This script assesses substances migrating from a **recycled PET bottle** into a
**body lotion**, using the `cosmetic` module. It is the first example that does
not stop at `CF`: cosmetics have no regulatory convention folding exposure into
the limit, so the chain continues to a **dose** and a **verdict**.

The example is built around one idea:

    We do not need to simulate accurately. We need to simulate robustly.
    Refinement is required only when a problem is found -- and then the exit
    is a branch, not an escalation: refine one rung, OR measure.

### Simulation Steps:
1. **The unknown substance** -- T0 bound, no identification needed.
2. **An identified substance that passes** -- toluene, stop there.
3. **An identified substance that fails** -- bisphenol A, and the branch.
4. **Refining one rung** -- M0 total transfer, then M3 kinetics.
5. **Turning the assessment around** -- what the material may contain.
6. **The gates** -- refusals rather than plausible wrong numbers.

### Expected Outcomes:
- A verdict with its **binding threshold, tier and provenance**, never a bare number.
- A demonstration that **M0 bounds M3**, so the cheap answer is the safe one.
- A **design answer** (maximum concentration, required decontamination), not just
  a compliance verdict.
- Gates firing where a silent approximation would have been wrong.

## What You Learn
- Assess a product with `substance % product` and `product.assess()`
- Read `cosmeticassessment`: dose, threshold, tier, margin, provenance
- Use the standard pipeline `substance % product << geometry >> packaging >> product`
- Climb the tiers only on failure, and see `refine()` refuse when there is no trigger
- Invert the question with `maxconcentration()` and `required_decontamination()`

## Provenance
The exposure classes, the simulant doctrine, the three product categories and the
Tier 1 defaults come from the CosPaTox consortium's published work. Any use must
cite the April 2024 deliverables (Guideline, Substance List, Scientific Dossier,
WCC Calculator) at https://cospatox.com/publication/ -- no DOI and no mandated
citation format exist; cite by name, consortium, month/year and URL.

---

@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2026-08-04
"""

# %% Path bootstrap -- run directly, no installation required
import os
import sys

_ex = os.path.dirname(os.path.abspath(__file__))
_root = os.path.dirname(_ex)
if _root not in sys.path:
    sys.path.insert(0, _root)

import numpy as np
if not hasattr(np, "trapz"):          # NumPy 2.0 removed np.trapz
    np.trapz = np.trapezoid

from patankar.loadpubchem import migrantToxtree
from patankar.geometry import Packaging3D
from patankar.layer import wPET
from patankar.cosmetic import (bodylotion, cosmeticgate, recyclate_C0,
                               required_decontamination, washinggel)


def rule(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def scalar(x):
    return float(np.asarray(x).ravel()[-1])


# %% 1 | The unknown substance -------------------------------------------------
# The hardest case in a recyclate is the substance nobody identified. It is also
# the cheapest to answer: at tier T0 the genotoxic TTC applies to *any* substance,
# so a single bound covers the whole unidentified population without naming one
# molecule. Nothing here needs PubChem, a geometry or a solver.

rule("1 | The unidentified substance -- one bound, no molecule named")

lotion = bodylotion()
print("product   : %s (%s, %s)" % (lotion.name, lotion._usage, lotion._target))
print("simulant  : inherited from the type -- EtOH 95 doctrine (hydrophobic)")
print("volume    : %.4g mL, body weight %.3g kg"
      % (scalar(lotion.volume) * 1e6, scalar(lotion.bodyweight)))
# NOTE -- SFPPy stores everything in SI, so `frequency` is per SECOND, not per
# day. Multiply by 86400 to read it back in the unit it was declared in.
print("exposure  : %.3g g applied, retention %.3g, area %.4g cm2, %.3g /day"
      % (scalar(lotion.amountapplied) * 1e3, scalar(lotion.retention),
         scalar(lotion.exposedarea) * 1e4, scalar(lotion.frequency) * 86400.0))

print("""
The commonest substance in a recyclate is the one nobody named: detected below
the limit, or not detected at all, with no structure, no CAS, no Cramer class.
Tier T0 exists for exactly that population -- the genotoxic TTC applies to ANY
molecule, so one bound covers all of them at once.

Declaring it is explicit. Leaving `substance` as None stays a refusal (G2), so
a forgotten injection can never become a silent T0 answer.""")

anon = bodylotion().declare_unidentified()
thr = anon.binding_threshold()
print("\n   threshold : %.4g %s  [%s]  %s / %s"
      % (thr.value, thr.units, thr.identity, thr.tier, thr.origin))
print("   Cmax      : %.4g mg/kg   -- what the PRODUCT may contain" % scalar(anon.maxconcentration()))
print("   DLmax     : %.4g mg/kg   -- what the LABORATORY must reach"
      % scalar(anon.maxdetectionlimit(censoring="DL")))

# A non-detect is left-censored data, not a number, and the convention depends
# on what is known (gate G4): DL/2 for an identified substance, a central
# estimate; DL for the T0 bound, because halving a BOUND is not the same
# operation as taking the midpoint of a censored distribution.
anon.CF0 = 0.01
print("\n   at CF = 0.01 mg/kg counted at DL:")
print("   ", str(anon.assess(censoring="DL")).replace("\n", "\n    "))

# The same bound, on a rinse-off product for a child. Retention 0.01 against 1
# for a leave-on lotion moves the tolerable concentration by more than an order
# of magnitude -- which is the whole reason exposure is modelled per category
# rather than folded into a single limit as the food convention does.
kid = washinggel().declare_unidentified()
print("\n   same bound, washing gel (rinse-off, child):")
print("     Cmax = %.4g mg/kg   (vs %.4g for the leave-on lotion)"
      % (scalar(kid.maxconcentration()), scalar(anon.maxconcentration())))
print("""
     Same toxicology, same threshold, 27x apart. Retention and body weight did
     that, not chemistry. A single food-style limit cannot express it.""")


# %% 2 | An identified substance that passes -----------------------------------
# Identifying a substance moves along the T axis: T0 (genotoxic TTC) -> T1
# (Cramer class) -> T2 (published TDI/SML). Higher tier = better INFORMED --
# which is not the same as less conservative, as step 3 shows. Here
# identification alone settles the case.

rule("2 | Toluene at DL/2 -- identified, and it passes")

toluene = migrantToxtree("toluene")
p = toluene % bodylotion()
p.CF0 = 0.5                                    # mg/kg, non-detect at DL = 1
a_tol = p.assess()
print(a_tol)
print("\nProvenance is part of the answer: %s" % a_tol.provenance)

# G8 -- stop on pass. Refinement without a trigger is a defect, not diligence.
print("\nAsking to refine a passing case:")
try:
    p.refine(a_tol)
except cosmeticgate as e:
    print("   REFUSED -- %s" % str(e).splitlines()[0])


# %% 3 | An identified substance that fails ------------------------------------
# Bisphenol A carries a T2 threshold three orders of magnitude below the Cramer
# TTC. At a detection limit of 50 mg/kg the same arithmetic that cleared toluene
# now fails -- which is the point: the tier that matters is the one that binds.

rule("3 | Bisphenol A at DL/2 = 25 mg/kg -- and it fails")

bpa = migrantToxtree("bisphenol A")
q = bpa % bodylotion()
q.CF0 = 25.0                                   # mg/kg, non-detect at DL = 50
a_bpa = q.assess()
print(a_bpa)

print("""
Note the gate on that result. The shipped EU 10/2011 snapshot still carries an
SML of 0.05 mg/kg for bisphenol A -- but Commission Regulation (EU) 2024/3190
of 19 December 2024 PROHIBITS BPA in food contact materials and repeals
Regulation (EU) 2018/213, the instrument that limit came from. There is nothing
left to back-convert, so G14 refuses the T2 route and the assessment falls back
to the TTC ladder.

That is the shape of the problem worth understanding: a regulatory database is
a dated snapshot, and a STALE limit is more dangerous than a missing one,
because it back-converts into a threshold that reads as authoritative. The
engine would rather drop a tier than quote a withdrawn number.

Cosmetic packaging is governed by Regulation (EC) 1223/2009, a different regime
again. Establish the applicable limit for the actual use before quoting one.""")

print("\nThe branch, offered because the case FAILED:")
print(q.refine(a_bpa))


# %% 4 | Refining exactly one rung ---------------------------------------------
# The measured entry above assumed the product already carries 25 mg/kg. The
# predictive entry asks a different question: given a recycled PET wall, how much
# CAN reach the product? M0 is the total-transfer bound; M3 is the kinetics.
# M0 must bound M3 -- if it does not, something is wrong with the setup.

rule("4 | The predictive entry -- and why M0 usually ends it")

print("""Step 3 assessed a product that had been ANALYSED: 25 mg/kg was in the
bottle already. This step asks the other question -- given a recycled PET wall,
how much CAN arrive? That is a different scenario, not a refinement of the same
one, and it is the question a formulator actually faces.
""")

bottle = Packaging3D('bottle',
                     body_radius=(3, 'cm'), body_height=(15, 'cm'),
                     neck_radius=(1.2, 'cm'), neck_height=(2, 'cm'))
C0 = recyclate_C0('PET', decontamination=0.90)
print("recyclate convention : C0 = %.3g mg/kg after 90%% decontamination" % C0)

# G20 -- a cosmetic swells the wall it touches, so the PLASTICIZED state is
# the conservative default and the glassy one is what would need arguing.
wall = wPET(l=(500, 'um'), C0=C0)
r = bodylotion()
bpa % r

CF_M0 = r.totaltransfer(wall, strict=False)
print("M0  CF = %.4g mg/kg   (every molecule ends in the product)" % scalar(CF_M0))
a_M0 = r.assess(CF=scalar(CF_M0))
print("    verdict: %s (margin %.3g)" % (a_M0.verdict, a_M0.margin))

# The standard SFPPy pipeline. `>> r` runs the solver and sets r.lastsimulation,
# which assess() then picks up ahead of CF0.
s = bodylotion()
bpa % s << bottle >> wall >> s
CF_M3 = scalar(s.lastsimulation.CF)
print("\nM3  CF = %.4g mg/kg   (Patankar solver, real contact conditions)"
      % CF_M3)
a_M3 = s.assess()
print("    verdict: %s (margin %.3g)" % (a_M3.verdict, a_M3.margin))
print("    provenance: %s" % a_M3.provenance)

print("\nM0 bounds M3: %s"
      % ("yes" if scalar(CF_M0) >= CF_M3 else "NO -- investigate"))

# The same wall, now against the UNIDENTIFIED population. Nothing about the
# transfer changes -- only what we are willing to claim about the migrant.
anon2 = bodylotion().declare_unidentified()
CF_anon = anon2.totaltransfer(wall, strict=False)
a_anon = anon2.assess(CF=scalar(CF_anon))
print("\nSame wall, same M0 bound, but for an UNIDENTIFIED substance:")
print("    CF = %.4g mg/kg against a T0 limit of %.3g -> %s"
      % (scalar(CF_anon), a_anon.threshold.value, a_anon.verdict.upper()))
print("""
    That is the useful result. A 90%-decontaminated recycled PET wall clears
    a NAMED substance comfortably, and FAILS the bound for the population
    nobody named. The material did not change; the claim did.

    The exit is the branch again -- identify substances to climb the T axis,
    or decontaminate further. Not "run a finer transfer model": the transfer
    was never the binding constraint.""")

print("""
Read step 4 carefully, because it is the whole doctrine in one result:

  M0 -- the crude, one-line, total-transfer bound -- ALREADY PASSED for BPA.

So M3 was never needed. It was run here only to show that the bound holds.
In practice the assessment stops at M0 and the solver is never called: gate
G8 would refuse the refinement anyway, because nothing failed.

Running M3 regardless would not have made the answer stronger. It would have
replaced a bound that needs no parameters with an estimate that needs D, k and
contact conditions -- three more things to defend, for a conclusion that was
already reached. Sophistication is not free.""")


# %% 5 | Turning the assessment around -----------------------------------------
# A compliance verdict answers "is this safe?". A design answer tells the
# formulator what to change. Both come from the same objects.

rule("5 | The design answer -- what may the material contain?")

d = bodylotion()
bpa % d
Cmax = d.maxconcentration(strict=False)
print("Maximum concentration of BPA tolerable IN THE PRODUCT:")
print("   Cmax = %.4g mg/kg" % scalar(Cmax))

# Back out the corresponding wall concentration, then the decontamination the
# recycling stream must achieve. PET density ~1350 kg/m3, wall 500 um.
rho_F = scalar(d.density)
V_F = scalar(d.volume)
A = scalar(d.surfacearea) if d.surfacearea is not None else scalar(bottle.surfacearea)
C0max = scalar(Cmax) * rho_F * V_F / (1350.0 * A * 500e-6)
print("   -> equivalent wall load C0max = %.4g mg/kg" % C0max)

eta = required_decontamination(C0max, "PET")
print("   -> required decontamination: %s"
      % ("none needed" if eta == 0 else "eta_min = %.3g %%" % (100 * eta)))
print("\nThat is a specification a recycling supplier can be held to.")


# %% 6 | The gates -------------------------------------------------------------
# A gate is the engine declining to return a number it cannot defend. Each names
# the rule it enforces and the way out. They are the reason a result can be put
# in front of a regulator.

rule("6 | Gates -- refusals, not approximations")


def show(label, fn):
    try:
        fn()
        print("   %-24s NOT FIRED (unexpected)" % label)
    except cosmeticgate as e:
        print("   %-24s %s" % (label, str(e).splitlines()[0]))


g = bodylotion()
bpa % g
g.CF0 = 1.0

show("G1  unrecognised unit", lambda: g.assess(CF=1.0, CFunit="mg per kilo",
                                               strict=True))
show("G4  censoring convention", lambda: g.assess(censoring="half", strict=True))
show("G7  local route", lambda: g.assess(route="local", strict=True))
show("G11 unknown family", lambda: recyclate_C0("PVC"))
show("G11 decontamination", lambda: recyclate_C0("PET", decontamination=90))

print("\nEach refusal is cheaper than the wrong number it prevented.")


# %% Closing -------------------------------------------------------------------

rule("Summary")
print("""
  censoring  a non-detect is data, not a number: DL/2 or DL   (step 1)
  T1         identification alone can settle a case           (step 2)
  T2         the binding tier is the one that binds           (step 3)
  M0         the crude bound passed, so M3 was never needed   (step 4)
  design     invert the question to get a specification       (step 5)
  gates      refuse rather than approximate                   (step 6)

  The result to remember from step 4: a one-line total-transfer bound
  settled the case. Every further decimal would have cost parameters to
  defend and bought no confidence.

  Cite CosPaTox (April 2024) for the exposure concepts:
  https://cospatox.com/publication/
""")
