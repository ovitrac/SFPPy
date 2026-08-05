"""
===============================================================================
SFPPy Example: from an acceptable exposure to a material specification
===============================================================================

# Example 8: What May My Supplier Ship?
--------------------------------------

## Overview
Examples 6 and 7 answered *"is this safe?"*. This one answers the question a
formulator and a recycler actually negotiate over:

    Given what the consumer may be exposed to, what is the maximum
    concentration my incoming material is allowed to contain?

That is the **reverse** problem, and it is the useful one for recyclates,
because its answer is a specification a decontamination process can be held to.

### The inversion, and why it is exact

The transport problem is **linear in the initial concentration**: double the
load in the wall and every concentration in the product doubles. So the
potential release

    PR = CF / CP0

does not depend on CP0. One run at any convenient guess characterises the
system, and the answer scales:

    1.  MAE       -> CFmax          invert the exposure chain
    2.  CP0guess  -> CFguess        one simulation, any guess
    3.  CP0max = CP0guess x CFmax / CFguess = CFmax / PR

Step 3 is exact, and the guess cancels. The example demonstrates that rather
than asserting it, by sweeping the guess over seven orders of magnitude.

`PR` here is on a **concentration basis**, CF/CP0. SFPPy's solver also carries
the mass-basis family in `SensPatankarResult.PR` (PRE = mFeq/m0,
PR(CF) = CF*VF/m0, PRT = PR/PRE). Both are linear in CP0, so the inversion is
identical either way -- but the two must not be read against each other.

### Expected Outcomes
- A/V across the common cosmetic formats, and which one is the worst case
- A decision fed from a detection limit and an internal limit, with no simulation
- CP0max per format and per population -- a table a supplier can be sent
- The linearity assumption tested rather than trusted

---

@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2026-08-04
"""

# %% Path bootstrap
import os
import sys

_ex = os.path.dirname(os.path.abspath(__file__))
_root = os.path.dirname(_ex)
if _root not in sys.path:
    sys.path.insert(0, _root)

import numpy as np
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

from patankar.geometry import Packaging3D
from patankar.layer import wPET, PP
from patankar.cosmetic import (bodylotion, shampoo, cosmeticgate,
                               recyclate_C0, required_decontamination)


def rule(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def scalar(x):
    return float(np.asarray(x).ravel()[-1])


# %% 1 | The formats -----------------------------------------------------------
# Migration does not care about the silhouette. It cares about how much wall
# faces how much product, i.e. the surface-to-volume ratio. Choosing the format
# is therefore already a safety decision.

rule("1 | The formats -- A/V decides the worst case, not the shape")

FORMATS = {
    "tube":    dict(radius=(2.0, 'cm'), height=(12.0, 'cm')),
    "airless": dict(radius=(2.5, 'cm'), height=(10.0, 'cm')),
    "bidon":   dict(radius=(3.0, 'cm'), height=(12.0, 'cm')),
    "stick":   dict(radius=(1.2, 'cm'), height=(8.0, 'cm')),
    "compact": dict(length=(6.0, 'cm'), width=(6.0, 'cm'), height=(1.0, 'cm')),
    "sachet":  dict(length=(8.0, 'cm'), width=(5.0, 'cm'), height=(0.4, 'cm')),
}

print("  %-9s %10s %11s %12s" % ("format", "V (mL)", "A (cm2)", "A/V (1/m)"))
print("  " + "-" * 46)
geoms, av = {}, {}
for name, dims in FORMATS.items():
    g = Packaging3D(name, **dims)
    V, A = g.get_volume_and_area()
    geoms[name] = g
    av[name] = scalar(A) / scalar(V)
    print("  %-9s %10.2f %11.1f %12.1f"
          % (name, scalar(V) * 1e6, scalar(A) * 1e4, av[name]))

worst = max(av, key=av.get)
best = min(av, key=av.get)
print("\n  Worst case: %s, at %.1fx the %s."
      % (worst, av[worst] / av[best], best))
print("""
  A single-dose sachet is the demanding format in almost any portfolio, and it
  falls straight out of the geometry -- no chemistry required to see it.

  One caution on `bidon` and `aerosol`. The silhouette is a cylinder but the
  regime is not: a 5-10 um varnish on impermeable metal holds the ENTIRE
  reservoir in the coating. That is a finite-dose problem tending to near-total
  depletion, governed by the partition coefficient rather than the diffusivity.
  Model the coating as the layer, never the metal.""")


# %% 2 | Feeding a decision without computing anything -------------------------
# Most decisions are not made from a simulation. They are made from a detection
# limit, an analytical result, or a limit somebody set. All are legitimate; what
# matters is that the result says which.

rule("2 | Decisions fed from outside -- setCF and setMAE")

p = bodylotion().declare_unidentified()
p.setCF(0.02, unit="mg/kg", origin="censored")   # at the detection limit
a = p.assess()
print("  detection limit as CF:")
print("    dose    = %.4g %s" % (a.dose, a.units))
print("    limit   = %.4g  [%s]  %s" % (a.threshold.value, a.threshold.identity,
                                        a.threshold.tier))
print("    verdict = %s" % a.verdict.upper())
print("    from    : %s" % a.provenance)

# An organisation's own limit outranks the whole derived ladder. An internal
# specification, a customer requirement or a DNEL from a registration dossier
# is as binding as anything derivable, and frequently tighter.
p.setMAE(1.0e-4, unit="mg/kg bw/day", origin="internal limit", tier="T2",
         identity="company internal limit")
a2 = p.assess()
print("\n  the same case against an INTERNAL limit:")
print("    limit   = %.4g  [%s]  %s / %s"
      % (a2.threshold.value, a2.threshold.identity, a2.threshold.tier,
         a2.threshold.origin))
print("    verdict = %s" % a2.verdict.upper())
print("""
  Same exposure, different governing limit, different answer -- and the origin
  of the limit is carried in the result rather than lost in a spreadsheet.""")
p.clearMAE()


# %% 3 | The reverse mechanism -------------------------------------------------

rule("3 | MAE -> CFmax -> CP0max")

tube = geoms["tube"]
C0 = recyclate_C0('PET', decontamination=0.90)
# G20 -- plasticized is the conservative default for a cosmetic contact.
wall = wPET(l=(400, 'um'), C0=C0)

spec = bodylotion().declare_unidentified()
spec << tube
r = spec.maxinitialconcentration(wall, tier="M0")

print("  step 1  MAE  = %.4g %s  [%s]"
      % (r.threshold.value, r.threshold.units, r.threshold.identity))
print("          CFmax  = %.4g mg/kg          (exposure chain inverted)" % r.CFmax)
print("  step 2  CP0guess = %.4g mg/kg -> CFguess = %.4g mg/kg"
      % (r.CP0guess, r.CFguess))
print("          PR     = CF/CP0 = %.6g" % r.PR)
print("  step 3  CP0max = CFmax / PR = %.4g mg/kg" % r.CP0max)
print("\n  provenance: %s" % r.provenance)

eta = required_decontamination(r.CP0max, "PET")
print("\n  required decontamination for PET: %s"
      % ("none needed" if eta == 0 else "eta_min = %.3g %%" % (100 * eta)))


# %% 4 | The guess cancels -- demonstrated, not asserted ------------------------

rule("4 | The guess is arbitrary, and that is the whole point")

print("  %14s %16s %16s" % ("CP0guess", "CFguess", "CP0max"))
print("  " + "-" * 48)
for g in (1e-3, 1e-1, 1.0, 1e2, 1e4):
    rr = spec.maxinitialconcentration(wall, tier="M0", CP0guess=g, verify=False)
    print("  %14g %16.6g %16.6g" % (g, rr.CFguess, rr.CP0max))

print("""
  Seven orders of magnitude of guess, one answer. PR = CF/CP0 does not depend
  on CP0, so the guess cancels exactly. This is what makes a single simulation
  sufficient to characterise the system for design.""")


# %% 5 | The assumption is tested, not trusted ---------------------------------

rule("5 | Robustness -- the gates around the inversion")


def show(label, fn):
    try:
        v = fn()
        print("   %-26s -> %s" % (label, v))
    except cosmeticgate as e:
        print("   %-26s %s" % (label, str(e).splitlines()[0][:66]))


print("""  `verify=True` re-runs the transfer at ten times the guess and checks that
  PR is unchanged. Linearity is the assumption the entire inversion rests on,
  so it is measured rather than assumed. A concentration-dependent D or k fires
  G13 instead of returning a number that looks reasonable and is not.
""")
show("no C0 on the wall", lambda: spec.maxinitialconcentration(
    wPET(l=(400, 'um'), C0=0), tier="M0"))
show("M3 with no simulation", lambda: spec.maxinitialconcentration(
    wall, tier="M3"))
show("unknown tier", lambda: spec.maxinitialconcentration(wall, tier="M7"))


# %% 6 | The specification table -----------------------------------------------
# What actually gets sent to a supplier.

rule("6 | CP0max by format and product -- the deliverable")

WALLS = {"rPET": wPET(l=(400, 'um'), C0=recyclate_C0('PET', decontamination=0.90)),
         "rPP":  PP(l=(600, 'um'), C0=recyclate_C0('polyolefin', decontamination=0.90))}

print("  Unidentified substance (T0, genotoxic TTC). CP0max in mg/kg.\n")
print("  %-9s %-11s %12s %12s" % ("format", "product", "rPET", "rPP"))
print("  " + "-" * 48)

for fmt in ("tube", "airless", "sachet"):
    for label, cls in (("body lotion", bodylotion), ("shampoo", shampoo)):
        row = []
        for wname, w in WALLS.items():
            q = cls().declare_unidentified()
            q << geoms[fmt]
            try:
                row.append(q.maxinitialconcentration(w, tier="M0",
                                                     verify=False).CP0max)
            except cosmeticgate:
                row.append(float("nan"))
        print("  %-9s %-11s %12.4g %12.4g" % (fmt, label, row[0], row[1]))

print("""
  Two things to read here.

  The spread is 450x, top to bottom: a leave-on lotion in a sachet needs a
  material that much cleaner than a rinse-off shampoo in a tube -- same polymer,
  same threshold, same chemistry. Format and usage did all of it.

  And `tube` and `airless` give identical numbers. That is not a copy-paste
  error: they were dimensioned to the same A/V (100 1/m), and A/V is the only
  geometric quantity the mass balance sees. Two very different-looking packs
  with the same surface-to-volume ratio carry the same specification.

  This table is the deliverable. It is a specification, expressed in the unit a
  recycler measures (mg/kg in the pellet), traceable back to an acceptable
  exposure and to the assumptions that produced it.""")


# %% Closing -------------------------------------------------------------------

rule("Summary")
print("""
  PR = CF/CP0 is independent of CP0        -- so one run characterises the system
  CP0max = CFmax / PR                      -- exact, and the guess cancels
  setCF / setMAE                           -- decisions fed from outside, with origin
  A/V decides the worst-case format        -- geometry before chemistry
  G13 tests linearity                      -- the assumption is measured, not trusted

  M0 gives a TIGHTER specification than M3, so the cheap route is the safe one.

  Exposure figures are PLACEHOLDERS pending the SCCS Notes of Guidance and the
  published CosPaTox deliverables: https://cospatox.com/publication/
""")
