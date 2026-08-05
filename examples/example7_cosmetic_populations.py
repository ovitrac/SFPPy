"""
===============================================================================
SFPPy Example: populations, product types, and one comparable unit
===============================================================================

# Example 7: Who Uses It, How, and Against Which Threshold
---------------------------------------------------------

## Overview
Example 6 assessed one product for one population. This one varies both, which
is the question a safety assessor actually has: **the same substance at the same
concentration is not the same risk** for an adult using a leave-on lotion and a
toddler being washed with a rinse-off gel.

Everything is reduced to a single comparable quantity:

    mg of substance per kg of body weight per day        [mg/kg bw/day]

That is the unit of a TTC and of a TDI, so exposure and threshold can be divided.
Getting every input into it is most of the work, and it is where the mistakes
live.

### The three conversions that matter

1. **Exposure -> dose.** The product mass reaching the skin per day, times the
   concentration, times dermal absorption, divided by body weight:

       E = CF x (amount x frequency x retention) x absorption / bodyweight

   `retention` is what separates leave-on (1) from rinse-off (~0.01). It is not
   a dilution of the migrant in the formulation -- it is the fraction of applied
   product that stays on the skin.

2. **SML -> TDI.** A specific migration limit is a concentration in *food*,
   derived from a dose through a convention. Running it backwards recovers the
   dose, and the convention must be stated:

       TDI = SML x (1 kg food per day) / (60 kg adult)

   `patankar.cosmetic` applies exactly this, and every threshold built that way
   carries the convention in its note.

3. **TTC -> nothing.** A TTC is already mg/kg bw/day. No conversion, no
   convention, no assumption about who is exposed.

### The asymmetry worth understanding

The 60 kg in conversion 2 is **not** an assumption about the person using the
cosmetic. It only *undoes* the food convention buried inside the SML. What comes
out is a dose per kg of body weight, which is intrinsically population-neutral --
so applying it to a 10 kg toddler is legitimate.

The population enters on the **exposure** side, through body weight, applied
amount and exposed area. Keeping the two sides separate is what makes the
comparison honest.

## Expected Outcomes
- A population x product matrix of doses, all in mg/kg bw/day
- The same matrix as margins against T0, T1 and T2 thresholds
- A demonstration that the binding constraint can be the population, not the
  chemistry

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

from patankar.loadpubchem import migrantToxtree
from patankar.cosmetic import bodylotion, shampoo, cosmeticgate


def rule(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def scalar(x):
    return float(np.asarray(x).ravel()[-1])


# %% The scenario grid ---------------------------------------------------------
#
# !! THESE FIGURES ARE ILLUSTRATIVE PLACEHOLDERS.
#
# The module's own class defaults are declared as placeholders pending the
# published CosPaTox deliverables and the SCCS Notes of Guidance, and so are the
# variations below. They are set here explicitly, per instance, precisely to show
# that they are inputs to be sourced -- not constants of nature. Replace them
# with the values of the applicable guidance before any external use.

POPULATIONS = {
    #  label        bodyweight  exposed area (whole body, applied fraction)
    "adult":    dict(bodyweight=(60.0, "kg"), exposedarea=(16000.0, "cm**2")),
    "child":    dict(bodyweight=(20.0, "kg"), exposedarea=(8000.0, "cm**2")),
    "toddler":  dict(bodyweight=(10.0, "kg"), exposedarea=(5000.0, "cm**2")),
    "infant":   dict(bodyweight=(5.0, "kg"),  exposedarea=(3000.0, "cm**2")),
}

# Amount applied scales roughly with body surface, so it is given per population
# too. Retention is the leave-on / rinse-off discriminator.
AMOUNT_G = {"adult": 8.0, "child": 4.0, "toddler": 2.5, "infant": 1.5}

USAGES = {
    "leave-on":  dict(cls=bodylotion, retention=1.00),
    "rinse-off": dict(cls=shampoo,    retention=0.01),
}


def build(population, usage, absorption=1.0):
    """A product instance parameterised for one population and one usage."""
    spec = dict(POPULATIONS[population])
    spec["amountapplied"] = (AMOUNT_G[population], "g")
    spec["retention"] = (USAGES[usage]["retention"], None)
    spec["dermalabsorption"] = (absorption, None)
    return USAGES[usage]["cls"](**spec)


# %% 1 | The exposure side -----------------------------------------------------

rule("1 | Exposure -- the same concentration, eight different doses")

CF = 1.0                                       # mg/kg of substance in the product
print("A single substance at CF = %.3g mg/kg in the product." % CF)
print("Tier 1 default: 100% dermal absorption -- deliberately penalising.\n")

print("  %-9s %-10s %10s %8s %14s" %
      ("usage", "population", "kg/day", "bw (kg)", "mg/kg bw/day"))
print("  " + "-" * 56)

doses = {}
for usage in USAGES:
    for pop in POPULATIONS:
        p = build(pop, usage)
        d = p.systemicdose(CF=CF)
        doses[(usage, pop)] = d
        print("  %-9s %-10s %10.4g %8.3g %14.4g"
              % (usage, pop, p.dailyproductmass, scalar(p.bodyweight), d))

ratio = doses[("leave-on", "infant")] / doses[("rinse-off", "adult")]
print("\n  Spread across the grid: %.0fx." % ratio)
print("  Nothing about the substance changed. Retention did most of it, body")
print("  weight the rest -- and the two compound.")


# %% 2 | The threshold side ----------------------------------------------------

rule("2 | Thresholds -- and the one conversion that carries a convention")

sub = migrantToxtree("bisphenol A")
p = sub % build("adult", "leave-on")
thr = p.binding_threshold()

print("substance : %s" % thr.identity)
print("value     : %.4g %s" % (thr.value, thr.units))
print("tier      : %s   origin: %s" % (thr.tier, thr.origin))
if thr.note:
    print("note      : %s" % thr.note)

print("""
Read the note. An SML is a limit on concentration IN FOOD; the dose it was
derived from is recovered by dividing out the food convention:

    TDI = SML x (1 kg food/day) / (60 kg adult)

The 60 kg is not a statement about who uses the cosmetic. It only undoes the
convention inside the SML. The result -- mg per kg of BODY WEIGHT per day -- is
population-neutral, which is why the same threshold can be applied to a 5 kg
infant without adjustment.

A TTC needs none of this: it is already mg/kg bw/day.""")


# %% 3 | Putting the two sides together ----------------------------------------

rule("3 | Margins -- exposure and threshold in the same unit")

# At CF = 1 mg/kg every row passes -- the grid is informative but undramatic.
# Raised to 5 mg/kg the SAME substance at the SAME concentration straddles the
# threshold, and where it falls is decided entirely by the population and usage.
CF_ASSESS = 5.0

print("substance: bisphenol A, at CF = %.3g mg/kg in the product" % CF_ASSESS)
print("threshold: %.4g mg/kg bw/day (%s)\n" % (thr.value, thr.tier))

print("  %-9s %-10s %14s %10s %9s" %
      ("usage", "population", "mg/kg bw/day", "margin", "verdict"))
print("  " + "-" * 56)

flipped = []
for usage in USAGES:
    for pop in POPULATIONS:
        q = sub % build(pop, usage)
        q.CF0 = CF_ASSESS
        a = q.assess()
        print("  %-9s %-10s %14.4g %10.3g %9s"
              % (usage, pop, a.dose, a.margin, a.verdict.upper()))
        flipped.append((usage, pop, a.verdict))

npass = sum(1 for _, _, v in flipped if v == "pass")
print("\n  %d of %d combinations pass." % (npass, len(flipped)))

print("""
  This is the point of the whole exercise. The chemistry is identical in every
  row: same substance, same concentration, same threshold. The verdict is not --
  and where it flips is decided by who uses the product and how, never by the
  transfer model.

  A safety assessment that reports one number for "the product" has silently
  chosen a population. Here the choice is explicit and auditable.

  Note which rows fail. A leave-on product on a small body weight is the
  demanding case, and no refinement of the migration model addresses it.""")


# %% 4 | Inverting it: the specification per population ------------------------

rule("4 | The design answer, per population")

print("Maximum concentration tolerable IN THE PRODUCT, mg/kg:\n")
print("  %-9s %-10s %14s %16s" %
      ("usage", "population", "Cmax (T2)", "Cmax (T0, unident.)"))
print("  " + "-" * 56)

for usage in USAGES:
    for pop in POPULATIONS:
        named = sub % build(pop, usage)
        anon = build(pop, usage).declare_unidentified()
        print("  %-9s %-10s %14.4g %16.4g"
              % (usage, pop, scalar(named.maxconcentration()),
                 scalar(anon.maxconcentration())))

print("""
  The right-hand column needs no substance at all -- it is the genotoxic TTC
  applied to the population nobody named. For a leave-on product on an infant it
  is the tightest specification in the table, and it is the one a recyclate has
  to meet before any identification work begins.""")


# %% 5 | A guard worth seeing --------------------------------------------------

rule("5 | Routes are not interchangeable")

try:
    (sub % build("infant", "leave-on")).assess(route="local", strict=True)
except cosmeticgate as e:
    print("   %s" % str(e).splitlines()[0])
    print("""
   Systemic exposure is mg/kg bw/day; a site-of-contact load is ug/cm2. They
   answer different questions and are never coerced into one another. Absence of
   a local threshold is an explicit state, not a permissive default.""")


# %% Closing -------------------------------------------------------------------

rule("Summary")
print("""
  exposure   E = CF x (amount x frequency x retention) x absorption / bw
  SML -> TDI divide out the food convention: 1 kg food/day, 60 kg adult
  TTC        already mg/kg bw/day -- no convention, no assumption

  Population enters through EXPOSURE. The threshold, once expressed per kg of
  body weight, is population-neutral. Keeping those two apart is what makes the
  comparison defensible.

  The exposure figures used here are PLACEHOLDERS. Source them from the SCCS
  Notes of Guidance and the published CosPaTox deliverables before external use:
  https://cospatox.com/publication/
""")
