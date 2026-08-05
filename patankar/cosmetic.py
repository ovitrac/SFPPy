#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
SFPPy Module: Cosmetic Exposure
===============================================================================
Defines **cosmetic and home-care products** as contacting media, and adds the
**exposure stage** that food never needed.

For food, EU practice folds the whole exposure chain into the limit itself:
6 dm²/kg food, 1 kg food per day, 60 kg body weight. The SML *is* the exposure
model, so a simulation may stop at `CF`. Regulation (EC) 1223/2009 carries no
such convention for cosmetic packaging, so the chain must be modelled
explicitly, per product category, and it terminates in a **dose**, not a
concentration:

    C_product --(amount applied, retention)--> dose --(area, absorption)-->
        systemic [mg/kg bw/day] -> DNEL / TTC
        local    [ug/cm2]       -> site-of-contact threshold

**Main Components:**
- **Base class `cosmetic`** (sibling of `realfood` and `simulant` under
  `foodproperty`) carrying the exposure chain and the assessment.
- **Categories** `bodylotion`, `shampoo`, `washinggel` -- each inheriting its
  prescribed simulant, so the simulant doctrine is expressed by the type and
  `_chemicalsubstance` is supplied for the `%` operator.
- **Gates** (`cosmeticgate`): enforced refusals that stop a computation rather
  than let it return a plausible wrong number. See `help_gates()`.
- **Recyclate helpers** `recyclate_C0` and `required_decontamination`.

**Two concentration entries, mutually exclusive:**
1. *Measured* -- the product was analysed. `CF0` is the measured value, or
   half the detection limit for a non-detect. No geometry, no layer, no solver.
2. *Predictive* -- `C0` is set on the packaging layer and the transfer is a
   mass balance over the real volumes (M0) or the Patankar solver (M3).

There is deliberately **no dilution factor**. It is not a parameter but a
consequence of geometry, wall thickness and densities, all of which SFPPy
already carries.

**Integration with SFPPy Modules:**
- Inherits the whole transfer machinery from `food.py` unchanged.
- `substance % product` works by inheritance (`food.foodphysics.__rmod__`).
- `product << packaging` inherits volume and surface area from `geometry.py`.
- Substances come from `loadpubchem` -- `migrant` or `migrantToxtree` -- and
  never as duck-typed objects: `migrant` owns identity (G12).

**Advanced Operators** -- the standard SFPPy pipeline applies unchanged:

    substance % product << geometry >> packaging >> product

`assess()` then reads the concentration in the order the pipeline produces it:
an explicit argument first, then `product.lastsimulation` (M3), then `CF0`
(the measured entry). The provenance of whichever was used is reported.

Example -- measured entry, no simulation:
```python
from patankar.loadpubchem import migrantToxtree
from patankar.cosmetic import bodylotion
m = migrantToxtree("toluene")
p = m % bodylotion()
p.CF0 = 0.5                                 # mg/kg, non-detect at DL = 1 mg/kg
p.assess()
```

Example -- predictive entry through the standard pipeline:
```python
from patankar.geometry import Packaging3D
from patankar.layer import wPET
from patankar.cosmetic import bodylotion, recyclate_C0
bottle = Packaging3D('bottle', body_radius=(3,'cm'), body_height=(15,'cm'))
# the medium swells the wall, so the PLASTICIZED state is the default (G20)
wall = wPET(l=(500,'um'), C0=recyclate_C0('PET', decontamination=0.90))
p = bodylotion()
m % p << bottle >> wall >> p                # runs; p.lastsimulation is set
p.assess()
```

**Provenance.** The concepts transposed here -- exposure classes, the simulant
doctrine, the three worked product categories and the Tier 1 defaults -- come
from the CosPaTox consortium's published work. Any use must cite the April 2024
deliverables (Guideline, Substance List, Scientific Dossier, WCC Calculator) at
https://cospatox.com/publication/ . No DOI and no consortium-mandated citation
format exist; cite by name, consortium, month/year and URL. The estimator layer
of the dossier is deliberately **not** transposed: it converts a test outcome
into an initial concentration and belongs to the laboratory workflow.

@version: 1.0
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2026-08-03
@rev: 2026-08-03

"""

# Dependencies
import numpy as np

from patankar.layer import check_units, NoUnits, layer
from patankar.food import foodproperty, ethanol, ethanol50

__all__ = ['as_migrant', 'bodylotion', 'cosmetic', 'cosmeticassessment',
           'cosmeticgate', 'help_cosmetic', 'help_gates', 'recyclate_C0',
           'required_decontamination', 'shampoo', 'threshold', 'unidentified',
           'washinggel']

__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.0"


# %% Gates -------------------------------------------------------------------

class cosmeticgate(Exception):
    """
    Raised when a computation is **refused** rather than approximated.

    A gate is not an error in the usual sense: it is the engine declining to
    return a number it cannot defend. Every gate names the rule it enforces and
    the way out. Gates may be downgraded to recorded warnings with
    ``strict=False``, in which case the assessment carries them in
    ``.gates`` and the result must be read as provisional.
    """

    def __init__(self, gate, message, remedy=None):
        self.gate = gate
        self.remedy = remedy
        full = "[%s] %s" % (gate, message)
        if remedy:
            full += "\n  -> %s" % remedy
        super().__init__(full)


GATES = {
    "G1-units":
        "A threshold or concentration is refused unless its unit is recognised. "
        "Source data carries a dozen spellings across a 1000x mg/ug divide; a "
        "silent numeric read mixes two comparable populations.",
    "G2-validity":
        "A substance without a molecular weight has no transfer model. UVCBs, "
        "polymers and mixtures are reported not-assessable, never estimated.",
    "G3-kprovenance":
        "A partition coefficient inherited from an affinity mixin is a taxonomy "
        "placeholder, not a measurement. Injecting the substance with % assigns "
        "the real value; falling back to the class default is recorded.",
    "G4-censoring":
        "A '< x' observation is left-censored data, not a number. It enters as "
        "DL/2 for an identified substance, or DL for the T0 bound on an "
        "uncharacterised one. Halving a bound is not halving an estimate.",
    "G5-ttc":
        "Where a source carries an adjudicated TTC, it wins over one recomputed "
        "from the Cramer class. The genotoxic override has already been applied "
        "upstream; recomputing can be 12000x permissive.",
    "G6-massbasis":
        "Every concentration declares a mass or volume basis. The solver never "
        "converts, so mg/kg in is mg/kg out only when the two densities agree.",
    "G7-route":
        "Systemic [mg/kg bw/day] and local [ug/cm2] are not commensurable and "
        "are never coerced into one another. Absence of a local threshold is an "
        "explicit state, not a default.",
    "G8-stoponpass":
        "Refinement without a trigger is a defect, not diligence. A tier that "
        "passes is not refined.",
    "G9-polymerstate":
        "The r- prefix in SFPPy means RUBBERY, not recycled. Using rPET below "
        "its glass transition models the polymer in the wrong state.",
    "G17-reactiondirection":
        "Assuming a reaction is complete is conservative for the PRODUCTS it "
        "forms and anti-conservative for the PARENT it consumes. An assessment "
        "that passes because its own analyte was destroyed is not an "
        "assessment. Complete conversion is applied only in the direction that "
        "raises exposure.",
    "G18-partitionbound":
        "Unitary partition is a convenience; total transfer is the guarantee. "
        "K=1 is not the worst case, and how much conservatism it gives away is "
        "decided by the polymer-to-product volume ratio -- negligible for a "
        "bottle, a third of the bound for a thick-walled stick or compact.",
    "G19-controlclass":
        "A parameter the customer controls -- design, contact time, contact "
        "temperature -- was left at a class default. Bounding a controlled "
        "parameter is not caution but a category error: it inflates the answer "
        "with a margin a specification sheet would have removed, and it buries "
        "the parameters that actually drive the result. Ask, do not guess.",
    "G21-tgnoop":
        "A declared glass transition that the selected diffusion model does "
        "not consume changes nothing. Recording an assumption the arithmetic "
        "ignores is worse than not declaring it: the result looks like it was "
        "bounded and was not.",
    "G20-plasticization":
        "A cosmetic swells the wall it touches. The glassy state is therefore "
        "the OPTIMISTIC one, and it is the claim that must be justified -- not "
        "the plasticized state. Assume plasticization is complete unless the "
        "medium is evidenced not to swell the polymer.",
    "G10-smlbasis":
        "An SML is back-converted to an implied TDI only when it is "
        "toxicologically derived. Group SMLs are not per-substance; QM- and "
        "analytically-derived limits carry no TDI at all.",
    "G11-recyclate":
        "A recyclate convention applies to a named polymer family at a "
        "declared decontamination efficiency. Neither is guessed.",
    "G15-nodata":
        "Absence of evidence is not evidence of absence. Nothing measured, "
        "nothing simulated and nothing supplied is the case that most deserves "
        "refusal, not the one that quietly passes because 0 > threshold is "
        "false. A genuine zero is a non-detect and belongs at DL or DL/2 (G4).",
    "G16-doublecensoring":
        "Two censoring conventions compose silently. A value entered "
        "pre-halved and then passed through censoring='DL/2' is halved twice: "
        "a dose four times too low, with a provenance that reads correctly. "
        "Store the raw detection limit; let assess() apply the convention once.",
    "G14-regstatus":
        "The shipped regulatory databases are dated snapshots. Where an "
        "instrument has since changed a substance's status, the stored limit "
        "is stale -- and a stale SML is worse than a missing one, because it "
        "back-converts into a TDI that reads as authoritative.",
    "G13-linearity":
        "Inverting an exposure into a material specification assumes CF is "
        "linear in CP0, so that PR = CF/CP0 is independent of the guess. A "
        "concentration-dependent D or k breaks it. The assumption is tested "
        "by rescaling the guess, not trusted.",
    "G12-substance":
        "A substance is consumed through loadpubchem -- migrant or "
        "migrantToxtree -- never as a duck-typed object. migrant owns identity: "
        "synonyms, registry numbers, molecular weight, the databases.",
}


# %% The unidentified substance ----------------------------------------------

class _unidentified:
    """
    The uncharacterised population -- an explicit declaration of *no identity*.

    In a recyclate the commonest case is the substance nobody named: detected
    below the limit, or not detected at all, with no structure, no CAS and no
    Cramer class. Tier T0 exists for exactly this population -- the genotoxic
    TTC applies to *any* molecule, so a single bound covers all of them without
    naming one.

    Attaching this sentinel says so **explicitly**. It is not a substance and
    never pretends to be one:

    - it carries no molecular weight, so it cannot enter a transfer model;
    - `migrant.__rmod__` rejects it, so `unidentified % product` raises and the
      sentinel can never reach the Patankar solver by accident -- tier M3 is
      structurally unreachable, which is correct: there is no `D` and no `k`
      for a molecule with no structure;
    - the M0 mass balance and the whole dose chain accept it, because neither
      reads the substance. `totaltransfer` is geometry and densities;
      `systemicdose` is amount, retention, frequency and body weight.

    Usage -- assignment, deliberately not the `%` operator, so that declaring
    the absence of identity never looks like injecting a substance::

        p = bodylotion()
        p.substance = unidentified
        p.CF0 = 1.0                        # at the detection limit
        p.assess(censoring="DL")           # T0, no molecule named

    G4 applies: a non-detect on an uncharacterised substance is counted at
    **DL**, not DL/2. Halving a bound is not the same operation as taking the
    midpoint of a censored distribution.

    Leaving `substance` as ``None`` remains a refusal (G2). A forgotten
    injection must not silently become a T0 answer; the declaration has to be
    made on purpose.
    """

    __slots__ = ()

    M = None                    # no molecular weight: no transfer model
    compound = "unidentified substance"
    CAS = None
    TTC = None                  # no adjudicated value
    CramerValue = None          # no Cramer class -> _ttc_of returns T0
    hasSML = False
    hasSMLgroup = False

    def __repr__(self):
        return "<unidentified substance -- T0, genotoxic TTC>"

    __str__ = __repr__


#: Singleton marking the uncharacterised population. See :class:`_unidentified`.
unidentified = _unidentified()


def _is_unidentified(sub):
    """True when `sub` is the explicit no-identity declaration."""
    return isinstance(sub, _unidentified)


def _identifier(sub):
    """
    Display identity, delegated to `migrant`.

    Chemical identity is `migrant`'s responsibility, not this module's:
    `migrant` resolves the query, owns the synonym list and normalises the
    registry numbers. So `compound` -- the resolved canonical name -- is read
    directly, and truncation goes through the object's own `dispmax`. Nothing
    here parses synonyms; doing so would duplicate (and eventually contradict)
    `loadpubchem`.

    One caution on `CAS`. It is a **list built from a set**, so its order is
    not stable across processes -- Python randomises string hashing per run,
    and toluene renders as either `108-88-3` or `25013-04-1` first depending on
    the seed. Indexing it would make the displayed identity vary between runs
    of the same script. Every registry number is therefore shown, in sorted
    order. The sort is a rendering choice for reproducibility; it asserts no
    ranking, because which number is principal is `loadpubchem`'s business.
    """
    if _is_unidentified(sub):
        return "unidentified substance"
    nm = getattr(sub, "compound", None)
    if not (isinstance(nm, str) and nm.strip()):
        return str(sub)                      # migrant renders itself
    nm = sub.dispmax(nm.strip(), 32) if hasattr(sub, "dispmax") else nm.strip()
    cas = getattr(sub, "CAS", None)
    if cas:
        items = sorted(cas) if isinstance(cas, (list, tuple, set)) else [str(cas)]
        return "%s (%s)" % (nm, ", ".join(str(c) for c in items))
    return nm


def as_migrant(query, toxtree=True):
    """
    Resolve a query into a `migrant` / `migrantToxtree`, or say why not.

    **Substances are consumed through `loadpubchem` only.** `migrant` manages
    identity -- synonyms, registry numbers, molecular weight, the databases --
    and this module never reconstructs any of it.

    `migrant()` *raises* `ValueError` for a query it cannot resolve, so the
    UVCB / polymer / mixture case is caught at **construction**, before an
    object exists to assess. This helper turns that into a reportable outcome
    rather than an exception in the middle of a sweep.

    Returns (substance, reason) with exactly one of the two non-None.
    """
    from patankar.loadpubchem import migrant, migrantToxtree
    cls = migrantToxtree if toxtree else migrant
    try:
        return cls(query), None
    except Exception as exc:
        return None, ("%s could not be resolved by loadpubchem (%s). It is "
                      "typically a UVCB, a polymer or a mixture: no single "
                      "structure, hence no molecular weight and no transfer "
                      "model." % (query, type(exc).__name__))


def help_gates():
    """Print the gate catalogue."""
    print("SFPPy cosmetic gates -- refusals, not approximations\n")
    for k in sorted(GATES):
        print("%-16s %s\n" % (k, GATES[k]))


def _raise(gate, message, remedy=None, strict=True, log=None):
    """Fire a gate: raise when strict, otherwise record and continue."""
    exc = cosmeticgate(gate, message, remedy)
    if strict:
        raise exc
    if log is not None:
        log.append(str(exc))
    return exc


# %% Units and provenance ----------------------------------------------------

# Accepted spellings for a systemic (body-weight) dose, normalised to mg/kg/day.
_SYSTEMIC_UNITS = {
    "mg/kg bw/day": 1.0, "mg/kg bw/d": 1.0, "mg/kgbw/day": 1.0,
    "mg/kg/day": 1.0, "mg/kg/d": 1.0,
    "ug/kg bw/day": 1e-3, "ug/kg bw/d": 1e-3, "ug/kgbw/day": 1e-3,
    "ug/kg/day": 1e-3, "ug/kg/d": 1e-3,
    "µg/kg bw/day": 1e-3, "µg/kg bw/d": 1e-3, "µg/kgbw/day": 1e-3,
    "µg/kg/day": 1e-3, "µg/kg/d": 1e-3,
}

# Accepted spellings for a local (site-of-contact) load, normalised to ug/cm2.
_LOCAL_UNITS = {
    "ug/cm2": 1.0, "µg/cm2": 1.0, "ug/cm²": 1.0, "µg/cm²": 1.0,
    "mg/cm2": 1e3, "mg/cm²": 1e3,
    "ng/cm2": 1e-3, "ng/cm²": 1e-3,
}

# Accepted spellings for a concentration in the product, normalised to mg/kg.
_CONC_UNITS = {
    "mg/kg": 1.0, "ppm": 1.0, "mg/g": 1e3, "ug/g": 1.0, "µg/g": 1.0,
    "g/kg": 1e3, "%": 1e4, "ug/kg": 1e-3, "µg/kg": 1e-3, "ppb": 1e-3,
}

# %% Polymer state ------------------------------------------------------------

# Glassy class -> its PLASTICIZED (swollen) sibling. The `w` prefix in SFPPy is
# the general swollen state, not "wet": whatever sorbs -- water, ethanol,
# surfactant, an oxidising solution -- depresses Tg and raises D. For a cosmetic
# the swollen state is the conservative default and the glassy one is the claim
# that must be argued (G20).
#
# Only pairs that SFPPy actually parameterises are listed. A polymer with no
# plasticized sibling does not fire G20: there is nothing to steer the user to,
# and inventing one would be worse than the silence.
PLASTICIZED_SIBLING = {
    "gPET": "wPET",
    "rigidPVC": "plasticizedPVC",
}

# Rubbery class -> the glassy class of the same polymer, used only to name the
# right alternative in G9's remedy.
GLASSY_OF = {
    "rPET": "gPET",
}

# %% Parameter control -------------------------------------------------------
#
# Provenance answers *where a number came from*. It does not answer *who owns
# it*, and the two demand opposite treatments.
#
#   CONTROLLED    the customer decides it and can evidence it -- package design,
#                 wall, contact time, contact temperature. It must be REALISTIC:
#                 specifically, the realistic worst case over the declared
#                 specification envelope. Bounding it is not caution, it is a
#                 category error -- it inflates the answer with a margin a
#                 specification sheet would have removed, and it buries the
#                 parameters that actually drive the result.
#
#   UNCONTROLLED  nobody measured it for this case -- D, k, identity, threshold,
#                 plasticization state, extent of reaction. It must be BOUNDED,
#                 and the direction of the bound must be stated.
#
# Hence the asymmetry enforced by G19: a defaulted CONTROLLED parameter is a
# defect, and the remedy is to ask the customer, not to guess.

CONTROLLED = "controlled"
UNCONTROLLED = "uncontrolled"

#: parameter -> (control class, what to ask for when it is defaulted)
PARAM_CONTROL = {
    "contacttime": (
        CONTROLLED,
        "the shelf life actually specified -- period-after-opening plus storage "
        "-- not a test convention"),
    "contacttemperature": (
        CONTROLLED,
        "the storage and transport temperature actually specified, including "
        "the excursion rather than the warehouse mean"),
    "surfacearea": (
        CONTROLLED,
        "the contact area of the real package; propagate it with "
        "`product << packaging3D` instead of accepting the class default"),
    "volume": (
        CONTROLLED,
        "the fill volume of the real package"),
    "amountapplied": (UNCONTROLLED, "consumer behaviour; bound it"),
    "frequency": (UNCONTROLLED, "consumer behaviour; bound it"),
    "retention": (UNCONTROLLED, "consumer behaviour; bound it"),
    "exposedarea": (UNCONTROLLED, "consumer behaviour; bound it"),
    "bodyweight": (UNCONTROLLED, "population parameter; bound it"),
    "dermalabsorption": (UNCONTROLLED,
                         "CosPaTox Tier 1 bounds it at 1 (100 % penetration)"),
}

#: Layer classes that are PET in one state or another. `boundingD` is scoped to
#: these, because the models it reconciles are PET parameterisations.
PET_FAMILY = ("gPET", "wPET", "rPET", "PET")


class _TProbe:
    """Minimal stand-in carrying only `T`, so `boundingD` can be asked about one
    layer of a multilayer without constructing a whole layer object."""
    __slots__ = ("T",)

    def __init__(self, T):
        self.T = T


#: Parameters a geometry supplies through `product << packaging3D`.
_GEOMETRY_PARAMS = ("surfacearea", "volume")

# How much of the total-transfer bound unitary partition may give away before
# it stops being a usable substitute for it (G18). Measured across cosmetic
# formats: a 400 mL bottle at 500 um gives away 4 %, a 5 mL sachet 9 %, a 30 mL
# jar with a 2 mm wall 29 %, a compact two thirds. The line falls naturally
# between the thin-walled, high-volume formats and the thick-walled ones.
UNITARY_PARTITION_MAXLOSS = 0.10

# The EU food-contact ACCELERATED TEST convention, which `foodphysics` supplies
# as its default contact condition. It is a test protocol, not a cosmetic shelf
# life, and inheriting it silently is the single most likely way for a cosmetic
# assessment to be built on the wrong contact conditions -- 10 days is three
# orders of magnitude short of a 36-month period-after-opening.
_TESTCONTACT_TIME_S = 864000.0      # 10 days
_TESTCONTACT_TEMP_C = 40.0

# Provenance of a numeric input. An M3 result on a defaulted parameter is
# weaker than the M0 bound it replaced -- the tag travels with the number.
MEASURED = "measured"
ADJUDICATED = "adjudicated"     # taken from a source that already decided
ADJUDICATED_FOOD = "adjudicated (food)"  # decided for INGESTION, proxy here (A5)
ESTIMATED = "estimated"         # from a model (Flory-Huggins, Piringer, ...)
DEFAULTED = "defaulted"         # a class default; the weakest tag
CENSORED = "censored"           # derived from a detection limit
SIMULATED = "simulated"         # from the Patankar solver (M3)


def _normalise_unit(u):
    """Collapse whitespace and non-breaking spaces in a unit string."""
    if u is None:
        return None
    return " ".join(str(u).replace(" ", " ").split())


def _convert(value, unit, table, gate, what, strict=True, log=None):
    """Convert `value` to the table's canonical unit, or fire the units gate."""
    if value is None:
        return None
    u = _normalise_unit(unit)
    if u is None:
        _raise(gate, "%s = %s supplied without a unit." % (what, value),
               "Supply the unit explicitly; recognised: %s"
               % ", ".join(sorted(set(table))[:6]) + " ...", strict, log)
        return None
    if u not in table:
        _raise(gate, "%s carries the unrecognised unit %r." % (what, u),
               "Normalise the spelling. Recognised: %s"
               % ", ".join(sorted(set(table))[:6]) + " ...", strict, log)
        return None
    return float(value) * table[u]


# %% Thresholds --------------------------------------------------------------

# Threshold of Toxicological Concern, ug/kg bw/day, indexed as in
# loadpubchem.migrantToxtree.TTC = [genotoxic, Cramer I, Cramer II, Cramer III]
TTC_UGKGBWDAY = (0.0025, 30.0, 9.0, 1.5)
TTC_TIER = ("T0", "T1", "T1", "T1")
TTC_NAME = ("genotoxic TTC", "Cramer class I TTC",
            "Cramer class II TTC", "Cramer class III TTC")

# The food convention hidden inside an SML: 1 kg food per day, 60 kg body
# weight. Used ONLY to run the conversion backwards (SML -> implied TDI).
FOOD_INTAKE_KGDAY = 1.0
FOOD_BODYWEIGHT_KG = 60.0

# `patankar.layer.layer` carries this as a class default and exposes no
# density(). Accepting it in the mass balance silently understates a PET wall by
# 26 %, non-conservatively -- so it is recognised and refused rather than used.
_LAYER_RHO_DEFAULT = 1000.0


class threshold:
    """
    A toxicological threshold with its identity, tier, route and provenance.

    A bare number is not enough. Which threshold was binding, how informed it
    was, and where it came from decide what the reader should do next: a
    substance failing at T0 because nothing is known about it is an
    identification problem, while one failing against its own SML-derived TDI
    is a formulation or material problem.

    Attributes
    ----------
    value    : float   -- mg/kg bw/day (systemic) or ug/cm2 (local)
    route    : str     -- "systemic" or "local"; never coerced (G7)
    identity : str     -- human-readable name of the threshold
    tier     : str     -- "T0", "T1" or "T2" on the T axis
    origin   : str     -- provenance tag
    """

    __slots__ = ("value", "route", "identity", "tier", "origin", "note")

    def __init__(self, value, route, identity, tier, origin, note=None):
        self.value = value
        self.route = route
        self.identity = identity
        self.tier = tier
        self.origin = origin
        self.note = note

    @property
    def units(self):
        return "mg/kg bw/day" if self.route == "systemic" else "ug/cm2"

    def __repr__(self):
        s = "<threshold %s %s [%s] (%s, %s)>" % (
            self.value, self.units, self.identity, self.tier, self.origin)
        return s

    __str__ = __repr__


# %% Regulatory currency -----------------------------------------------------
#
# SFPPy ships a snapshot of EU 10/2011 Annex I, the US FDA FCN inventory and the
# Chinese GB positive lists. A snapshot has a date, and regulation moves. Where
# an instrument has changed a substance's status AFTER the shipped snapshot, the
# value in the database is stale -- and a stale SML is worse than a missing one,
# because it back-converts into a TDI that looks authoritative.
#
# This table is deliberately small, explicit and auditable: one entry per
# substance whose status is known to have changed, each naming the instrument,
# its date, and what it did. It is not a regulatory database and does not try to
# be one. Extend it when an instrument is known; never infer an entry.

REGULATORY_STATUS = {
    # Bisphenol A. Commission Regulation (EU) 2024/3190 of 19 December 2024
    # prohibits BPA and its salts in the manufacture of food contact materials
    # and articles, and REPEALS Regulation (EU) 2018/213 -- the instrument that
    # carried the 0.05 mg/kg specific migration limit. The limit no longer
    # exists, so it cannot be back-converted into an implied TDI.
    "80-05-7": dict(
        status="prohibited",
        scope="EU food contact materials",
        instrument="Commission Regulation (EU) 2024/3190 of 19 December 2024",
        date="2024-12-19",
        supersedes="Regulation (EU) 2018/213 (repealed), SML 0.05 mg/kg",
        note="Limited derogations exist (polysulfone filtration membranes; "
             "large-capacity epoxy varnishes on tanks and vessels), with "
             "transition periods to 2026-07 and 2028-01. Cosmetic packaging "
             "falls under Regulation (EC) 1223/2009, a different regime: "
             "check it separately.",
    ),
}


def _regulatory_status(sub):
    """
    Look up a known post-snapshot change of regulatory status.

    Keyed on CAS. `migrant.CAS` is a list built from a set, so every entry is
    tested rather than the first one -- which is not stable across processes.
    """
    cas = getattr(sub, "CAS", None)
    if not cas:
        return None
    items = cas if isinstance(cas, (list, tuple, set)) else [cas]
    for c in items:
        hit = REGULATORY_STATUS.get(str(c).strip())
        if hit is not None:
            return hit
    return None


def _sml_to_tdi(sub, strict=True, log=None):
    """
    Recover the implied TDI from a specific migration limit.

    An SML is a concentration limit in food, derived from a toxicological dose
    through the food convention. Running that backwards is legitimate -- but
    ONLY when the limit is toxicologically derived. Group SMLs are not
    per-substance, and QM- or analytically-derived limits carry no TDI (G10).
    """
    if not getattr(sub, "hasSML", False) or sub.SML is None:
        return None

    # G14 -- the shipped database is a snapshot. If the limit has since been
    # withdrawn, there is nothing to back-convert.
    reg = _regulatory_status(sub)
    if reg is not None and reg["status"] == "prohibited":
        _raise("G14-regstatus",
               "%s is prohibited in %s by %s; the SML carried by the shipped "
               "database snapshot (%s) has been withdrawn and cannot be "
               "back-converted into an implied TDI."
               % (_identifier(sub), reg["scope"], reg["instrument"],
                  reg["supersedes"]),
               "Fall back to the TTC ladder, and establish the applicable "
               "limit under the regime that actually governs this use. "
               + reg["note"], False, log)
        return None

    if getattr(sub, "hasSMLgroup", False):
        _raise("G10-smlbasis",
               "%s is regulated in a group of %d substances; its SML is not "
               "per-substance." % (_identifier(sub),
                                   getattr(sub, "nSMLgroup", 0)),
               "Fall back to the TTC ladder and record that the SML was not "
               "usable.", strict, log)
        return None
    sml = float(np.asarray(sub.SML).ravel()[0])
    tdi = sml * FOOD_INTAKE_KGDAY / FOOD_BODYWEIGHT_KG
    # A5 -- the adjudication was performed for a FOOD INGESTION scenario and is
    # then used as a systemic limit for DERMAL exposure to a cosmetic. That is
    # defensible as a proxy and indefensible as an adjudicated limit for this
    # route, so the origin says which. One word, and it is the difference
    # between a proxy and a claim.
    return threshold(tdi, "systemic",
                     "implied TDI from SML = %g mg/kg food" % sml, "T2",
                     ADJUDICATED_FOOD,
                     note="back-converted with the food convention "
                          "(%g kg food/day, %g kg bw); adjudicated for FOOD "
                          "INGESTION and used here as a proxy for a cosmetic "
                          "route -- not a dermal adjudication"
                          % (FOOD_INTAKE_KGDAY, FOOD_BODYWEIGHT_KG))


def _ttc_of(sub, strict=True, log=None):
    """
    Systemic threshold from the TTC ladder.

    If the substance already carries an adjudicated TTC, that value wins and a
    disagreement with the Cramer class is reported rather than resolved (G5):
    upstream sources apply the genotoxic override, so a substance may be
    Cramer class I and still carry the T0 value.
    """
    ttc = getattr(sub, "TTC", None)
    cramer = getattr(sub, "CramerValue", None)
    if ttc is None:
        if cramer is None:
            return threshold(TTC_UGKGBWDAY[0] * 1e-3, "systemic",
                             TTC_NAME[0], "T0", DEFAULTED,
                             note="no Cramer class available; the unknown "
                                  "substance is treated as genotoxic")
        return threshold(TTC_UGKGBWDAY[cramer] * 1e-3, "systemic",
                         TTC_NAME[cramer], TTC_TIER[cramer], ESTIMATED)
    ttc = float(np.asarray(ttc).ravel()[0])
    tier, name, origin = "T1", "TTC (adjudicated)", ADJUDICATED
    if abs(ttc - TTC_UGKGBWDAY[0]) < 1e-12:
        tier, name = "T0", "genotoxic TTC (adjudicated)"
    if cramer is not None:
        expected = TTC_UGKGBWDAY[cramer]
        if abs(ttc - expected) > 1e-12:
            _raise("G5-ttc",
                   "adjudicated TTC = %g ug/kg bw/day disagrees with Cramer "
                   "class %s (which would give %g)." % (ttc, cramer, expected),
                   "The adjudicated value is kept. Do not recompute the TTC "
                   "from the Cramer class.", False, log)
    return threshold(ttc * 1e-3, "systemic", name, tier, origin)


# %% Assessment result -------------------------------------------------------

class cosmeticassessment:
    """
    The outcome of one substance against one product.

    Carries the dose, the **binding threshold with its identity and tier**, the
    margin, the provenance of every number that went in, and any gate that
    fired. `verdict` is one of "pass", "exceed" or "not-assessable" -- the
    third is a first-class outcome, never a silent omission.
    """

    __slots__ = ("substance", "product", "route", "dose", "threshold",
                 "verdict", "provenance", "gates", "CF", "note", "asks")

    def __init__(self, substance, product, route, dose, thr, verdict,
                 provenance=None, gates=None, CF=None, note=None, asks=None):
        self.substance = substance
        self.product = product
        self.route = route
        self.dose = dose
        self.threshold = thr
        self.verdict = verdict
        self.provenance = provenance or {}
        self.gates = gates or []
        self.CF = CF                       # kept so a failing screen can seed CF0
        self.note = note
        #: Controlled parameters still at a class default: (name, value, ask).
        #: A question sheet for the customer, not a diagnostic -- see G19.
        self.asks = asks or []

    @property
    def margin(self):
        """Margin of exposure: threshold / dose. Below 1 means exceeded."""
        if self.dose is None or not self.dose or self.threshold is None:
            return None
        return self.threshold.value / self.dose

    @property
    def units(self):
        return "mg/kg bw/day" if self.route == "systemic" else "ug/cm2"

    def __repr__(self):
        name = _identifier(self.substance)
        if self.verdict == "not-assessable":
            out = ["<cosmeticassessment %s in %s -- NOT ASSESSABLE>"
                   % (name, self.product)]
            if self.note:
                out.append("   reason : %s" % self.note)
        else:
            m = self.margin
            out = ["<cosmeticassessment %s in %s>" % (name, self.product),
                   "   route  : %s" % self.route,
                   "   dose   : %.4g %s" % (self.dose, self.units),
                   "   limit  : %.4g %s  [%s] (%s)"
                   % (self.threshold.value, self.threshold.units,
                      self.threshold.identity, self.threshold.tier),
                   "   margin : %s" % ("n/a" if m is None else "%.3g" % m),
                   "   verdict: %s" % self.verdict.upper()]
        if self.provenance:
            out.append("   from   : %s"
                       % ", ".join("%s=%s" % kv
                                   for kv in sorted(self.provenance.items())))
        for g in self.gates:
            out.append("   ! %s" % g.splitlines()[0])
        if self.asks:
            out.append("   ask    : %s (controlled, still defaulted)"
                       % ", ".join(p for p, _, _ in self.asks))
        print("\n".join(out))
        return str(self.verdict)

    __str__ = __repr__


# %% The cosmetic medium -----------------------------------------------------

class cosmetic(foodproperty):
    """
    ===========================================================================
    SFPPy Module: Cosmetic Exposure
    ===========================================================================
    A cosmetic or home-care product as a contacting medium, **plus** the
    exposure stage that converts a concentration into a dose.

    A sibling of `realfood` and `simulant` under `foodproperty`, so it inherits
    the whole transfer machinery -- `h`, `k`, `contacttime`, `surfacearea`,
    `volume`, property propagation, the `%`, `<<` and `>>` operators -- without
    a single change to `food.py`.

    ---------------------------------------------------------------------------
    **Exposure parameters** (SI internally; give any unit as a tuple)
    ---------------------------------------------------------------------------
    - `amountapplied`    : mass of product per application (kg)
    - `frequency`        : applications per unit time (1/s)
    - `retention`        : fraction remaining on skin (1 = leave-on)
    - `exposedarea`      : skin area receiving the product (m2)
    - `bodyweight`       : body weight of the target (kg)
    - `dermalabsorption` : absorbed fraction (CosPaTox Tier 1 default = 1)

    ---------------------------------------------------------------------------
    **The two concentration entries**
    ---------------------------------------------------------------------------
    *Measured* -- set `CF0` to the analysed value, or to half the detection
    limit for a non-detect, and run nothing. `CF0` already exists on
    `foodlayer`, so no new field is needed.

    *Predictive* -- set `C0` on the packaging layer and either take the M0
    mass-balance bound (`totaltransfer`) or run `senspatankar` for M3.

    ---------------------------------------------------------------------------
    **Doctrine**
    ---------------------------------------------------------------------------
    Robust, not accurate. Start at the most conservative tier and refine only
    when a problem is found -- and when one is found, the exit is a branch, not
    an escalation: refine one rung, *or* measure. A tier-M3 result on
    unvalidated parameters is less trustworthy than a tier-M0 bound, because it
    looks authoritative while resting on numbers nobody measured.

    Example
    -------
    ```python
    from patankar.loadpubchem import migrantToxtree
    from patankar.cosmetic import shampoo
    p = migrantToxtree("toluene") % shampoo()
    p.CF0 = 0.5                     # mg/kg, non-detect at DL = 1 mg/kg
    p.assess()
    ```
    """

    level = "property"
    description = "cosmetic or home-care product"
    name = "generic cosmetic product"

    # ---- exposure chain (SI) ----------------------------------------------
    # Values below are neutral placeholders; the concrete categories declare
    # the real ones. They are NOT a usable scenario on their own.
    amountapplied, amountappliedUnits = check_units((1.0, "g"))
    frequency, frequencyUnits = check_units((1.0, "1/day"))
    retention, retentionUnits = check_units((1.0, NoUnits))
    exposedarea, exposedareaUnits = check_units((100.0, "cm**2"))
    bodyweight, bodyweightUnits = check_units((60.0, "kg"))
    dermalabsorption, dermalabsorptionUnits = check_units((1.0, NoUnits))

    # ---- doctrine ----------------------------------------------------------
    _usage = None            # "leave-on" or "rinse-off"
    _target = "adult"        # named human target
    _tier1defaults = True    # 100 % penetration, per CosPaTox Tier 1

    #: Does this product swell the wall it touches? True for every cosmetic --
    #: emulsions, hydroalcoholic bases, surfactant solutions and oxidising
    #: solutions all sorb into a polymer and depress its glass transition.
    #: Set False **only** with evidence for this medium and this polymer; the
    #: declaration is recorded in the provenance either way (G20).
    plasticizing = True

    #: Use `boundingD` -- the envelope of Welle and rubbery Piringer -- as the
    #: diffusivity of every PET layer this product touches. **On by default**,
    #: because the per-substance model switching it replaces destroys the
    #: molecular-size discrimination a functional-barrier study depends on.
    #: Setting it False is recorded in the provenance, never silent.
    useboundingD = True

    # ---- contact conditions ------------------------------------------------
    # `foodphysics` supplies the EU food-contact ACCELERATED TEST convention --
    # 10 days at 40 degC (`food.testcontact`). That is a test protocol for a
    # food, and inheriting it silently is the most likely way for a cosmetic
    # assessment to end up built on the wrong contact conditions: 10 days is
    # more than a factor 50 short of an ordinary cosmetic shelf life.
    #
    # **No equivalent standard exists for cosmetic packaging.** In its absence
    # the values below are a REPRESENTATIVE PLACEHOLDER chosen for cosmetics --
    # ambient storage at 23 degC over an 18-month shelf life -- and deliberately
    # NOT a transposition of the food protocol. They carry no regulatory
    # authority and are not presented as a convention.
    #
    # Both are CONTROLLED parameters (see PARAM_CONTROL), so G19 asks for the
    # real specification whichever default is in place. The purpose of the
    # default is to be honest about the product, not to excuse not asking.
    contacttime, contacttimeUnits = check_units((18.0, "months"))
    contacttemperature, contacttemperatureUnits = check_units((23.0, "degC"))

    # `foodphysics.__init__` copies every public class attribute onto the
    # instance; read-only properties must be exempted here as food.py does for
    # its own (see food.foodphysics._lowLevelPredictionPropertyList).
    _lowLevelPredictionPropertyList = (
        foodproperty._lowLevelPredictionPropertyList
        + ["usage", "target", "dailyproductmass"])

    # Parameters this module owns; the foodphysics registry does not know them,
    # so tuple-to-SI conversion is done here.
    _exposureparams = {
        "amountapplied": "kg",
        "frequency": "1/s",
        "retention": NoUnits,
        "exposedarea": "m**2",
        "bodyweight": "kg",
        "dermalabsorption": NoUnits,
    }

    def __init__(self, **kwargs):
        """Convert exposure parameters to SI, then defer to `foodphysics`."""
        for key in list(kwargs):
            if key in self._exposureparams and isinstance(kwargs[key], tuple):
                kwargs[key], _ = check_units(kwargs[key])
        self._provenance = {}
        self._gatelog = []
        self._noidentity = False
        self._MAE = None            # explicitly set acceptable exposure
        self._CForigin = None       # provenance of an explicitly set CF
        self._CFnote = None
        # A6 -- parameters the caller supplied explicitly. Everything else is a
        # class default, and for a CONTROLLED parameter that is a defect (G19).
        self._explicitparams = set(k for k in kwargs if k in PARAM_CONTROL)
        # A9 -- the wall state a transfer computation actually assumed.
        self._wallstate = None
        # A9/boundingD -- which layers had their D replaced by the envelope.
        self._Dorigin = None
        super().__init__(**kwargs)

    # ---- parameter control (A6) --------------------------------------------

    def __lshift__(self, other):
        """
        `product << packaging3D` -- inherit the real geometry.

        Overridden only to record that surface area and volume then came from a
        real package rather than from a class default, so that G19 stops asking
        for them. The geometry itself is `geometry.py`'s business and is not
        touched here.
        """
        before = tuple(self._scalar(getattr(self, p, None))
                       for p in _GEOMETRY_PARAMS)
        out = super().__lshift__(other)
        for p, was in zip(_GEOMETRY_PARAMS, before):
            if self._scalar(getattr(self, p, None)) != was:
                self._explicitparams.add(p)
        return out

    def update(self, **kwargs):
        """Record which controlled parameters an update supplied, then defer."""
        self._explicitparams.update(k for k in kwargs if k in PARAM_CONTROL)
        return super().update(**kwargs)

    def __rshift__(self, other):
        """
        `product >> packaging` -- propagate, then fix the PET diffusivity.

        This is where `boundingD` becomes the default: the substance has just
        been attached to the layer by `foodphysics.__rshift__`, and the solver
        has not yet run, so it is the one point at which the per-substance model
        switching can be replaced by a single parameterisation before it reaches
        a result. See `applyboundingD`.
        """
        out = super().__rshift__(other)
        try:
            self.applyboundingD(other)
        except Exception:
            pass                       # never let a refinement break a pipeline
        return out

    def _isdefaulted(self, param):
        """
        True when `param` still carries its class default.

        Two signals, and either is sufficient: the parameter was never named in
        a constructor call, an `update()` or a geometry propagation; **and** its
        value is still identical to the class attribute. Requiring both keeps a
        caller who re-supplies the default value on purpose from being accused
        of not having decided.
        """
        if param in self._explicitparams:
            return False
        mine = self._scalar(getattr(self, param, None))
        theirs = self._scalar(getattr(type(self), param, None))
        if mine is None or theirs is None:
            return False
        return mine == theirs

    def controlreport(self):
        """
        Split every parameter by **who controls it**, and flag the defaults.

        Returns a dict with three lists of ``(parameter, value, note)``:

        - ``controlled``   -- the customer decides these and can evidence them.
          They must be **realistic**: the realistic worst case over the declared
          specification envelope, not a bound and not a typical value.
        - ``uncontrolled`` -- nobody measured them for this case. They must be
          **bounded**, and the direction of the bound stated.
        - ``asks``         -- the controlled parameters still sitting at a class
          default. Each carries what to ask the customer for. This list is what
          G19 refuses on, and it is the useful output of the whole method: it is
          a question sheet, not a diagnostic.
        """
        out = {"controlled": [], "uncontrolled": [], "asks": []}
        for param in sorted(PARAM_CONTROL):
            klass, ask = PARAM_CONTROL[param]
            value = self._scalar(getattr(self, param, None))
            defaulted = self._isdefaulted(param)
            note = "defaulted" if defaulted else "supplied"
            out["controlled" if klass == CONTROLLED
                else "uncontrolled"].append((param, value, note))
            if defaulted and klass == CONTROLLED:
                out["asks"].append((param, value, ask))
        return out

    def _gate_control(self, strict=False, log=None):
        """
        G19 -- a controlled parameter left at a class default.

        Not fatal by default: an assessment that has not yet been given the real
        contact conditions is still worth computing, provided nobody mistakes
        the result for one that used them. So this records rather than raises
        unless `strict`, and the ask list travels with the assessment.

        The food accelerated-test convention gets a sentence of its own, because
        inheriting 10 days at 40 degC into a cosmetic is a specific and
        recognisable error rather than a generic missing input.
        """
        report = self.controlreport()
        if not report["asks"]:
            return report
        t = self._scalar(getattr(self, "contacttime", None))
        T = self._scalar(getattr(self, "contacttemperature", None))
        inherited_food_test = (t == _TESTCONTACT_TIME_S
                               and T == _TESTCONTACT_TEMP_C)
        names = ", ".join(p for p, _, _ in report["asks"])
        msg = ("%d controlled parameter(s) still carry a class default: %s. "
               "These are decided and documented by the customer."
               % (len(report["asks"]), names))
        if inherited_food_test:
            msg += (" The contact conditions are the EU food-contact "
                    "ACCELERATED TEST convention (10 days at 40 degC), which "
                    "is a test protocol for a food, not a cosmetic shelf life.")
        remedy = "Ask for: " + "; ".join(
            "%s -- %s" % (p, ask) for p, _, ask in report["asks"])
        _raise("G19-controlclass", msg, remedy, strict, log)
        return report

    # ---- identity ----------------------------------------------------------

    def declare_unidentified(self):
        """
        Declare that **no identity is claimed** -- the T0 case.

        The commonest substance in a recyclate is the one nobody named. Tier T0
        exists for it: the genotoxic TTC applies to any molecule, so one bound
        covers the whole uncharacterised population without naming a single
        structure.

        This is a deliberate declaration, not a fallback. Leaving `substance`
        as ``None`` remains a refusal (G2), so a forgotten injection can never
        become a silent T0 answer.

        `food.foodlayer.substance` type-checks against `migrant`, and rightly
        so: the sentinel is not a substance. The declaration is therefore
        carried here rather than assigned there, and `substance` itself is left
        untouched.

        Returns `self`, so it chains::

            p = bodylotion().declare_unidentified()
            p.CF0 = 1.0                    # at the detection limit
            p.assess(censoring="DL")       # T0, no molecule named
        """
        self._noidentity = True
        return self

    def isunidentified(self):
        """
        True when the uncharacterised population has been declared.

        A method rather than a property: `foodphysics.__init__` walks the MRO
        and `setattr`s every public, non-callable class attribute onto the
        instance, so a read-only property of this name would break
        construction. Callables are skipped, which is why this is one.
        """
        return bool(getattr(self, "_noidentity", False))

    @property
    def _assessedsubstance(self):
        """The substance to assess: the sentinel when no identity is claimed."""
        if self.isunidentified():
            return unidentified
        return self.substance

    # ---- helpers -----------------------------------------------------------

    @staticmethod
    def _scalar(x, default=None):
        if x is None:
            return default
        a = np.asarray(x, dtype=float).ravel()
        return float(a[0]) if a.size else default

    @property
    def usage(self):
        """"leave-on" or "rinse-off"; None on the abstract base."""
        return self._usage

    @property
    def target(self):
        """Named human target the class defaults describe."""
        return self._target

    @property
    def dailyproductmass(self):
        """
        Mass of product reaching the skin per day, kg/day.

        amount applied x frequency x retention. The retention factor is what
        separates a leave-on from a rinse-off product, and it is a distinct
        quantity from any dilution of the migrant into the formulation.
        """
        a = self._scalar(self.amountapplied)
        f = self._scalar(self.frequency)
        r = self._scalar(self.retention)
        return a * f * r * 86400.0          # frequency is stored in 1/s

    # ---- gates -------------------------------------------------------------

    def _gate_substance(self, sub, strict=True, log=None):
        """
        G12 -- substances come from `loadpubchem`, never from elsewhere.

        `migrant` resolves the query, owns the synonym list and the registry
        numbers, carries the molecular weight and the regulatory databases, and
        `migrantToxtree` adds the Cramer classification and the TTC. Accepting a
        duck-typed object would mean reimplementing that identity layer here,
        where it would silently drift out of step with `loadpubchem`.

        The `unidentified` sentinel passes: it is not a duck-typed substance
        but an explicit declaration that no identity is claimed, which is the
        very case tier T0 exists for.
        """
        from patankar.loadpubchem import migrant
        if _is_unidentified(sub):
            return True
        if not isinstance(sub, migrant):
            _raise("G12-substance",
                   "the attached substance is a %s, not a migrant."
                   % type(sub).__name__,
                   "Build it with `migrant(...)` or `migrantToxtree(...)` from "
                   "patankar.loadpubchem, or use cosmetic.as_migrant(query).",
                   strict, log)
            return False
        return True

    def _gate_validity(self, sub, strict=True, log=None):
        """
        G2 -- a substance without a molecular weight has no transfer model.

        `Dpiringer` and the Flory-Huggins partition models both need a
        molecular weight and a structure. UVCBs, polymers and mixtures have
        neither, so they lie outside the validity domain of the transfer model
        whether or not a database holds them. They are reported
        not-assessable, which is a first-class outcome.

        The `unidentified` sentinel passes -- not because it has a molecular
        weight, but because the T0-M0 route never asks for one. The M0 bound is
        a mass balance over geometry and densities, and the dose chain reads
        only product parameters. M3 stays unreachable for it by construction
        (`migrant.__rmod__` rejects the sentinel), so no transfer model is ever
        built without a structure.
        """
        if _is_unidentified(sub):
            return True
        M = getattr(sub, "M", None)
        M = self._scalar(M) if M is not None else None
        if M is None or not np.isfinite(M) or M <= 0:
            _raise("G2-validity",
                   "%s has no usable molecular weight; it is a UVCB, a polymer "
                   "or a mixture." % _identifier(sub),
                   "Report it as not-assessable. Do not estimate a transfer "
                   "for a substance with no single structure.", strict, log)
            return False
        return True

    def _gate_kprovenance(self, strict=False, log=None):
        """
        G3 -- record whether `k` was computed or inherited.

        `food.foodlayer.k0` resolves through `_compute_kmodel()` and falls back
        to the class attribute silently when the substance has no evaluable
        model or `chemicalsubstance` is unset. Injecting the substance with the
        `%` operator is what supplies the real value; the fallback is a
        taxonomy placeholder and must not be reported as a computed partition.
        """
        computed = None
        if hasattr(self, "_compute_kmodel"):
            try:
                computed = self._compute_kmodel()
            except Exception:
                computed = None
        if computed:
            return ESTIMATED
        _raise("G3-kprovenance",
               "k = %s is the affinity-mixin class default, not a computed "
               "partition coefficient." % self._scalar(getattr(self, "k", None)),
               "Inject the substance with `substance %% product` so the "
               "Flory-Huggins model can evaluate, or supply k explicitly.",
               strict, log)
        return DEFAULTED

    def _gate_route(self, route, sub, strict=True, log=None):
        """G7 -- refuse a local assessment when no local threshold exists."""
        if route == "systemic":
            return True
        local = getattr(sub, "DNELlocal", None)
        if local is None:
            _raise("G7-route",
                   "no site-of-contact threshold is available for %s; a "
                   "systemic value cannot stand in for it."
                   % _identifier(sub),
                   "Report the local route as not-assessable. Systemic "
                   "[mg/kg bw/day] and local [ug/cm2] are not commensurable.",
                   strict, log)
            return False
        return True

    # ---- the exposure chain ------------------------------------------------

    def systemicdose(self, CF=None, CFunit="mg/kg", strict=True, log=None):
        """
        Systemic dose in mg/kg bw/day from a concentration in the product.

            E = CF * (amount x frequency x retention) * absorption / bodyweight

        `CF` defaults to `CF0`, which is where the *measured* entry lives.
        """
        log = self._gatelog if log is None else log
        CF = self._scalar(self.CF0) if CF is None else CF
        c = _convert(CF, CFunit, _CONC_UNITS, "G1-units",
                     "product concentration", strict, log)
        if c is None:
            return None
        alpha = self._scalar(self.dermalabsorption, 1.0)
        bw = self._scalar(self.bodyweight)
        if not bw or bw <= 0:
            _raise("G1-units", "body weight is not set.",
                   "Give bodyweight=(value,'kg').", True, log)
        # mg/kg x kg/day / kg = mg/kg bw/day
        return c * self.dailyproductmass * alpha / bw

    def localdose(self, CF=None, CFunit="mg/kg", strict=True, log=None):
        """
        Local (site-of-contact) load in ug/cm2 **per application**.

            L = CF * amount applied * retention / exposed area

        Deliberately per application rather than per day: a surface load is not
        additive over time in the way a systemic dose is.
        """
        log = self._gatelog if log is None else log
        CF = self._scalar(self.CF0) if CF is None else CF
        c = _convert(CF, CFunit, _CONC_UNITS, "G1-units",
                     "product concentration", strict, log)
        if c is None:
            return None
        a = self._scalar(self.amountapplied)
        r = self._scalar(self.retention)
        area = self._scalar(self.exposedarea)
        if not area or area <= 0:
            _raise("G1-units", "exposed area is not set.",
                   "Give exposedarea=(value,'cm**2').", True, log)
        # (mg/kg) x kg / m2 = mg/m2 -> ug/cm2 is x 1e3 / 1e4 = x 0.1
        return c * a * r / area * 0.1

    def maxconcentration(self, thr=None, strict=True, log=None):
        """
        The product's own SML-equivalent, mg/kg.

        Inverting the systemic chain at a threshold yields the number EU
        practice hands the food case by convention -- here computed per product
        and per target rather than assumed:

            Cmax = threshold x bodyweight / (amount x frequency x retention x absorption)
        """
        log = self._gatelog if log is None else log
        if thr is None:
            thr = self.binding_threshold(strict=strict, log=log)
        if thr is None:
            return None
        if thr.route != "systemic":
            _raise("G7-route",
                   "maxconcentration inverts the systemic chain; a local "
                   "threshold cannot be used.",
                   "Use the local route explicitly.", strict, log)
            return None
        alpha = self._scalar(self.dermalabsorption, 1.0)
        bw = self._scalar(self.bodyweight)
        denom = self.dailyproductmass * alpha
        if denom <= 0:
            return None
        return thr.value * bw / denom

    def maxdetectionlimit(self, censoring="DL/2", strict=True, log=None):
        """
        The detection limit the analytical method must achieve, mg/kg.

        For an unidentified substance the useful question is not whether it
        passes but what the laboratory must reach. A non-detect is counted at
        half the detection limit by default -- the standard treatment of
        left-censored exposure data -- so the tolerable DL is twice the maximum
        tolerable concentration.

        G4: for the **T0 bound** on an uncharacterised substance, pass
        ``censoring="DL"``. Halving a bound is not the same operation as taking
        the midpoint of a censored distribution, and the CosPaTox dossier
        counts non-detects at DL with no moderating factor.
        """
        log = self._gatelog if log is None else log
        if censoring not in ("DL/2", "DL"):
            _raise("G4-censoring", "unknown censoring convention %r." % censoring,
                   "Use 'DL/2' (identified substance) or 'DL' (T0 bound).",
                   True, log)
        cmax = self.maxconcentration(strict=strict, log=log)
        if cmax is None:
            return None
        return cmax * (2.0 if censoring == "DL/2" else 1.0)

    # ---- transfer ----------------------------------------------------------

    def totaltransfer(self, packaging, C0=None, C0unit="mg/kg",
                      strict=True, log=None):
        """
        M0 -- total transfer bound, mg/kg of product. No solver.

            CF = C0 * (rho_P * A * l) / (rho_F * V_F)

        This is a mass balance over the real volumes, and it reproduces the
        CosPaTox dilution estimator exactly: substituting L = rho_F V_F /
        (rho_P A l) returns CF = C0 / L. The dilution factor is therefore not a
        parameter of this module -- it is a consequence of geometry, wall and
        densities, all of which SFPPy already carries.

        G6: `C0` is a **mass** concentration in the polymer and `CF` a mass
        concentration in the product. The solver works on volumes and never
        converts, so the two densities enter explicitly here.

        G9: refuses a polymer used outside its state -- the `r` prefix in SFPPy
        means *rubbery*, not recycled.

        G20: refuses the GLASSY state for a product that swells the wall. For a
        cosmetic the swollen state is the conservative default; the glassy one
        is the claim that must be argued.
        """
        log = self._gatelog if log is None else log
        if not isinstance(packaging, layer):
            raise TypeError("packaging must be a patankar.layer, not %s"
                            % type(packaging).__name__)
        self._gate_plasticization(packaging, strict=strict, log=log)
        self._gate_polymerstate(packaging, strict=strict, log=log)

        c0 = self._scalar(packaging.C0) if C0 is None else C0
        c0 = _convert(c0, C0unit, _CONC_UNITS, "G1-units",
                      "initial concentration in the polymer", strict, log)
        if c0 is None:
            return None

        A = self._scalar(self.surfacearea)
        VF = self._scalar(self.volume)
        rhoF = self._scalar(self.density)
        l = self._scalar(packaging.l)
        # G6 -- the polymer density enters the mass balance directly, so a
        # silent default is a silent error in the answer. Falling back to
        # 1000 kg/m3 for PET (1350) understates the wall mass by 26 %, i.e.
        # NON-conservatively, in a module that gates the spelling of a unit.
        # The layer is asked; if it cannot say, the gate fires.
        rhoP = None
        try:                                   # real polymers override density(T)
            d = packaging.density()
            rhoP = self._scalar(d[0] if isinstance(d, tuple) else d, None)
        except Exception:
            pass
        if rhoP is None:
            # Base `layer` exposes no density() but carries rho = 1000 kg/m3 as
            # a CLASS default. Accepting it here would silently reinstate the
            # very number this gate exists to refuse, so the instance value is
            # only used when it differs from the class default, i.e. when
            # somebody actually declared it.
            declared = self._scalar(getattr(packaging, "rho", None), None)
            if declared is not None and declared != _LAYER_RHO_DEFAULT:
                rhoP = declared
        if not rhoP:
            _raise("G6-massbasis",
                   "%s does not expose a density, and the mass balance will "
                   "not assume one: a default of 1000 kg/m3 understates a PET "
                   "wall by 26 %%, non-conservatively."
                   % type(packaging).__name__,
                   "Give the layer a density -- `rho=(1350,'kg/m**3')` -- or "
                   "use a polymer class that defines `density()`.", True, log)
        if not (A and VF and rhoF and l):
            _raise("G6-massbasis",
                   "the mass balance needs surface area, product volume, the "
                   "product density and the wall thickness; one is missing.",
                   "Propagate the geometry with `product << packaging3D` and "
                   "give the layer a thickness.", True, log)
        mP = rhoP * A * l                      # mass of packaging wall, kg
        mF = rhoF * VF                         # mass of product, kg
        return c0 * mP / mF

    # ---- conservative coupling (A8) ----------------------------------------

    def volumeratio(self, packaging):
        """
        Polymer-to-product volume ratio $V_P/V_F$ -- dimensionless.

        The number that decides whether a partition shortcut is safe. Note that
        it is **not** the surface-to-volume ratio: $A/V$ decides which format is
        worst for the dose, $V_P/V_F$ decides whether unitary partition may
        stand in for total transfer (G18).
        """
        A = self._scalar(self.surfacearea)
        VF = self._scalar(self.volume)
        l = self._scalar(getattr(packaging, "l", None))
        if not (A and VF and l):
            return None
        return (A * l) / VF

    def unitarypartition(self, packaging, strict=False, log=None):
        """
        Remove thermodynamic retention: set $K_{P/F} = 1$.

        **What it is for.** A charged or strongly hydrophilic species has no
        evaluable Flory-Huggins partition -- the model is parameterised for
        neutral organics -- so rather than invent a coefficient, maximise the
        affinity for the medium and let the wall retain nothing. Combined with
        a plasticized layer (G20) and the bare species' molar mass in Piringer,
        this bounds a migrant whose speciation is unknown **without modelling
        the speciation at all**: the substance is already in the polymer, and
        the question is how it leaves, not how it got there.

        **Why it is not the worst case.** $K_{P/F} = 1$ says the substance has
        no preference. A genuinely hydrophilic ion has $K_{P/F} \\ll 1$ and
        releases more. What that costs is decided entirely by the volume ratio,
        through the equilibrium release fraction

            PR_E = 1 / (1 + K_PF * V_P/V_F)

        so the conservatism given away against the total-transfer bound is
        ``1 - PR_E``. For a bottle it is a few per cent and unitary partition is
        operationally M0; for a thick-walled stick or compact it is a third or
        more, and M0 must be used instead. G18 draws that line rather than
        leaving it to judgement.

        Returns a dict carrying the ratio, the release fraction, the loss and
        the tier it recommends.
        """
        log = self._gatelog if log is None else log
        r = self.volumeratio(packaging)
        if r is None:
            _raise("G6-massbasis",
                   "the volume ratio needs surface area, product volume and a "
                   "wall thickness; one is missing.",
                   "Propagate the geometry with `product << packaging3D` and "
                   "give the layer a thickness.", True, log)
            return None
        PR = 1.0 / (1.0 + r)                    # K_PF = 1
        loss = 1.0 - PR
        self.k = 1.0
        try:
            packaging.k = 1.0
        except Exception:                        # layer may refuse; K is a ratio
            pass                                 # and the medium side suffices
        tier = "M1"
        if loss > UNITARY_PARTITION_MAXLOSS:
            tier = "M0"
            _raise("G18-partitionbound",
                   "V_P/V_F = %.3g, so unitary partition releases only %.1f %% "
                   "of the wall inventory and gives away %.1f %% of the "
                   "total-transfer bound." % (r, 100 * PR, 100 * loss),
                   "Use totaltransfer() -- M0 is unconditional. Unitary "
                   "partition is a convenience for thin-walled, high-volume "
                   "formats (V_P/V_F <~ %.2g), not a general substitute."
                   % UNITARY_PARTITION_MAXLOSS, strict, log)
        return {"volumeratio": r, "releasefraction": PR, "lossvsM0": loss,
                "tier": tier}

    def completereaction(self, C0parent, products=None, consume_parent=False,
                         C0unit="mg/kg", strict=True, log=None):
        """
        Bound a reactive contact by assuming the reaction has **gone to
        completion** -- without any kinetics.

        **The move.** `migration.py` integrates a Fickian balance with no
        volumetric source or sink term, so a medium that attacks the wall
        cannot be represented as chemistry. It can be represented as
        *inventory*: assume every convertible bond has already converted, so the
        degradation products stand at their stoichiometric ceiling at $t = 0$.
        A source-term problem becomes an initial-concentration problem, which
        M0 already handles exactly, and `radigen` is needed only if the bound
        that results fails.

        **The direction rule (G17).** Completion is conservative for what it
        *forms* and anti-conservative for what it *consumes*. If the medium
        oxidises the very substance under assessment, assuming completion makes
        that substance vanish and returns an unbounded margin -- an assessment
        that passes by destroying its own analyte. So completion is applied in
        the exposure-raising direction only: full formation of products, never
        full consumption of a parent.

        **What it does not bound.** A stoichiometric ceiling for a product
        nobody listed bounds nothing. The residual population is exactly the
        uncharacterised one, so a reactive case does not close here -- it
        reduces to the T0 case, and `declare_unidentified()` answers it. That is
        stated in the returned note rather than left to be discovered.

        Parameters
        ----------
        C0parent : float
            Initial concentration of the parent substance in the polymer.
        products : dict, optional
            ``{name: ceiling}`` -- stoichiometric maximum of each degradation
            product, in the same unit. Enumerated by the resin chemistry, not
            by this module.
        consume_parent : bool
            Request depletion of the parent by the reaction. Refused (G17).
        """
        log = self._gatelog if log is None else log
        if consume_parent:
            _raise("G17-reactiondirection",
                   "complete conversion was requested in the direction that "
                   "CONSUMES the substance under assessment.",
                   "Apply completion only where it raises exposure: full "
                   "formation of products, never full consumption of a parent. "
                   "If depletion of the parent is real and material, it is a "
                   "measurement, not an assumption.", strict, log)
            return None
        c0 = _convert(C0parent, C0unit, _CONC_UNITS, "G1-units",
                      "parent concentration in the polymer", strict, log)
        if c0 is None:
            return None
        listed = {}
        for name, value in (products or {}).items():
            v = _convert(value, C0unit, _CONC_UNITS, "G1-units",
                         "ceiling for %s" % name, strict, log)
            if v is not None:
                listed[name] = v
        total = c0 + sum(listed.values())
        return {
            "parent": c0,
            "products": listed,
            "inventory": total,
            "unit": "mg/kg",
            "note": ("complete-conversion inventory: parent retained in full "
                     "(G17) plus the stoichiometric ceiling of %d enumerated "
                     "product(s). Products NOT enumerated are not bounded by "
                     "this number -- assess them as the uncharacterised "
                     "population with declare_unidentified()." % len(listed)),
        }

    def _gate_polymerstate(self, packaging, strict=True, log=None):
        """
        G9 -- the `r` prefix means RUBBERY, not recycled.

        `rPET` is PET above its glass transition (~76 degC); a recycled PET
        bottle at ambient temperature is the GLASSY polymer carrying a
        different `C0` -- or, in contact with a cosmetic, the PLASTICIZED one
        (see G20). Using `rPET` below Tg overestimates migration by orders of
        magnitude while looking conservative rather than wrong.

        Three states, three different things, and only the first two are about
        temperature:

        =========  ==============================================  ===========
        prefix     meaning                                         set by
        =========  ==============================================  ===========
        ``g``      glassy, T < Tg                                  temperature
        ``r``      RUBBERY, T > Tg                                 temperature
        ``w``      PLASTICIZED (swollen), Tg depressed by sorption  the medium
        =========  ==============================================  ===========

        A recyclate is none of them: it is the same polymer with a different
        `C0`. `w` is *not* "wet" in the sense of "water only" -- it is the
        general swollen state, whatever does the swelling (G20).
        """
        cls = type(packaging).__name__
        Tg = getattr(packaging, "Tg", None)
        if Tg is None:
            return True
        Tgv = Tg[0] if isinstance(Tg, tuple) else Tg
        T = self._scalar(getattr(packaging, "T", None))
        if T is None:
            return True
        if cls.startswith("r") and T < Tgv:
            # A swelling medium changes what this means. Below Tg the rubbery
            # parameterisation is not a state error but the recommended
            # PLASTICIZATION SURROGATE: no Tg-aware model exists for a swollen
            # wall, so the rubbery parameters stand in for one. Refuse it only
            # where the contact was declared non-swelling, i.e. where the r-
            # prefix really was mistaken for "recycled".
            if self.plasticizing:
                return True
            sib = PLASTICIZED_SIBLING.get(GLASSY_OF.get(cls, ""), None)
            _raise("G9-polymerstate",
                   "%s is the RUBBERY state (Tg = %g degC) but the contact "
                   "temperature is %g degC, and this contact was declared "
                   "NON-swelling. The r- prefix does not mean recycled."
                   % (cls, Tgv, T),
                   "Use the glassy class for a genuinely dry contact"
                   + (", or %s if the wall does swell after all" % sib
                      if sib else "")
                   + ". A recyclate is the same polymer with a different C0, "
                     "not a different state.",
                   strict, log)
            return False
        return True

    @staticmethod
    def _Dof(classname, sub, l=None, T=None, Tg=None):
        """
        Evaluate D for `sub` in a fresh layer of `classname`. None if it cannot.

        A fresh instance, deliberately: the caller's layer must not acquire a
        substance or a Tg as a side effect of being checked.
        """
        import patankar.layer as _lay
        cls = getattr(_lay, classname, None)
        if cls is None:
            return None
        kw = {}
        if l is not None:
            kw["l"] = l
        if T is not None:
            kw["T"] = (T, "degC")
        if Tg is not None:
            kw["Tg"] = (Tg, "degC")
        try:
            probe = cls(**kw)
            probe.substance = sub
            return float(np.asarray(probe.D).ravel()[0])
        except Exception:
            return None

    def _stateD(self, packaging, sub):
        """
        D in the glassy and the plasticized state, for THIS substance.

        Returned as ``(glassyname, Dglassy, plasticname, Dplastic)`` with any
        member None when it could not be evaluated.

        **Why this is measured rather than assumed.** SFPPy does not evaluate
        both states with the same diffusion model. The hole-free-volume model
        (`DFV`) is selected only for toluene; `Dwelle` covers `gPET` but *not*
        `wPET`; everything else falls back to Piringer, where `wPET` carries the
        rubbery-PET parameters. So a glassy-to-plasticized comparison can cross
        a model boundary, and across a model boundary the ordering is not
        guaranteed by physics. It has to be checked.
        """
        cls = type(packaging).__name__
        glassy = cls if cls in PLASTICIZED_SIBLING else None
        plastic = PLASTICIZED_SIBLING.get(cls, None)
        if plastic is None:                      # already plasticized?
            for g, w in PLASTICIZED_SIBLING.items():
                if cls == w:
                    glassy, plastic = g, w
                    break
        if glassy is None or plastic is None:
            return (None, None, None, None)
        l = getattr(packaging, "l", None)
        T = self._scalar(getattr(packaging, "T", None))
        return (glassy, self._Dof(glassy, sub, l, T),
                plastic, self._Dof(plastic, sub, l, T))

    def boundingD(self, packaging, sub=None):
        """
        The bounding diffusivity for a **swollen** wall: the envelope of the
        models available, not any one of them.

        **Why an envelope and not a choice.** No diffusion model in SFPPy is
        validated for a plasticized polymer, and the two candidates fail in
        opposite directions:

        - **Welle** (`Dwelle`) is PET-specific and resolves molecular size
          through the van der Waals volume rather than through $M$ alone, which
          makes it the better *structural* model. But its values come from **dry
          PET**, it carries **no Tg dependence at all**, and it was established
          by fitting lag times in GC-MS -- a methodology whose inversion is
          contested. Its purpose in the literature was to establish PET as a
          high barrier, which is the opposite of the question a safety
          assessment asks. Used as-is it collapses steeply with molecular size:
          at 23 degC it puts Irganox 1076 at ~1e-26 m2/s, i.e. "nothing ever
          migrates".
        - **Piringer with the RUBBERY parameters** carries no molecular-shape
          information -- the glassy-to-rubbery step is a constant offset in
          $A_P$ (3.1 -> 6.4), so it multiplies D by e^3.3 ~ 27 for every
          substance alike -- but it is the only available stand-in for a wall
          whose mobility has been raised by sorption. **Using the rubbery
          parameterisation below Tg is deliberate**, and it is the recommended
          practice here: plasticization of PET in contact with real formulations
          has not been studied by the people who set safety limits, so the
          rubbery state is the defensible surrogate rather than an error.

        Neither bounds the other everywhere. Measured at 23 degC, rubbery
        Piringer exceeds Welle by 2 to 7 orders of magnitude for anything above
        about M = 45, and **falls below it for the smallest migrants** --
        formaldehyde comes out ~15x faster under Welle. So the bound is the
        larger of the two, taken per substance.

        Returns ``(D, modelname)``, with D None when neither model evaluates.
        """
        sub = self._assessedsubstance if sub is None else sub
        if sub is None or _is_unidentified(sub):
            return (None, None)
        T = self._scalar(getattr(packaging, "T", None))
        if T is None:
            return (None, None)
        candidates = []
        # Welle, structural but dry
        try:
            from patankar.property import Dwelle
            V = self._scalar(getattr(sub, "volumeDwelle", None))
            if V:
                candidates.append((float(np.asarray(
                    Dwelle.evaluate("gPET", V, T)).ravel()[0]),
                    "Dwelle (dry PET)"))
        except Exception:
            pass
        # Piringer with the rubbery parameters, the plasticization surrogate
        try:
            from patankar.property import Dpiringer
            M = self._scalar(getattr(sub, "M", None))
            if M:
                candidates.append((float(np.asarray(
                    Dpiringer.evaluate("rPET", M, T)).ravel()[0]),
                    "Piringer rubbery (plasticization surrogate)"))
        except Exception:
            pass
        if not candidates:
            return (None, None)
        return max(candidates, key=lambda c: c[0])

    def applyboundingD(self, packaging, force=False):
        """
        Install `boundingD` as the diffusivity of every PET layer. **Default.**

        **The regression this repairs.** SFPPy selects a different diffusion
        model per (substance, polymer) pair: `DFV` for toluene only, `Dwelle`
        for `gPET`, Piringer elsewhere. In `wPET` at 40 degC that put toluene at
        4.94e-16 and limonene at 5.24e-16 -- within 6 %, and **the wrong way
        round**, the larger molecule coming out faster. A functional-barrier
        study lives entirely on the size discrimination those two numbers were
        supposed to carry, so mixed models do not merely add scatter: they
        delete the effect being studied.

        Applying one rule to every substance restores it -- toluene 1.04e-15
        against limonene 5.24e-16, a factor 2 in the right direction, both from
        the same parameterisation.

        Only PET layers are touched; every other layer keeps the diffusivity its
        own model gives it. Set `useboundingD = False` on the product to opt out,
        which is recorded rather than silent.
        """
        if not (self.useboundingD or force):
            return None
        touched = self._boundingDvector(packaging, report=True)
        if not touched:
            return None
        # A closure, not a frozen array: the substance attached to the wall can
        # change after this point, and a D that silently kept the old migrant's
        # value would be the same class of defect this replaces.
        packaging.Dmodel = lambda: self._boundingDvector(packaging)
        self._Dorigin = {"applied": touched,
                         "note": "one parameterisation for every substance; "
                                 "see cosmetic.boundingD"}
        return touched

    def _boundingDvector(self, packaging, report=False):
        """
        Per-layer diffusivity with PET layers replaced by `boundingD`.

        Returns the vector, or -- with ``report=True`` -- the list of
        ``(index, class, D, model)`` actually substituted. `None` when nothing
        applies, so the caller can leave the layer's own model alone.

        The migrant is read from **the layer**, falling back to the product.
        A diffusivity is a property of the substance in that layer, so binding
        it to the product's migrant would let a wall reused with a different
        substance keep the previous one's D -- silently, and looking computed.
        """
        sub = getattr(packaging, "substance", None) or self._assessedsubstance
        if sub is None or _is_unidentified(sub):
            return None
        try:
            classes = list(packaging.layerclass_history)
            base = np.asarray(packaging._compute_Dmodel(), dtype=float).ravel()
            Tv = np.asarray(packaging.T, dtype=float).ravel()
        except Exception:
            return None
        if len(base) != len(classes):
            return None
        out = base.copy()
        touched = []
        for i, cname in enumerate(classes):
            if cname not in PET_FAMILY:
                continue
            D, model = self.boundingD(
                _TProbe(Tv[i] if i < len(Tv) else Tv[0]), sub)
            if D is None:
                continue
            out[i] = D
            touched.append((i, cname, D, model))
        if not touched:
            return None
        return touched if report else out

    def _gate_tgnoop(self, packaging, strict=False, log=None):
        """
        G21 -- a declared Tg the diffusion model never reads.

        `Tg` is declarable on the plasticized classes so that a medium swelling
        more than water can push the bound further. That only works when the
        selected D model consumes Tg -- which, in SFPPy today, means the
        hole-free-volume model, and that model is selected **for toluene only**
        (`loadpubchem.Dmodel_extensions["DFV"]`). Every other substance falls
        back to Piringer, which is parameterised by molar mass and ignores Tg
        entirely.

        So a declared Tg is silently inert for almost every migrant. Recording
        an assumption the arithmetic ignores is worse than not declaring it, and
        this gate says so rather than letting the result look bounded.
        """
        if getattr(packaging, "_Tguser", None) is None:
            return True
        sub = self._assessedsubstance
        if sub is None or _is_unidentified(sub):
            return True
        cls = type(packaging).__name__
        l = getattr(packaging, "l", None)
        T = self._scalar(getattr(packaging, "T", None))
        Tgv = self._scalar(packaging._Tguser)
        D1 = self._Dof(cls, sub, l, T, Tg=Tgv)
        D2 = self._Dof(cls, sub, l, T, Tg=Tgv - 10.0)
        if D1 is None or D2 is None or D1 != D2:
            return True
        _raise("G21-tgnoop",
               "Tg = %g degC was declared on %s, but D is unchanged when Tg "
               "moves 10 K: the diffusion model selected for %s does not read "
               "Tg." % (Tgv, cls, _identifier(sub)),
               "The Tg lever acts through the hole-free-volume model, which is "
               "selected for toluene only. For this substance, bound the "
               "transport by supplying D explicitly, or stay with the "
               "total-transfer bound, which needs no D at all.", strict, log)
        return False

    def _gate_plasticization(self, packaging, strict=True, log=None):
        """
        G20 -- for a cosmetic, the glassy state is the OPTIMISTIC one.

        Food-contact practice reaches for the glassy class because a dry food
        does not swell the wall. A cosmetic does: an emulsion, a hydroalcoholic
        base, a surfactant solution and an oxidising solution all sorb into the
        polymer and depress its glass transition. Sorption is fast relative to
        the shelf life, so **assuming plasticization is complete costs nothing
        and needs no kinetics** -- which is why it is the default here and the
        glassy state is what must be argued for.

        The size of what is at stake is not marginal. Through the hole-free
        volume model (`property.py`, Zhu, Welle & Vitrac, *Soft Matter*, 2019),
        the plasticized and glassy parameterisations of PET differ by about two
        orders of magnitude in D at 40 degC.

        Declare `plasticizing = False` on the product **only** with evidence
        that the medium does not swell this polymer. The declaration is
        recorded in the provenance either way.
        """
        # A9 -- record the wall state actually assumed, so that every result can
        # state it. A plasticization bound whose Tg is not reported is not a
        # stated assumption, it is a hidden one.
        cls = type(packaging).__name__
        Tg = getattr(packaging, "Tg", None)
        Tgv = (Tg[0] if isinstance(Tg, tuple) else Tg) if Tg is not None else None
        self._wallstate = {
            "class": cls,
            "Tg_degC": None if Tgv is None else float(Tgv),
            "plasticized": cls in PLASTICIZED_SIBLING.values(),
            "medium_swells": bool(self.plasticizing),
        }
        self._gate_tgnoop(packaging, strict=strict, log=log)
        if not self.plasticizing:
            return True
        sib = PLASTICIZED_SIBLING.get(cls, None)
        if sib is None:
            return True

        # Verify the claim rather than assert it. The two states are not always
        # evaluated by the same diffusion model, and across a model boundary the
        # ordering is a property of the models, not of the physics: for a very
        # small migrant the glassy state can come out FASTER, in which case
        # steering to the swollen class would be anti-conservative -- the exact
        # error this gate exists to prevent.
        sub = self._assessedsubstance
        Dg = Dw = None
        if sub is not None and not _is_unidentified(sub):
            _, Dg, _, Dw = self._stateD(packaging, sub)
        if Dg is not None and Dw is not None and Dw <= Dg:
            _raise("G20-plasticization",
                   "%s (glassy) is modelled FASTER than %s (plasticized) for "
                   "%s -- D = %.3g vs %.3g m2/s -- so the swollen state is not "
                   "the conservative one here. The two states are evaluated by "
                   "different diffusion models, and the ordering crossed over."
                   % (cls, sib, _identifier(sub), Dg, Dw),
                   "Keep %s, which is the bounding state for this substance, "
                   "or avoid the question entirely with totaltransfer(): M0 "
                   "needs no D and is unconditional." % cls,
                   strict, log)
            return False

        verified = "" if Dw is None or Dg is None else (
            " (checked: D = %.3g vs %.3g m2/s)" % (Dw, Dg))
        _raise("G20-plasticization",
               "%s is the GLASSY state, but %s swells the wall it touches; the "
               "glassy state is the optimistic one, not the safe one.%s"
               % (cls, self.name, verified),
               "Use %s -- the same polymer in its swollen state -- or declare "
               "`plasticizing = False` on the product with evidence that this "
               "medium does not swell %s." % (sib, cls),
               strict, log)
        return False

    # ---- assessment --------------------------------------------------------

    # ---- explicit entries: feeding a decision without computing it ---------

    def setCF(self, value, unit="mg/kg", origin=MEASURED, note=None):
        """
        Set the concentration in the product **without computing it**.

        Most decisions are not made from a simulation. They are made from a
        detection limit, an analytical result, a supplier declaration, or a
        concentration back-calculated from an acceptable exposure. All of those
        are legitimate inputs; what matters is that the result says which one
        it was.

        `origin` is free text, but a short vocabulary keeps results comparable:

        ==================  ====================================================
        ``measured``        an analytical determination on the product
        ``censored``        at or below a detection limit (see `assess`'s
                            ``censoring=`` for the DL/2 vs DL convention, G4)
        ``MAE-derived``     back-calculated from an acceptable exposure, i.e.
                            the output of `maxconcentration`
        ``supplier``        declared by the material or formulation supplier
        ``regulatory``      fixed by a limit rather than observed
        ``assumed``         a working hypothesis; the weakest tag
        ==================  ====================================================

        Overrides any earlier value and takes precedence over `CF0`. A
        simulation, if one is run afterwards, still wins -- the pipeline is the
        more specific statement.
        """
        c = _convert(value, unit, _CONC_UNITS, "G1-units",
                     "product concentration", True, self._gatelog)
        self.CF0 = c
        self._CForigin = str(origin)
        self._CFnote = note
        return self

    def setMAE(self, value, unit="mg/kg bw/day", route="systemic",
               origin=ADJUDICATED, tier="T2", identity=None, note=None):
        """
        Set the **Maximum Acceptable Exposure** explicitly.

        The threshold ladder (T2 SML-derived TDI, T1 Cramer TTC, T0 genotoxic
        TTC) covers what can be *derived* from the substance. It does not cover
        what an organisation *decides*: an internal limit, a customer
        specification, a DNEL taken from a registration dossier, a value agreed
        in a dossier review. Those are equally binding and often tighter.

        Setting an MAE here makes it the binding threshold, ahead of the whole
        ladder, and every result reports it as such with the origin given.

        The unit is the one that makes exposure and threshold comparable --
        ``mg/kg bw/day`` on the systemic route. On the local route it is
        ``ug/cm2``; the two are never coerced into one another (G7).
        """
        table = _SYSTEMIC_UNITS if route == "systemic" else _LOCAL_UNITS
        v = _convert(value, unit, table, "G1-units",
                     "maximum acceptable exposure", True, self._gatelog)
        self._MAE = threshold(v, route,
                              identity or "MAE (set explicitly)",
                              tier, str(origin), note=note)
        return self

    def clearMAE(self):
        """Drop an explicitly set MAE and fall back to the threshold ladder."""
        self._MAE = None
        return self

    # ---- the reverse question: what may the MATERIAL contain? --------------

    def potentialrelease(self, packaging, CP0guess=None, tier="M0",
                         strict=True, log=None):
        """
        Potential release ``PR = CF / CP0``, and the concentration that produced it.

        `PR` is the quantity that makes the reverse problem tractable, because
        **it does not depend on** ``CP0``. The transport problem is linear in
        the initial concentration: double the load in the wall and every
        concentration in the product doubles. So one simulation at any
        convenient guess characterises the system, and the answer scales.

        Returned as ``(PR, CP0guess, CFguess, provenance)``.

        Two tiers, both of them simulations in the sense that matters -- a
        model is run and a number comes out:

        ``M0``
            the total-transfer mass balance of `totaltransfer`. Analytic, no
            solver, and exactly linear by construction.
        ``M3``
            the Patankar solution. Requires the pipeline to have been run
            (``substance % product << geometry >> packaging >> product``); the
            stored ``lastsimulation`` is used.

        This is SFPPy's own potential-release quantity expressed on a
        concentration basis. `migration.SensPatankarResult.PR` carries the
        mass-basis family (``PRE = mFeq/m0``, ``PR(CF) = CF*VF/m0``,
        ``PRT = PR/PRE``); both are linear in ``CP0`` and the inversion below
        is identical either way.
        """
        log = self._gatelog if log is None else log
        if tier not in ("M0", "M3"):
            _raise("G8-stoponpass", "unknown transfer tier %r." % tier,
                   "Use 'M0' (mass balance) or 'M3' (solver).", True, log)

        if tier == "M0":
            cp0 = self._scalar(packaging.C0) if CP0guess is None else CP0guess
            if not cp0 or cp0 <= 0:
                _raise("G6-massbasis",
                       "the guess concentration in the wall is zero or unset; "
                       "PR cannot be characterised from it.",
                       "Give the layer a C0, or pass CP0guess explicitly.",
                       True, log)
            cf = self.totaltransfer(packaging, C0=cp0, strict=strict, log=log)
            prov = "M0 (mass balance)"
        else:
            if not (self.hassimulation and self.lastsimulation is not None):
                _raise("G8-stoponpass",
                       "tier M3 needs a solved simulation and none is stored.",
                       "Run the pipeline first: "
                       "`substance %% product << geometry >> packaging >> product`, "
                       "or ask for tier='M0'.", True, log)
            cp0 = self._scalar(packaging.C0) if CP0guess is None else CP0guess
            if not cp0 or cp0 <= 0:
                _raise("G6-massbasis",
                       "the wall carries no initial concentration; PR is "
                       "undefined.", "Set C0 on the layer.", True, log)
            cf = float(np.asarray(self.lastsimulation.CF).ravel()[-1])
            prov = "M3 (Patankar solver)"

        cf = None if cf is None else float(np.asarray(cf).ravel()[-1])
        if cf is None or not np.isfinite(cf):
            _raise("G6-massbasis",
                   "the transfer calculation returned no usable concentration.",
                   "Check the geometry, the wall thickness and both densities.",
                   True, log)
        return cf / cp0, cp0, cf, prov

    def maxinitialconcentration(self, packaging, MAE=None, route="systemic",
                                CP0guess=None, tier="M0", verify=True,
                                rtol=1e-6, strict=True, log=None):
        """
        The reverse question: **what may the MATERIAL contain?**

        Compliance asks whether a product is safe. Design asks what a supplier
        may ship. For a recyclate the second is the useful one, because it is a
        specification a decontamination process can be held to.

        Three steps, and the middle one is where the linearity earns its keep::

            1.  MAE  ->  CFmax     invert the exposure chain (maxconcentration)
            2.  CP0guess -> CFguess    one run at any convenient guess
            3.  CP0max = CP0guess * CFmax / CFguess = CFmax / PR

        Step 3 is exact because ``PR = CF/CP0`` does not depend on ``CP0``. The
        guess is genuinely arbitrary: it cancels.

        `MAE` may be given directly (mg/kg bw/day), or left to `binding_threshold`
        -- which itself honours an MAE previously set with `setMAE`.

        Robustness
        ----------
        `verify=True` re-runs the transfer at ten times the guess and checks
        that ``CF`` scales by ten. That is the assumption the whole inversion
        rests on, so it is tested rather than trusted. A concentration-dependent
        diffusivity or partition coefficient breaks it, and the check fires
        G13 rather than returning a number that looks reasonable and is not.

        Returns a namespace carrying ``CP0max``, ``CFmax``, ``PR``, the guess
        actually used, the tier, and the provenance of each input.
        """
        from types import SimpleNamespace
        log = self._gatelog if log is None else log

        # --- 1. MAE -> CFmax
        if MAE is None:
            thr = self.binding_threshold(route=route, strict=strict, log=log)
        else:
            thr = MAE if isinstance(MAE, threshold) else threshold(
                float(MAE), route, "MAE (given)", "T2", ADJUDICATED)
        if thr is None:
            return None
        CFmax = self.maxconcentration(thr=thr, strict=strict, log=log)
        if CFmax is None or not np.isfinite(CFmax) or CFmax <= 0:
            _raise("G7-route",
                   "the acceptable exposure did not invert to a usable "
                   "concentration in the product.",
                   "Check the exposure parameters and the route.", True, log)

        # --- 2. one run at the guess
        PR, cp0, cfg, prov = self.potentialrelease(
            packaging, CP0guess=CP0guess, tier=tier, strict=strict, log=log)

        if PR <= 0:
            # No transfer at all: nothing the wall contains can reach the
            # product. Reporting "infinite" is honest; silently returning a
            # huge number is not.
            _raise("G6-massbasis",
                   "the potential release is zero, so no initial concentration "
                   "is limiting.",
                   "Check that the wall is in contact and carries a thickness; "
                   "otherwise the specification is genuinely unbounded.",
                   False, log)
            return SimpleNamespace(CP0max=float("inf"), CFmax=CFmax, PR=0.0,
                                   CP0guess=cp0, CFguess=cfg, tier=tier,
                                   threshold=thr, provenance=prov,
                                   note="no transfer: unbounded")

        # --- 3. verify the linearity the inversion rests on
        if verify and tier == "M0":
            factor = 10.0
            PR2, _, _, _ = self.potentialrelease(
                packaging, CP0guess=cp0 * factor, tier=tier,
                strict=strict, log=log)
            if not np.isclose(PR2, PR, rtol=max(rtol, 1e-9)):
                _raise("G13-linearity",
                       "PR = CF/CP0 changed from %.6g to %.6g when the guess "
                       "was scaled by %g; the transfer is not linear in CP0 "
                       "and the inversion does not hold."
                       % (PR, PR2, factor),
                       "Solve the inverse problem directly, or use a model "
                       "with concentration-independent D and k.", strict, log)

        CP0max = CFmax / PR
        return SimpleNamespace(
            CP0max=CP0max, CFmax=CFmax, PR=PR,
            CP0guess=cp0, CFguess=cfg, tier=tier, threshold=thr,
            provenance="%s; MAE=%s (%s, %s)"
                       % (prov, thr.identity, thr.tier, thr.origin),
            note=None)

    def binding_threshold(self, sub=None, route="systemic", strict=True,
                          log=None):
        """
        The most informed threshold available for the substance.

        Ladder, ordered by information: T2 the substance's own SML converted to
        an implied TDI where that conversion is legitimate; T1 the Cramer TTC;
        T0 the genotoxic TTC. Higher tier means **better informed** -- the
        ladder is explicitly artificial.

        It does **not** mean less conservative. Conservatism is what knowledge
        replaces, and where the replacement lands is an outcome, not a
        guarantee: bisphenol A moves from a Cramer III TTC of 1.5e-3 to an
        SML-derived TDI of 8.33e-4, i.e. being better informed moved the
        threshold DOWN. Do not promise a direction the ladder does not owe.

        With the `unidentified` sentinel the ladder starts and ends at T0: the
        genotoxic TTC needs no identity, which is what makes one bound cover
        the whole uncharacterised population.
        """
        log = self._gatelog if log is None else log

        # An explicitly set MAE outranks the whole ladder. What an organisation
        # decides -- an internal limit, a customer specification, a DNEL from a
        # dossier -- is as binding as anything derivable, and often tighter.
        mae = getattr(self, "_MAE", None)
        if mae is not None:
            if mae.route != route:
                _raise("G7-route",
                       "the MAE was set on the %s route; it cannot be read on "
                       "the %s route." % (mae.route, route),
                       "Set an MAE for this route, or assess the other one.",
                       strict, log)
                return None
            return mae

        sub = self._assessedsubstance if sub is None else sub
        if sub is None:
            _raise("G2-validity", "no substance is attached to this product.",
                   "Inject one: `substance %% product`. For the "
                   "uncharacterised population, declare it explicitly: "
                   "`product.declare_unidentified()`.", True, log)
        if _is_unidentified(sub):
            if route != "systemic":
                _raise("G7-route",
                       "no site-of-contact threshold exists for an "
                       "unidentified substance; a systemic value cannot stand "
                       "in for it.",
                       "Identify the substance, or assess the systemic route.",
                       strict, log)
                return None
            return threshold(TTC_UGKGBWDAY[0] * 1e-3, "systemic",
                             TTC_NAME[0], "T0", DEFAULTED,
                             note="no identity claimed; the uncharacterised "
                                  "population is treated as genotoxic")
        if route != "systemic":
            if not self._gate_route(route, sub, strict, log):
                return None
        t2 = _sml_to_tdi(sub, strict=False, log=log)
        if t2 is not None:
            return t2
        return _ttc_of(sub, strict=strict, log=log)

    def assess(self, CF=None, CFunit="mg/kg", route="systemic",
               censoring=None, strict=False):
        """
        Assess the attached substance against this product.

        Returns a `cosmeticassessment` carrying the dose, the **binding
        threshold with its identity and tier**, the margin, the provenance of
        every input and any gate that fired. `strict=False` by default so that
        a sweep over many substances records refusals instead of stopping; set
        `strict=True` for a single defensible case.
        """
        # Gates raised while the transfer was computed -- G9, G20, G18, G6 --
        # belong to the assessment that consumes the result. `totaltransfer`
        # and `unitarypartition` record them on the instance, so seed the log
        # from there: a refusal that fired during transfer and then vanished
        # from the verdict would be worse than one that never fired at all.
        log = list(getattr(self, "_gatelog", []) or [])
        sub = self._assessedsubstance
        prov = {}
        if sub is None:
            _raise("G12-substance", "no substance is attached to this product.",
                   "Inject one: `substance %% product`, or declare the "
                   "uncharacterised population with "
                   "`product.declare_unidentified()`.", True, log)

        # G12 -- substances are consumed through loadpubchem only.
        if not self._gate_substance(sub, strict=False, log=log):
            return cosmeticassessment(sub, self.name, route, None, None,
                                      "not-assessable", prov, log,
                                      note="the substance is not a migrant (G12)")

        # G2 -- validity domain. A UVCB has no transfer model.
        if not self._gate_validity(sub, strict=False, log=log):
            return cosmeticassessment(sub, self.name, route, None, None,
                                      "not-assessable", prov, log,
                                      note="no molecular weight: UVCB, polymer "
                                           "or mixture (G2)")

        # G3 -- where did k come from?
        prov["k"] = self._gate_kprovenance(strict=False, log=log)

        # G19 -- who controls the parameters that went in. Recorded, not fatal:
        # the arithmetic is still worth doing, provided the result cannot be
        # mistaken for one built on the customer's real contact conditions.
        report = self._gate_control(strict=False, log=log)
        asks = report["asks"]
        prov["control"] = ("complete" if not asks
                           else "%d defaulted" % len(asks))

        # Which diffusivity produced the number, when the envelope replaced the
        # model-selected one. A D that came from somewhere other than the
        # layer's own model must say so.
        do = getattr(self, "_Dorigin", None)
        if do and do.get("applied"):
            _, cname, Dv, model = do["applied"][0]
            prov["D"] = "%s bounding %.3g m2/s [%s]" % (cname, Dv, model)
        elif not self.useboundingD:
            prov["D"] = "model-selected (boundingD disabled)"

        # A9 -- state the wall state assumed, whenever a transfer produced CF.
        ws = getattr(self, "_wallstate", None)
        if ws:
            prov["wall"] = "%s%s" % (
                ws["class"],
                "" if ws["Tg_degC"] is None else " (Tg=%g degC)" % ws["Tg_degC"])

        # Concentration, in the order the pipeline produces it:
        #   explicit argument  >  the last simulation  >  CF0 (measured entry)
        if CF is None:
            if self.hassimulation and self.lastsimulation is not None:
                CF = float(np.asarray(self.lastsimulation.CF).ravel()[-1])
                prov["CF"] = SIMULATED
            else:
                CF = self._scalar(self.CF0)
                if CF:
                    prov["CF"] = getattr(self, "_CForigin", None) or MEASURED
                else:
                    prov["CF"] = DEFAULTED
        else:
            prov["CF"] = MEASURED

        # G15 -- absence of evidence is not evidence of absence.
        #
        # `CF0` inherits the foodlayer default of array([0]). With a substance
        # attached, nothing measured and nothing simulated, the arithmetic runs
        # to 0 > threshold == False and returns PASS. In a module whose whole
        # argument is refusal rather than approximation, that is the one case
        # that most deserves refusal.
        #
        # A genuine zero is not a real case: a measured zero is a non-detect,
        # and a non-detect enters at DL or DL/2 by G4.
        if CF is None or float(CF) == 0.0:
            _raise("G15-nodata",
                   "no concentration is available: nothing was measured, "
                   "nothing was simulated, and none was supplied.",
                   "Measure it and use `setCF(value, origin=...)`; run the "
                   "pipeline; pass `CF=`; or, for a non-detect, give the "
                   "detection limit and the convention "
                   "(`assess(censoring='DL')`).", False, log)
            return cosmeticassessment(sub, self.name, route, None, None,
                                      "not-assessable", prov, log, CF=CF,
                                      note="no concentration data (G15)",
                                      asks=asks)

        if censoring is not None:
            if censoring not in ("DL/2", "DL"):
                _raise("G4-censoring",
                       "unknown censoring convention %r." % censoring,
                       "Use 'DL/2' or 'DL'.", True, log)
            # G16 -- the two censoring conventions compose silently.
            #
            # A value entered pre-halved (`CF0 = DL/2`) and then passed through
            # `censoring='DL/2'` is halved twice: a dose four times too low,
            # with `CF=censored (DL/2)` in the provenance looking entirely
            # correct. Only one of the two may carry the convention.
            if str(prov.get("CF", "")).startswith(CENSORED):
                _raise("G16-doublecensoring",
                       "the concentration is already recorded as censored "
                       "(%s) and `censoring=%r` would apply the convention a "
                       "second time." % (prov["CF"], censoring),
                       "Canonical form: store the RAW detection limit and let "
                       "`assess(censoring=...)` apply the convention once. "
                       "Do not pre-halve.", False, log)
            else:
                CF = CF * (0.5 if censoring == "DL/2" else 1.0)
                prov["CF"] = "%s (%s)" % (CENSORED, censoring)

        thr = self.binding_threshold(sub, route=route, strict=strict, log=log)
        if thr is None:
            return cosmeticassessment(sub, self.name, route, None, None,
                                      "not-assessable", prov, log, CF=CF,
                                      note="no threshold on the %s route (G7)"
                                           % route, asks=asks)
        prov["threshold"] = thr.origin

        # G4 -- a censored T0 bound must not be halved
        if censoring == "DL/2" and thr.tier == "T0":
            _raise("G4-censoring",
                   "the T0 bound was halved as if it were a central estimate.",
                   "For an uncharacterised substance assessed at the genotoxic "
                   "TTC, count the non-detect at DL, not DL/2.", False, log)

        dose = (self.systemicdose(CF, CFunit, strict, log) if route == "systemic"
                else self.localdose(CF, CFunit, strict, log))
        if dose is None:
            return cosmeticassessment(sub, self.name, route, None, thr,
                                      "not-assessable", prov, log, CF=CF,
                                      note="the dose could not be computed (G1)",
                                      asks=asks)
        verdict = "exceed" if dose > thr.value else "pass"
        return cosmeticassessment(sub, self.name, route, dose, thr, verdict,
                                  prov, log, CF=CF, asks=asks)

    def refine(self, assessment):
        """
        G8 -- stop on pass.

        Refinement without a trigger is a defect, not diligence. When a case
        fails, the exit is a **branch, not an escalation**: refine one rung, or
        verify experimentally. Both are offered; neither is automatic.
        """
        if assessment.verdict == "pass":
            raise cosmeticgate(
                "G8-stoponpass",
                "the case passes at tier %s with a margin of %.3g; there is "
                "nothing to refine." % (assessment.threshold.tier,
                                        assessment.margin or float("nan")),
                "Refine only when a problem is found.")
        return ("Two exits, of equal standing:\n"
                "  (a) refine ONE rung -- name the assumption relaxed;\n"
                "  (b) verify experimentally -- measure the product.\n"
                "A tier-M3 result on defaulted parameters is weaker than the "
                "M0 bound it replaced.")


# %% Product categories ------------------------------------------------------
# The three worked references of the CosPaTox dossier. Each inherits its
# prescribed simulant, so the simulant doctrine is carried by the type: EtOH 95
# for hydrophobic products, EtOH 50 for hydrophilic products, detergents and
# home care. Inheriting the simulant is also what supplies `_chemicalsubstance`,
# without which `_compute_kmodel` cannot evaluate and `k` silently falls back to
# the affinity-mixin default (G3).
#
# The exposure figures below are PLACEHOLDERS pending the values of the
# published CosPaTox deliverables and the SCCS Notes of Guidance. They are
# declared as class defaults precisely so they can be overridden per instance,
# and they are tagged DEFAULTED in every result that uses them.

class bodylotion(cosmetic, ethanol):
    """Body lotion -- leave-on, hydrophobic, EtOH 95 doctrine."""
    name = "body lotion"
    description = "leave-on body lotion (CosPaTox use case, adult)"
    level = "user"
    _usage = "leave-on"
    _target = "adult"
    amountapplied, amountappliedUnits = check_units((8.0, "g"))
    frequency, frequencyUnits = check_units((1.0, "1/day"))
    retention, retentionUnits = check_units((1.0, NoUnits))
    exposedarea, exposedareaUnits = check_units((16000.0, "cm**2"))
    bodyweight, bodyweightUnits = check_units((60.0, "kg"))
    volume, volumeUnits = check_units((200, "mL"))


class shampoo(cosmetic, ethanol50):
    """Shampoo -- rinse-off, hydrophilic, EtOH 50 doctrine."""
    name = "shampoo"
    description = "rinse-off shampoo (CosPaTox use case, adult)"
    level = "user"
    _usage = "rinse-off"
    _target = "adult"
    amountapplied, amountappliedUnits = check_units((10.0, "g"))
    frequency, frequencyUnits = check_units((1.0, "1/day"))
    retention, retentionUnits = check_units((0.01, NoUnits))
    exposedarea, exposedareaUnits = check_units((1440.0, "cm**2"))
    bodyweight, bodyweightUnits = check_units((60.0, "kg"))
    volume, volumeUnits = check_units((250, "mL"))


class washinggel(cosmetic, ethanol50):
    """Washing gel -- rinse-off, hydrophilic, EtOH 50 doctrine, child target."""
    name = "washing gel"
    description = "rinse-off washing gel (CosPaTox use case, child)"
    level = "user"
    _usage = "rinse-off"
    _target = "child"
    amountapplied, amountappliedUnits = check_units((5.0, "g"))
    frequency, frequencyUnits = check_units((1.0, "1/day"))
    retention, retentionUnits = check_units((0.01, NoUnits))
    exposedarea, exposedareaUnits = check_units((8000.0, "cm**2"))
    bodyweight, bodyweightUnits = check_units((10.0, "kg"))
    volume, volumeUnits = check_units((250, "mL"))


# %% Recyclate initial concentration -----------------------------------------

# Conventional worst-case contamination of a post-consumer recyclate, mg/kg,
# before decontamination. Polyolefins sorb apolar migrants far more readily
# than PET and decontaminate far less effectively, hence the order of magnitude
# between them.
#
# !! The PET figure must be traced to its source (understood to be the EFSA
# challenge-test convention carried into EU 2022/1616) before it is quoted in
# any deliverable. It is a regulatory-adjacent number, not a model output.
C_REF_MGKG = {
    "PET": 3.0,
    "polyolefin": 30.0,
}


def recyclate_C0(family="PET", decontamination=0.0):
    """
    Initial concentration of a recyclate, mg/kg.

        C0 = Cref * (1 - eta)

    A recyclate is **not** a polymer state but a contamination state: the same
    polymer carrying a history. So no new layer class is created -- only `C0`
    is set on an existing one (see G9 on the `r` prefix).

    `decontamination` is the efficiency certified by the recycling process, and
    is therefore the single lever available to the customer.
    """
    if family not in C_REF_MGKG:
        raise cosmeticgate("G11-recyclate",
                           "unknown polymer family %r." % family,
                           "Use one of: %s" % ", ".join(sorted(C_REF_MGKG)))
    eta = float(decontamination)
    if not 0.0 <= eta < 1.0:
        raise cosmeticgate("G11-recyclate",
                           "decontamination efficiency %g is outside [0, 1)." % eta,
                           "Give a fraction, not a percentage.")
    return C_REF_MGKG[family] * (1.0 - eta)


def required_decontamination(C0max, family="PET"):
    """
    The decontamination efficiency the recycler must certify.

        eta_min = 1 - C0max / Cref

    The counterpart of `cosmetic.maxdetectionlimit`: rather than a verdict on a
    material, a specification the supplier can act on. Returns 0 when the
    untreated recyclate already complies, and None when no efficiency suffices.
    """
    if family not in C_REF_MGKG:
        raise cosmeticgate("G11-recyclate",
                           "unknown polymer family %r." % family,
                           "Use one of: %s" % ", ".join(sorted(C_REF_MGKG)))
    ref = C_REF_MGKG[family]
    if C0max is None or C0max <= 0:
        return None
    if C0max >= ref:
        return 0.0
    return 1.0 - C0max / ref


# %% Help --------------------------------------------------------------------

def help_cosmetic():
    """Print the cosmetic class tree and the exposure chain."""
    print(__doc__)
    print("Categories:")
    for c in (bodylotion, shampoo, washinggel):
        print("  %-12s %-10s %-6s bw=%g kg  retention=%g"
              % (c.name, c._usage, c._target,
                 np.asarray(c.bodyweight).ravel()[0],
                 np.asarray(c.retention).ravel()[0]))
    print("\nGates:")
    help_gates()


if __name__ == '__main__':
    from patankar.loadpubchem import migrantToxtree
    from patankar.layer import gPET, wPET, rPET, LDPE
    from patankar.geometry import Packaging3D

    print("=" * 72)
    print("SFPPy cosmetic module -- self test")
    print("=" * 72)

    # --- 1. the three categories, measured entry (non-detect at DL = 1 mg/kg)
    m = migrantToxtree("toluene")
    print("\n[1] Measured entry -- CF0 = DL/2 with DL = 1 mg/kg\n")
    for product in (bodylotion, shampoo, washinggel):
        p = m % product()
        p.CF0 = 0.5
        print("--- %s (%s, %s) ---" % (p.name, p.usage, p.target))
        repr(p.assess())                       # SFPPy repr idiom: prints
        print("   Cmax   : %.4g mg/kg" % p.maxconcentration(strict=False))
        print("   DLmax  : %.4g mg/kg (DL/2 convention)"
              % p.maxdetectionlimit(strict=False))
        print("   DLmax  : %.4g mg/kg (forced DL, T0 bound)"
              % p.maxdetectionlimit(censoring="DL", strict=False))
        print()

    # --- 2. predictive entry, M0 bound then the standard SFPPy pipeline
    print("[2] Predictive entry -- recycled PET bottle\n")
    C0 = recyclate_C0("PET", decontamination=0.90)
    print("    C0 = 3 mg/kg x (1 - 0.90) = %.3g mg/kg" % C0)
    bottle = Packaging3D('bottle',
                         body_radius=(3, 'cm'), body_height=(15, 'cm'),
                         neck_radius=(1, 'cm'), neck_height=(3, 'cm'))

    # M0 -- mass balance, no solver
    p = m % bodylotion()
    p = p << bottle                            # inherit volume and surface area
    # A cosmetic swells the wall: the plasticized state is the conservative
    # default and the glassy one is what has to be argued (G20).
    CF0bound = p.totaltransfer(wPET(l=(500, "um"), C0=C0), strict=False)
    print("    V_F = %.4g L, A = %.4g dm2"
          % (float(np.asarray(p.volume).ravel()[0]) * 1e3,
             float(np.asarray(p.surfacearea).ravel()[0]) * 1e2))
    print("    M0  CF = %.4g mg/kg  (total transfer bound)" % CF0bound)
    repr(p.assess(CF=CF0bound))

    # M3 -- the standard pipeline: substance % product << geometry >> layer >> product
    print("\n    M3 through the standard pipeline:")
    q = bodylotion()
    wall = wPET(l=(500, "um"), C0=C0)
    m % q << bottle >> wall >> q               # runs; q.lastsimulation is set
    print("    contact: %g days at %g degC"
          % (float(np.asarray(q.contacttime).ravel()[0]) / 86400.0,
             float(np.asarray(q.contacttemperature).ravel()[0])))
    print("    M3  CF = %.4g mg/kg  (Patankar solver)"
          % float(np.asarray(q.lastsimulation.CF).ravel()[-1]))
    repr(q.assess())                           # CF read from lastsimulation
    print("    M0 is a bound on M3: %s"
          % ("yes" if CF0bound >= float(np.asarray(q.lastsimulation.CF).ravel()[-1])
             else "NO -- investigate"))

    cmax = p.maxconcentration(strict=False)
    C0max = cmax * float(np.asarray(p.density).ravel()[0]) \
        * float(np.asarray(p.volume).ravel()[0]) \
        / (1350.0 * float(np.asarray(p.surfacearea).ravel()[0]) * 500e-6)
    eta = required_decontamination(C0max, "PET")
    print("    required decontamination: eta_min = %s"
          % ("none needed" if eta == 0 else "%.3g %%" % (100 * eta)))

    # --- 3. the gates
    print("\n[3] Gates -- refusals, not approximations\n")

    def _try(label, fn):
        try:
            fn()
            print("    %-22s NOT FIRED (unexpected)" % label)
        except cosmeticgate as e:
            print("    %-22s fired: %s" % (label, str(e).splitlines()[0]))

    # G9 now refuses rPET below Tg only where the contact was declared
    # NON-swelling. For a cosmetic the rubbery parameterisation below Tg is the
    # recommended plasticization surrogate, not a state error.
    _dry = m % bodylotion()
    _dry.plasticizing = False
    _try("G9 polymer state (dry)", lambda: _dry.totaltransfer(
        rPET(l=(500, "um"), C0=C0, T=(25, "degC")), strict=True))
    print("    %-22s allowed for a swelling medium: %s"
          % ("rPET below Tg",
             "yes -- the plasticization surrogate (G9)"
             if p._gate_polymerstate(rPET(l=(500, "um"), C0=C0, T=(25, "degC")),
                                     strict=False) else "no"))
    _try("G20 plasticization", lambda: p.totaltransfer(
        gPET(l=(500, "um"), C0=C0, T=(23, "degC")), strict=True))
    _try("G19 control class", lambda: bodylotion()._gate_control(strict=True))
    _try("G21 Tg no-op", lambda: (migrantToxtree("bisphenol A") % bodylotion())
         .totaltransfer(wPET(l=(500, "um"), C0=C0, T=(23, "degC"),
                             Tg=(30, "degC")), strict=True))
    # G20 inverted: for a very small migrant the GLASSY state comes out faster,
    # because the two states are not evaluated by the same diffusion model.
    _try("G20 ordering crossed", lambda: (migrantToxtree("formaldehyde")
                                          % bodylotion())
         .totaltransfer(gPET(l=(500, "um"), C0=C0, T=(23, "degC")), strict=True))
    _try("G17 reaction direction", lambda: p.completereaction(
        1.0, consume_parent=True, strict=True))
    _try("G18 partition bound", lambda: bodylotion(
        volume=(30, "mL"), surfacearea=(60, "cm**2")).unitarypartition(
            wPET(l=(2, "mm"), C0=C0), strict=True))
    _try("G1 units", lambda: p.systemicdose(1.0, "mg per kilo", strict=True))
    _try("G4 censoring", lambda: p.maxdetectionlimit(censoring="half"))
    _try("G11 recyclate family", lambda: recyclate_C0("PVC"))
    _try("G11 decontamination", lambda: recyclate_C0("PET", decontamination=90))

    a = (m % shampoo()); a.CF0 = 1e-6
    _try("G8 stop on pass", lambda: a.refine(a.assess()))

    _try("G12 not a migrant", lambda: shampoo()._gate_substance("toluene"))

    print("\n    G2 validity domain -- a UVCB never becomes a migrant at all:")
    sub, why = as_migrant("68038-01-7")        # Trialkylamin
    print("      as_migrant -> %s" % ("resolved" if sub else "REFUSED"))
    print("      %s" % why)

    print("\n" + "=" * 72)
