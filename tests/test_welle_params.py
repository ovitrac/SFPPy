"""
Regression tests for the Welle diffusivity model (patankar.property.Dwelle).

The parameters are checked against the original publications, and the model is
validated against the 182 diffusion coefficients measured in GPPS and HIPS by
Welle (2021), below and above Tg.

Sources
-------
- Ewender J., Welle F. (2022). A new method for the prediction of diffusion
  coefficients in poly(ethylene terephthalate) - Validation data.
  Packag. Technol. Sci. 35(5):405-413, Table 1 (gPET).
- Welle F. (2021). Diffusion coefficients and activation energies of diffusion
  of organic molecules in polystyrene below and above glass transition
  temperature. Polymers 13(8):1317, Table 8 (PS, rPS, HIPS, rHIPS),
  Tables 3-5 (measured D), Tables 6-7 (molecular volumes), Tg = 100 degC.

@project: SFPPy - Safe Food Packaging in Python
@author: Olivier Vitrac
@license: MIT
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from patankar.property import Dwelle

# =============================================================================
# Published parameters (a in 1/K, b in cm2/s, c in A3, d in 1/K)
# =============================================================================
PUBLISHED = {
    "gPET":  {"a": 1.93e-3, "b": 2.37e-6, "c": 11.1,  "d": 1.50e-4},  # Ewender & Welle 2022, Table 1
    "PS":    {"a": 2.59e-3, "b": 7.38e-9, "c": 55.71, "d": 2.73e-5},  # Welle 2021, Table 8, GPPS < Tg
    "rPS":   {"a": 2.44e-3, "b": 6.46e-8, "c": 25.51, "d": 7.55e-5},  # Welle 2021, Table 8, GPPS > Tg
    "HIPS":  {"a": 2.55e-3, "b": 9.21e-9, "c": 73.28, "d": 2.04e-5},  # Welle 2021, Table 8, HIPS < Tg
    "rHIPS": {"a": 2.46e-3, "b": 2.07e-7, "c": 45.00, "d": 3.57e-5},  # Welle 2021, Table 8, HIPS > Tg
}

# =============================================================================
# Measured diffusion coefficients, Welle (2021) - D in cm2/s, T in degC
# =============================================================================
VOLUME = {  # molinspiration van der Waals volumes, A3 (Tables 6-7)
    "n-octane": 146.6, "n-decane": 180.2, "n-dodecane": 213.8, "n-tetradecane": 247.3,
    "n-hexadecane": 281.0, "n-octadecane": 314.6, "styrene": 111.8, "acetone": 64.7,
    "ethyl acetate": 90.6, "toluene": 100.6, "chlorobenzene": 97.6,
    "phenyl cyclohexane": 174.0, "benzophenone": 174.4, "methanol": 37.2,
    "ethanol": 54.0, "1-propanol": 70.8, "1-butanol": 87.6,
}
MEASURED = {
    "GPPS": {
        # Table 3 (desorption)
        "n-octane":      {80: 3.4e-12, 85: 1.4e-11, 90: 3.8e-11, 95: 1.1e-10, 100: 3.8e-10, 105: 8.6e-10, 110: 2.0e-9, 115: 4.0e-9},
        "n-decane":      {80: 2.8e-13, 85: 2.2e-12, 90: 8.6e-12, 95: 3.3e-11, 100: 1.6e-10, 105: 4.0e-10, 110: 1.1e-9, 115: 2.4e-9},
        "n-dodecane":    {80: 3.0e-14, 85: 4.3e-13, 90: 2.4e-12, 95: 1.2e-11, 100: 7.5e-11, 105: 2.1e-10, 110: 6.3e-10, 115: 1.5e-9},
        "n-tetradecane": {80: 4.7e-15, 85: 1.1e-13, 90: 8.9e-13, 95: 6.3e-12, 100: 4.2e-11, 105: 1.3e-10, 110: 4.1e-10, 115: 9.9e-10},
        "n-hexadecane":  {80: 1.5e-15, 85: 4.1e-14, 90: 4.3e-13, 95: 3.4e-12, 100: 2.6e-11, 105: 7.6e-11, 110: 2.6e-10, 115: 5.2e-10},
        "n-octadecane":  {80: 1.6e-15, 85: 1.9e-14, 90: 2.8e-13, 95: 1.6e-12, 100: 1.8e-11},
        "styrene":       {80: 1.5e-11, 85: 3.3e-11, 90: 7.1e-11, 95: 1.8e-10, 100: 4.1e-10, 105: 9.7e-10, 110: 1.6e-9, 115: 3.3e-9},
        # Table 4 (desorption)
        "acetone":            {85: 9.8e-9, 90: 1.1e-8, 95: 1.2e-8, 100: 1.0e-8, 105: 1.2e-8},
        "ethyl acetate":      {85: 1.1e-9, 90: 1.8e-9, 95: 3.3e-9, 100: 5.2e-9, 105: 9.8e-9},
        "toluene":            {85: 1.0e-10, 90: 2.1e-10, 95: 4.9e-10, 100: 1.1e-9, 105: 2.4e-9},
        "chlorobenzene":      {85: 2.0e-10, 90: 3.9e-10, 95: 8.8e-10, 100: 1.9e-9, 105: 4.0e-9},
        "phenyl cyclohexane": {85: 2.1e-13, 90: 1.2e-12, 95: 9.5e-12, 100: 4.1e-11, 105: 1.2e-10},
        "benzophenone":       {85: 1.4e-12, 90: 5.7e-12, 95: 2.8e-11, 100: 9.6e-11, 105: 1.5e-10},
        "styrene (sheet 2)":  {85: 2.9e-11, 90: 6.4e-11, 95: 1.5e-10, 100: 3.5e-10, 105: 8.3e-10},
        # Table 5 (permeation, biaxially oriented film)
        "methanol":   {0: 1.2e-9, 25: 2.2e-9, 40: 3.1e-9, 60: 2.7e-9, 70: 2.1e-9, 80: 2.4e-9, 90: 2.4e-9},
        "ethanol":    {0: 2.8e-11, 25: 1.7e-10, 40: 3.8e-10, 60: 1.0e-9, 70: 1.0e-9, 80: 1.4e-9, 90: 1.6e-9},
        "1-propanol": {40: 2.5e-11, 60: 1.0e-10, 70: 2.0e-10, 80: 3.5e-10, 90: 5.6e-10},
        "1-butanol":  {40: 2.7e-12, 60: 1.7e-11, 70: 4.1e-11, 80: 8.5e-11, 90: 1.7e-10},
    },
    "HIPS": {
        # Table 3 (n-decane: analytical artefacts)
        "n-octane":      {80: 2.5e-12, 90: 2.0e-11, 95: 1.8e-10, 100: 4.0e-10, 105: 1.9e-9, 110: 2.9e-9},
        "n-dodecane":    {80: 8.6e-14, 90: 1.6e-12, 95: 6.5e-12, 100: 3.0e-11, 105: 2.0e-10, 110: 4.9e-10},
        "n-tetradecane": {80: 6.7e-16, 90: 3.1e-14, 95: 1.4e-12, 100: 9.7e-12, 105: 8.2e-11, 110: 2.3e-10},
        "n-hexadecane":  {80: 1.8e-16, 90: 8.1e-15, 95: 5.0e-13, 100: 4.6e-12, 105: 4.6e-11, 110: 1.4e-10},
        # Table 4 (acetone: analytical artefacts)
        "ethyl acetate":      {75: 5.1e-10, 80: 6.9e-10, 85: 1.3e-9, 90: 2.2e-9, 100: 6.3e-9, 105: 1.3e-8, 110: 2.3e-8, 115: 5.4e-8},
        "toluene":            {75: 2.9e-11, 80: 4.6e-11, 85: 1.2e-10, 90: 2.4e-10, 100: 1.1e-9, 105: 3.6e-9, 110: 7.8e-9, 115: 2.0e-8},
        "chlorobenzene":      {75: 4.5e-11, 80: 6.9e-11, 85: 1.7e-10, 90: 3.3e-10, 100: 1.4e-9, 105: 4.2e-9, 110: 8.7e-9, 115: 2.2e-8},
        "phenyl cyclohexane": {80: 6.4e-15, 85: 4.0e-14, 90: 1.5e-13, 100: 5.1e-12, 105: 3.4e-11, 110: 1.4e-10, 115: 5.3e-10},
        "benzophenone":       {80: 4.5e-14, 85: 1.7e-13, 90: 5.5e-13, 100: 1.0e-11, 105: 4.9e-11, 110: 1.1e-10, 115: 2.3e-10},
        "styrene":            {75: 1.8e-11, 80: 2.9e-11, 85: 7.9e-11, 90: 4.0e-11, 100: 8.2e-10, 105: 2.7e-9, 110: 6.2e-9, 115: 1.7e-8},
    },
}
TG = 100.0  # degC, GPPS and HIPS (Welle 2021)


def _key(family, T):
    """Welle parameter set: below Tg (PS, HIPS) or at/above Tg (rPS, rHIPS)."""
    if family == "GPPS":
        return "PS" if T < TG else "rPS"
    return "HIPS" if T < TG else "rHIPS"


def _records():
    for family, substances in MEASURED.items():
        for name, data in substances.items():
            V = VOLUME[name.split(" (")[0]]
            for T, D in data.items():
                yield family, name, V, float(T), D * 1e-4  # cm2/s -> m2/s


def _log_ratios(key):
    """log10(D_predicted / D_measured) for all records of one parameter set."""
    out = []
    with np.errstate(over="ignore", under="ignore"):
        for family, _, V, T, Dmeas in _records():
            if _key(family, T) == key:
                Dpred = Dwelle.evaluate(polymer=key, Vvdw=V, T=T)
                out.append(np.log10(Dpred / Dmeas) if Dpred > 0 and np.isfinite(Dpred) else np.nan)
    return np.array(out)


class TestWelleParameters:
    """The implemented parameters equal the published ones."""

    @pytest.mark.parametrize("key", sorted(PUBLISHED))
    def test_parameters_match_publications(self, key):
        for p in "abcd":
            assert_allclose(Dwelle.welle_data[key][p], PUBLISHED[key][p], rtol=1e-12,
                            err_msg=f"Dwelle.welle_data[{key!r}][{p!r}]")

    def test_measurement_dataset_size(self):
        assert sum(1 for _ in _records()) == 182


class TestWelleReferenceValues:
    """Point values of D, in m2/s."""

    @staticmethod
    def _closed_form(key, V, T):
        """D = 1e-4 b (V/c)^((a - 1/T)/d), in m2/s, from the published parameters."""
        p = PUBLISHED[key]
        return 1e-4 * p["b"] * (V / p["c"]) ** ((p["a"] - 1 / (T + 273.15)) / p["d"])

    @pytest.mark.parametrize("key,V,T,expected", [
        ("rHIPS", 150, 100, 1.2457e-14),   # was 0.0 with d = 2.07e-7
        ("gPET", 100, 40, 2.1582e-18),     # was 2.0671e-18 with b = 2.27e-6
    ])
    def test_reference_values(self, key, V, T, expected):
        assert_allclose(self._closed_form(key, V, T), expected, rtol=1e-4)
        assert_allclose(Dwelle.evaluate(polymer=key, Vvdw=V, T=T), expected, rtol=1e-4)

    @pytest.mark.parametrize("T", [100, 110, 120, 130])
    def test_rHIPS_finite_positive_and_decreasing_with_volume(self, T):
        V = np.linspace(40, 400, 37)
        D = np.array([Dwelle.evaluate(polymer="rHIPS", Vvdw=v, T=T) for v in V])
        assert np.all(np.isfinite(D)) and np.all(D > 0)
        assert np.all(np.diff(D) < 0)
        assert np.all((D > 1e-20) & (D < 1e-8))


class TestWelleAgainstMeasurements:
    """Validation against the 182 diffusion coefficients of Welle (2021).

    Welle's correlations are central estimators (median log10 ratio close to 0);
    they are not upper bounds. Outside the domain V > c (methanol and ethanol in
    GPPS, V < c = 55.71 A3), the correlation over-predicts by several decades.
    """

    @pytest.mark.parametrize("key", ["rPS", "HIPS", "rHIPS"])
    def test_within_one_decade(self, key):
        r = _log_ratios(key)
        assert np.all(np.isfinite(r)), f"{key}: non-finite predictions"
        assert np.mean(np.abs(r) < 1) >= 0.90, f"{key}: {np.mean(np.abs(r) < 1):.0%} within one decade"
        assert abs(np.median(r)) < 0.3, f"{key}: median log10 ratio {np.median(r):+.2f}"

    def test_PS_within_one_decade_when_V_above_c(self):
        c = PUBLISHED["PS"]["c"]
        r = np.array([np.log10(Dwelle.evaluate(polymer="PS", Vvdw=V, T=T) / D)
                      for fam, _, V, T, D in _records() if _key(fam, T) == "PS" and V > c])
        assert np.all(np.abs(r) < 1)
