"""
Composition-order tests for food classes.

A concrete medium is built by multiple inheritance, mixing a food property
class (`realfood`, `simulant`, ...) with a contact condition (`ambient`,
`frozen`, ...). Both `foodlayer` and every contact condition declare
`contacttime` and `contacttemperature`, so **which one wins is decided by the
MRO**, i.e. by the order the bases are listed:

    class F(ambient, realfood, ...)   -> ambient wins    (contact-first)
    class F(realfood, ambient, ...)   -> foodlayer wins  (food-first)

That is ordinary Python attribute resolution, not a quirk: for a food-first
declaration C3 places the whole `realfood -> foodproperty -> foodlayer` chain
ahead of `ambient`, so `foodlayer`'s defaults outrank the contact condition.

These tests pin that contract, and in particular pin the **invariant that
`contacttime` and `contacttemperature` always come from the same class**.
That invariant was silently violated for years: `foodlayer` declared the
temperature as `contactemperature` (missing a `t`), so it never competed for
the name and the contact condition's temperature survived in either order,
while its contact time did not. Nothing in the suite noticed. A future
misspelling of either name breaks `test_time_and_temperature_agree` below.

@project: SFPPy - Safe Food Packaging in Python
@author: Olivier Vitrac
@license: MIT
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from patankar.food import (foodlayer, realfood, ambient, frozen, hotfilled,
                           perfectlymixed, fat, parametersWithUnits,
                           paramaterNamesWithUnits)


CONTACT_CONDITIONS = [ambient, frozen, hotfilled]


def declared(cls, attr):
    """Value a class declares itself for `attr` (not inherited), as a scalar."""
    return float(np.asarray(cls.__dict__[attr]).ravel()[0])


def value(obj, attr):
    """Instance value of `attr` as a scalar."""
    return float(np.asarray(getattr(obj, attr)).ravel()[0])


def make(contact, contact_first):
    """Build a concrete medium mixing `realfood` with a contact condition."""
    bases = ((contact, realfood, perfectlymixed, fat) if contact_first
             else (realfood, contact, perfectlymixed, fat))
    name = "%s_%s" % (contact.__name__, "first" if contact_first else "last")
    return type(name, bases, {"name": name, "level": "user"})()


class TestContactConditionPrecedence:
    """Listing the contact condition first is what makes its values win."""

    @pytest.mark.parametrize("contact", CONTACT_CONDITIONS)
    @pytest.mark.parametrize("attr", ["contacttime", "contacttemperature"])
    def test_contact_first_wins(self, contact, attr):
        """contact-first: the contact condition supplies the value."""
        obj = make(contact, contact_first=True)
        assert_allclose(value(obj, attr), declared(contact, attr), rtol=1e-12)

    @pytest.mark.parametrize("contact", CONTACT_CONDITIONS)
    @pytest.mark.parametrize("attr", ["contacttime", "contacttemperature"])
    def test_food_first_takes_foodlayer_defaults(self, contact, attr):
        """food-first: foodlayer outranks the contact condition."""
        obj = make(contact, contact_first=False)
        assert_allclose(value(obj, attr), declared(foodlayer, attr), rtol=1e-12)

    @pytest.mark.parametrize("contact", CONTACT_CONDITIONS)
    @pytest.mark.parametrize("contact_first", [True, False])
    def test_time_and_temperature_agree(self, contact, contact_first):
        """
        The harmonisation invariant.

        Whichever class wins must win for BOTH contact time and contact
        temperature. If one of the two is ever misspelled again it will stop
        competing in the MRO, the two will be sourced from different classes,
        and this test fails.
        """
        obj = make(contact, contact_first=contact_first)
        winner = contact if contact_first else foodlayer
        assert_allclose(value(obj, "contacttime"),
                        declared(winner, "contacttime"), rtol=1e-12)
        assert_allclose(value(obj, "contacttemperature"),
                        declared(winner, "contacttemperature"), rtol=1e-12)

    @pytest.mark.parametrize("contact", CONTACT_CONDITIONS)
    def test_order_actually_matters(self, contact):
        """
        Guard against the contract quietly becoming order-independent.

        The chosen conditions differ from the foodlayer defaults, so the two
        orders must give different answers. If they ever agree, either the
        precedence rule changed or a default drifted.
        """
        first = make(contact, contact_first=True)
        last = make(contact, contact_first=False)
        assert (value(first, "contacttime") != value(last, "contacttime")
                or value(first, "contacttemperature")
                != value(last, "contacttemperature"))


class TestExplicitOverride:
    """Constructor keywords beat class defaults, whatever the order."""

    @pytest.mark.parametrize("contact_first", [True, False])
    def test_kwargs_win_over_both(self, contact_first):
        """
        This is the path the food-tree widget uses: it passes contacttime and
        contacttemperature explicitly, which is why the GUI is insensitive to
        the composition order.
        """
        bases = ((ambient, realfood, perfectlymixed, fat) if contact_first
                 else (realfood, ambient, perfectlymixed, fat))
        cls = type("Explicit", bases, {"name": "explicit", "level": "user"})
        obj = cls(contacttime=(3, "days"), contacttemperature=(7, "degC"))
        assert_allclose(value(obj, "contacttime"), 3 * 86400.0, rtol=1e-12)
        assert_allclose(value(obj, "contacttemperature"), 7.0, rtol=1e-12)


class TestUnitsNaming:
    """Units labels follow key+'Units' and are protected from reassignment."""

    @pytest.mark.parametrize("key", sorted(parametersWithUnits))
    def test_units_label_is_reachable(self, key):
        """
        __repr__ looks the label up as f"{key}Units"; a label declared under
        any other spelling is dead and the display silently falls back to the
        registry string.
        """
        obj = foodlayer()
        if not hasattr(obj, key):          # k/k0 are mutually exclusive
            pytest.skip("%s not defined on a bare foodlayer" % key)
        assert hasattr(obj, key + "Units"), \
            "%sUnits is missing: the label is unreachable from __repr__" % key

    @pytest.mark.parametrize("unitname", sorted(set(paramaterNamesWithUnits)))
    def test_units_cannot_be_reassigned(self, unitname):
        """Units are SI and fixed; kwargs must not overwrite them."""
        obj = foodlayer(**{unitname: "furlong"})
        assert getattr(obj, unitname, None) != "furlong"


class TestDeprecatedAliases:
    """The renamed attributes keep working, with a DeprecationWarning."""

    @pytest.mark.parametrize("old,new", [
        ("contactemperature", "contacttemperature"),
        ("CF0units", "CF0Units"),
        ("contacttime_units", "contacttimeUnits"),
    ])
    def test_alias_reads_through(self, old, new):
        """The old name still resolves to the canonical attribute."""
        obj = foodlayer()
        with pytest.warns(DeprecationWarning):
            got = getattr(obj, old)
        expected = getattr(obj, new)
        if isinstance(expected, np.ndarray):
            assert_allclose(np.asarray(got).ravel(),
                            np.asarray(expected).ravel(), rtol=1e-12)
        else:
            assert got == expected

    @pytest.mark.parametrize("old,canonical", [
        ("CF0units", "CF0Units"),
        ("contacttime_units", "contacttimeUnits"),
        ("contactemperatureUnits", "contacttemperatureUnits"),
    ])
    def test_deprecated_units_write_does_not_raise(self, old, canonical):
        """
        Assigning a deprecated units label warns and is ignored, rather than
        raising: before the rename these names were silently accepted, so a
        hard failure would be a regression for existing code.
        """
        obj = foodlayer()
        with pytest.warns(DeprecationWarning):
            setattr(obj, old, "furlong")
        assert getattr(obj, canonical) != "furlong"
