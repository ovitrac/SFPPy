# -*- coding: utf-8 -*-
"""
chemspipy.objects
~~~~~~~~~~~~~~~~~

Objects returned by ChemSpiPy API methods.

"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import warnings


from .utils import memoized_property


class Compound(object):
    """ A class for retrieving and caching details about a specific ChemSpider record.

    The purpose of this class is to provide access to various parts of the ChemSpider API that return information about
    a compound given its ChemSpider ID. Information is loaded lazily when requested, and cached for future access.
    """

    def __init__(self, cs, record_id):
        """

        :param ChemSpider cs: ``ChemSpider`` session.
        :param int|string record_id: Compound record ID.
        """
        self._cs = cs
        self._record_id = int(record_id)
        # TODO: Allow optional initialize  with a record-type response from the API (kwarg or class method from_dict?).

    def __eq__(self, other):
        return isinstance(other, Compound) and self.csid == other.csid

    def __repr__(self):
        return 'Compound(%r)' % self.csid

    def _repr_png_(self):
        """For IPython notebook, display 2D image."""
        return self.image

    @property
    def record_id(self):
        """Compound record ID.

        :rtype: int
        """
        return self._record_id

    @property
    def csid(self):
        """ChemSpider ID.

        .. deprecated:: 2.0.0
           Use :py:attr:`~chemspipy.objects.Compound.record_id` instead.

        :rtype: int
        """
        warnings.warn('Use record_id instead of csid.', DeprecationWarning)
        return self._record_id

    @property
    def image_url(self):
        """Return the URL of a PNG image of the 2D chemical structure.

        :rtype: string
        """
        return 'http://www.chemspider.com/ImagesHandler.ashx?id=%s' % self.record_id

    @memoized_property
    def _details(self):
        """Request compound info and cache the result."""
        return self._cs.get_details(self.record_id)

    @property
    def molecular_formula(self):
        """Return the molecular formula for this Compound.

        :rtype: string
        """
        return self._details['formula']

    @property
    def smiles(self):
        """Return the SMILES for this Compound.

        :rtype: string
        """
        return self._details['smiles']

    # TODO: Convert tool to get inchi?

    @property
    def stdinchi(self):
        """Return the Standard InChI for this Compound.

        .. deprecated:: 2.0.0
           Use :py:attr:`~chemspipy.objects.Compound.inchi` instead.

        :rtype: string
        """
        warnings.warn('Use inchi instead of stdinchi.', DeprecationWarning)
        return self.inchi

    @property
    def stdinchikey(self):
        """Return the Standard InChIKey for this Compound.

        .. deprecated:: 2.0.0
           Use :py:attr:`~chemspipy.objects.Compound.inchikey` instead.

        :rtype: string
        """
        warnings.warn('Use inchikey instead of stdinchikey.', DeprecationWarning)
        return self.inchikey

    @property
    def inchi(self):
        """Return the InChI for this Compound.

        :rtype: string
        """
        return self._cs.convert(self.mol_2d, 'Mol', 'InChI')

    @property
    def inchikey(self):
        """Return the InChIKey for this Compound.

        :rtype: string
        """
        return self._cs.convert(self.mol_2d, 'Mol', 'InChIKey')

    @property
    def average_mass(self):
        """Return the average mass of this Compound.

        :rtype: float
        """
        return self._details['averageMass']

    @property
    def molecular_weight(self):
        """Return the molecular weight of this Compound.

        :rtype: float
        """
        return self._details['molecularWeight']

    @property
    def monoisotopic_mass(self):
        """Return the monoisotopic mass of this Compound.

        :rtype: float
        """
        return self._details['monoisotopicMass']

    @property
    def nominal_mass(self):
        """Return the nominal mass of this Compound.

        :rtype: float
        """
        return self._details['nominalMass']

    @property
    def common_name(self):
        """Return the common name for this Compound.

        :rtype: string
        """
        return self._details['commonName']

    @memoized_property
    def mol_2d(self):
        """Return the MOL file for this Compound with 2D coordinates.

        :rtype: string
        """
        return self._details['mol2D']

    @memoized_property
    def mol_3d(self):
        """Return the MOL file for this Compound with 3D coordinates.

        :rtype: string
        """
        return self._details['mol3D']

    @memoized_property
    def image(self):
        """Return a 2D depiction of this Compound.

        :rtype: bytes
        """
        return self._cs.get_image(self.record_id)

    @memoized_property
    def external_references(self):
        """Return external references for this Compound.

        :rtype: list[dict]
        """
        return self._cs.get_external_references(self.record_id)
