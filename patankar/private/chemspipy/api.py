# -*- coding: utf-8 -*-
"""
chemspipy.api
~~~~~~~~~~~~~

Core API for interacting with ChemSpider web services.

"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import base64
import logging
import sys
import warnings
import zlib

import requests

from . import __version__, errors
from .objects import Compound
from .search import Results


log = logging.getLogger(__name__)


#: Default API URL.
API_URL = 'https://api.rsc.org'
#: Default API version.
API_VERSION = 'v1'

#: 2D coordinate dimensions
MOL2D = '2d'
#: 3D coordinate dimensions
MOL3D = '3d'
#: Both coordinate dimensions
BOTH = 'both'

#: Ascending sort direction
ASCENDING = 'ascending'
#: Descending sort direction
DESCENDING = 'descending'

#: Record ID sort order
RECORD_ID = 'record_id'
#: CSID sort order (same as RECORD_ID, kept for backwards compatibility)
CSID = 'csid'
#: Mass defect sort order
MASS_DEFECT = 'mass_defect'
#: Molecular weight sort order
MOLECULAR_WEIGHT = 'molecular_weight'
#: Reference count sort order
REFERENCE_COUNT = 'reference_count'
#: Datasource count sort order
DATASOURCE_COUNT = 'datasource_count'
#: Pubmed count sort order
PUBMED_COUNT = 'pubmed_count'
#: RSC count sort order
RSC_COUNT = 'rsc_count'


#: Map sort directions to strings required by REST API.
DIRECTIONS = {
    ASCENDING: 'ascending',
    DESCENDING: 'descending'
}

#: Map sort orders to strings required by REST API.
ORDERS = {
    RECORD_ID: 'recordId',
    CSID: 'recordId',
    MASS_DEFECT: 'massDefect',
    MOLECULAR_WEIGHT: 'molecularWeight',
    REFERENCE_COUNT: 'referenceCount',
    DATASOURCE_COUNT: 'dataSourceCount',
    PUBMED_COUNT: 'pubMedCount',
    RSC_COUNT: 'rscCount'
}

#: All available compound details fields.
FIELDS = [
    'SMILES', 'Formula', 'AverageMass', 'MolecularWeight', 'MonoisotopicMass', 'NominalMass', 'CommonName',
    'ReferenceCount', 'DataSourceCount', 'PubMedCount', 'RSCCount', 'Mol2D', 'Mol3D'
]


class ChemSpider(object):
    """Provides access to the ChemSpider API.

    Usage::

        >>> from chemspipy import ChemSpider
        >>> cs = ChemSpider('<YOUR-API-KEY>')

    """

    def __init__(self, api_key, user_agent=None, api_url=API_URL, api_version=API_VERSION):
        """

        :param string api_key: Your ChemSpider API key.
        :param string user_agent: (Optional) Identify your application to ChemSpider servers.
        :param string api_url: (Optional) API server. Default https://api.rsc.org.
        :param string api_version: (Optional) API version. Default v1.
        """
        log.debug('Initializing ChemSpider')
        self.api_url = api_url
        self.http = requests.session()
        self.http.headers['User-Agent'] = user_agent if user_agent else 'ChemSpiPy/{} Python/{} '.format(
            __version__, sys.version.split()[0]
        )
        self.api_key = api_key
        self.api_version = api_version

    def __repr__(self):
        return f'ChemSpider() with api_key: {self.api_key}'

    def request(self, method, api, namespace, endpoint, params=None, json=None):
        """Make a request to the ChemSpider API.

        :param string method: HTTP method.
        :param string api: Top-level API, e.g. compounds.
        :param string namespace: API namespace, e.g. filter, lookups, records, or tools.
        :param string endpoint: Web service endpoint URL.
        :param dict params: Query parameters to add to the URL.
        :param dict json: JSON data to send in the request body.
        :return: Web Service response JSON.
        :rtype: dict
        """
        # Construct request URL
        url = '{}/{}/{}/{}/{}'.format(self.api_url, api, self.api_version, namespace, endpoint)

        # Set apikey header
        headers = {'apikey': self.api_key}

        log.debug('{} : {} : {} : {}'.format(url, headers, params, json))

        # Make request
        r = self.http.request(method, url, params=params, json=json, headers=headers)

        # Raise exception for HTTP errors
        if not r.ok:
            err = {
                400: errors.ChemSpiPyBadRequestError,
                401: errors.ChemSpiPyAuthError,
                404: errors.ChemSpiPyNotFoundError,
                405: errors.ChemSpiPyMethodError,
                413: errors.ChemSpiPyPayloadError,
                429: errors.ChemSpiPyRateError,
                500: errors.ChemSpiPyServerError,
                503: errors.ChemSpiPyUnavailableError
            }.get(r.status_code, errors.ChemSpiPyHTTPError)
            raise err(message=r.reason, http_code=r.status_code)

        log.debug('Request duration: {}'.format(r.elapsed))
        return r.json()

    def get(self, api, namespace, endpoint, params=None):
        """Convenience method for making GET requests.

        :param string api: Top-level API, e.g. compounds.
        :param string namespace: API namespace, e.g. filter, lookups, records, or tools.
        :param string endpoint: Web service endpoint URL.
        :param dict params: Query parameters to add to the URL.
        :return: Web Service response JSON.
        :rtype: dict
        """
        return self.request('GET', api=api, namespace=namespace, endpoint=endpoint, params=params)

    def post(self, api, namespace, endpoint, json=None):
        """Convenience method for making POST requests.

        :param string api: Top-level API, e.g. compounds.
        :param string namespace: API namespace, e.g. filter, lookups, records, or tools.
        :param string endpoint: Web service endpoint URL.
        :param dict json: JSON data to send in the request body.
        :return: Web Service response content.
        :rtype: dict or string
        """
        return self.request('POST', api=api, namespace=namespace, endpoint=endpoint, json=json)

    def get_compound(self, csid):
        """Return a Compound object for a given ChemSpider ID.

        :param string|int csid: ChemSpider ID.
        :return: The Compound with the specified ChemSpider ID.
        :rtype: :class:`~chemspipy.objects.Compound`
        """
        return Compound(self, csid)

    def get_compounds(self, csids):
        """Return a list of Compound objects, given a list ChemSpider IDs.

        :param list[string|int] csids: List of ChemSpider IDs.
        :return: List of Compounds with the specified ChemSpider IDs.
        :rtype: list[:class:`~chemspipy.objects.Compound`]
        """
        return [Compound(self, csid) for csid in csids]

    def search(self, query, order=None, direction=ASCENDING, raise_errors=False):
        """Search ChemSpider for the specified query and return the results.

        The accepted values for ``order`` are: :data:`~chemspipy.api.RECORD_ID`, :data:`~chemspipy.api.MASS_DEFECT`,
        :data:`~chemspipy.api.MOLECULAR_WEIGHT`, :data:`~chemspipy.api.REFERENCE_COUNT`,
        :data:`~chemspipy.api.DATASOURCE_COUNT`, :data:`~chemspipy.api.PUBMED_COUNT` or
        :data:`~chemspipy.api.RSC_COUNT`.

        :param string|int query: Search query.
        :param string order: (Optional) Field to sort the result by.
        :param string direction: (Optional) :data:`~chemspipy.api.ASCENDING` or :data:`~chemspipy.api.DESCENDING`.
        :param bool raise_errors: (Optional) If True, raise exceptions. If False, store on Results ``exception``
                                  property.
        :return: Search Results list.
        :rtype: :class:`~chemspipy.search.Results`
        """
        return Results(self, self.filter_name, (query, order, direction), raise_errors=raise_errors)

    def get_datasources(self):
        """Get the list of datasources in ChemSpider.

        Many other endpoints let you restrict which sources are used to lookup the requested query. Restricting the
        sources makes queries faster.

        :return: List of datasources.
        :rtype: list[string]
        """
        response = self.get(api='compounds', namespace='lookups', endpoint='datasources')
        return response['dataSources']

    def get_details(self, record_id, fields=FIELDS):
        """Get details for a compound record.

        The available fields are listed in :data:`~chemspipy.api.FIELDS`.

        :param int record_id: Record ID.
        :param list[string] fields: (Optional) List of fields to include in the result.
        :return: Record details.
        :rtype: dict
        """
        params = {'fields': ','.join(fields)}
        endpoint = '{}/details'.format(record_id)
        response = self.get(api='compounds', namespace='records', endpoint=endpoint, params=params)
        return response

    def get_details_batch(self, record_ids, fields=FIELDS):
        """Get details for a list of compound records.

        The available fields are listed in :data:`~chemspipy.api.FIELDS`.

        :param list[int] record_ids: List of record IDs (up to 100).
        :param list[string] fields: (Optional) List of fields to include in the results.
        :return: List of record details.
        :rtype: list[dict]
        """
        json = {'recordIds': record_ids, 'fields': fields}
        response = self.post(api='compounds', namespace='records', endpoint='batch', json=json)
        return response['records']

    def get_external_references(self, record_id, datasources=None):
        """Get external references for a compound record.

        Optionally filter the results by data source. Use :meth:`~chemspipy.api.ChemSpider.get_datasources` to get the
        available datasources.

        :param int record_id: Record ID.
        :param list[string] datasources: (Optional) List of datasources to restrict the results to.
        :return: External references.
        :rtype: list[dict]
        """
        params = {}
        if datasources is not None:
            params['dataSources'] = ','.join(datasources)
        endpoint = '{}/externalreferences'.format(record_id)
        response = self.get(api='compounds', namespace='records', endpoint=endpoint, params=params)
        return response['externalReferences']

    def get_image(self, record_id):
        """Get image for a compound record.

        :param int record_id: Record ID.
        :return: Image.
        :rtype: bytes
        """
        endpoint = '{}/image'.format(record_id)
        response = self.get(api='compounds', namespace='records', endpoint=endpoint)
        return base64.b64decode(response['image'])

    def get_mol(self, record_id):
        """Get MOLfile for a compound record.

        :param int record_id: Record ID.
        :return: MOLfile.
        :rtype: string
        """
        endpoint = '{}/mol'.format(record_id)
        response = self.get(api='compounds', namespace='records', endpoint=endpoint)
        return response['sdf']

    def filter_element(self, include_elements, exclude_elements=None, include_all=False, complexity=None, isotopic=None,
                       order=None, direction=None):
        """Search compounds by element.

        Set include_all to true to only consider records that contain all of the elements in ``include_elements``,
        otherwise all records that contain any of the elements will be returned.

        A compound with a complexity of 'multiple' has more than one disconnected system in it or a metal atom or ion.

        The accepted values for ``order`` are: :data:`~chemspipy.api.RECORD_ID`, :data:`~chemspipy.api.MASS_DEFECT`,
        :data:`~chemspipy.api.MOLECULAR_WEIGHT`, :data:`~chemspipy.api.REFERENCE_COUNT`,
        :data:`~chemspipy.api.DATASOURCE_COUNT`, :data:`~chemspipy.api.PUBMED_COUNT` or
        :data:`~chemspipy.api.RSC_COUNT`.

        :param list[string] include_elements: List of up to 15 elements to search for compounds containing.
        :param list[string] exclude_elements: List of up to 100 elements to exclude compounds containing.
        :param bool include_all: (Optional) Whether to only include compounds that have all include_elements.
        :param string complexity: (Optional) 'any', 'single', or 'multiple'
        :param string isotopic: (Optional) 'any', 'labeled', or 'unlabeled'.
        :param string order: (Optional) Field to sort the result by.
        :param string direction: (Optional) :data:`~chemspipy.api.ASCENDING` or :data:`~chemspipy.api.DESCENDING`.
        :return: Query ID that may be passed to ``filter_status`` and ``filter_results``.
        :rtype: string
        """
        json = {
            'includeElements': include_elements,
            'excludeElements': exclude_elements,
            'options': {'includeAll': include_all, 'complexity': complexity, 'isotopic': isotopic},
            'orderBy': ORDERS.get(order),
            'orderDirection': DIRECTIONS.get(direction)
        }
        response = self.post(api='compounds', namespace='filter', endpoint='element', json=json)
        return response['queryId']

    def filter_formula(self, formula, datasources=None, order=None, direction=None):
        """Search compounds by formula.

        Optionally filter the results by data source. Use :meth:`~chemspipy.api.ChemSpider.get_datasources` to get the
        available datasources.

        The accepted values for ``order`` are: :data:`~chemspipy.api.RECORD_ID`, :data:`~chemspipy.api.MASS_DEFECT`,
        :data:`~chemspipy.api.MOLECULAR_WEIGHT`, :data:`~chemspipy.api.REFERENCE_COUNT`,
        :data:`~chemspipy.api.DATASOURCE_COUNT`, :data:`~chemspipy.api.PUBMED_COUNT` or
        :data:`~chemspipy.api.RSC_COUNT`.

        :param string formula: Molecular formula.
        :param list[string] datasources: (Optional) List of datasources to restrict the results to.
        :param string order: (Optional) Field to sort the result by.
        :param string direction: (Optional) :data:`~chemspipy.api.ASCENDING` or :data:`~chemspipy.api.DESCENDING`.
        :return: Query ID that may be passed to ``filter_status`` and ``filter_results``.
        :rtype: string
        """
        json = {
            'formula': formula,
            'dataSources': datasources,
            'orderBy': ORDERS.get(order),
            'orderDirection': DIRECTIONS.get(direction)
        }
        response = self.post(api='compounds', namespace='filter', endpoint='formula', json=json)
        return response['queryId']

    def filter_formula_batch(self, formulas, datasources=None, order=None, direction=None):
        """Search compounds with a list of formulas.

        Optionally filter the results by data source. Use :meth:`~chemspipy.api.ChemSpider.get_datasources` to get the
        available datasources.

        The accepted values for ``order`` are: :data:`~chemspipy.api.RECORD_ID`, :data:`~chemspipy.api.MASS_DEFECT`,
        :data:`~chemspipy.api.MOLECULAR_WEIGHT`, :data:`~chemspipy.api.REFERENCE_COUNT`,
        :data:`~chemspipy.api.DATASOURCE_COUNT`, :data:`~chemspipy.api.PUBMED_COUNT` or
        :data:`~chemspipy.api.RSC_COUNT`.

        :param list[string] formulas: Molecular formula.
        :param list[string] datasources: (Optional) List of datasources to restrict the results to.
        :param string order: (Optional) Field to sort the result by.
        :param string direction: (Optional) :data:`~chemspipy.api.ASCENDING` or :data:`~chemspipy.api.DESCENDING`.
        :return: Query ID that may be passed to ``filter_formula_batch_status`` and ``filter_formula_batch_results``.
        :rtype: string
        """
        json = {
            'formulas': formulas,
            'dataSources': datasources,
            'orderBy': ORDERS.get(order),
            'orderDirection': DIRECTIONS.get(direction)
        }
        response = self.post(api='compounds', namespace='filter', endpoint='formula/batch', json=json)
        return response['queryId']

    def filter_formula_batch_status(self, query_id):
        """Get formula batch filter status using a query ID that was returned by a previous filter request.

        :param string query_id: Query ID from a previous formula batch filter request.
        :return: Status dict with 'status', 'count', and 'message' fields.
        :rtype: dict
        """
        endpoint = 'formula/batch/{}/status'.format(query_id)
        response = self.get(api='compounds', namespace='filter', endpoint=endpoint)
        return response

    def filter_formula_batch_results(self, query_id):
        """Get formula batch filter results using a query ID that was returned by a previous filter request.

        Each result is a dict containing a ``formula`` key and a ``results`` key.

        :param string query_id: Query ID from a previous formula batch filter request.
        :return: List of results.
        :rtype: list[dict]
        """
        endpoint = 'formula/batch/{}/results'.format(query_id)
        response = self.get(api='compounds', namespace='filter', endpoint=endpoint)
        return response['batchResults']

    def filter_inchi(self, inchi):
        """Search compounds by InChI.

        :param string inchi: InChI.
        :return: Query ID that may be passed to ``filter_status`` and ``filter_results``.
        :rtype: string
        """
        json = {'inchi': inchi}
        response = self.post(api='compounds', namespace='filter', endpoint='inchi', json=json)
        return response['queryId']

    def filter_inchikey(self, inchikey):
        """Search compounds by InChIKey.

        :param string inchikey: InChIKey.
        :return: Query ID that may be passed to ``filter_status`` and ``filter_results``.
        :rtype: string
        """
        json = {'inchikey': inchikey}
        response = self.post(api='compounds', namespace='filter', endpoint='inchikey', json=json)
        return response['queryId']

    def filter_intrinsicproperty(self, formula=None, molecular_weight=None, nominal_mass=None, average_mass=None,
                                 monoisotopic_mass=None, molecular_weight_range=None, nominal_mass_range=None,
                                 average_mass_range=None, monoisotopic_mass_range=None, complexity=None, isotopic=None,
                                 order=None, direction=None):
        """Search compounds by intrinsic property, such as formula and mass.

        At least one of formula, molecular_weight, nominal_mass, average_mass, monoisotopic_mass must be specified.

        A compound with a complexity of 'multiple' has more than one disconnected system in it or a metal atom or ion.

        The accepted values for ``order`` are: :data:`~chemspipy.api.RECORD_ID`, :data:`~chemspipy.api.MASS_DEFECT`,
        :data:`~chemspipy.api.MOLECULAR_WEIGHT`, :data:`~chemspipy.api.REFERENCE_COUNT`,
        :data:`~chemspipy.api.DATASOURCE_COUNT`, :data:`~chemspipy.api.PUBMED_COUNT` or
        :data:`~chemspipy.api.RSC_COUNT`.

        :param string formula: Molecular formula.
        :param float molecular_weight: Molecular weight.
        :param float nominal_mass: Nominal mass.
        :param float average_mass: Average mass.
        :param float monoisotopic_mass: Monoisotopic mass.
        :param float molecular_weight_range: Molecular weight range.
        :param float nominal_mass_range: Nominal mass range.
        :param float average_mass_range: Average mass range.
        :param float monoisotopic_mass_range: Monoisotopic mass range.
        :param string complexity: (Optional) 'any', 'single', or 'multiple'
        :param string isotopic: (Optional) 'any', 'labeled', or 'unlabeled'.
        :param string order: (Optional) Field to sort the result by.
        :param string direction: (Optional) :data:`~chemspipy.api.ASCENDING` or :data:`~chemspipy.api.DESCENDING`.
        :return: Query ID that may be passed to ``filter_status`` and ``filter_results``.
        :rtype: string
        """
        json = {
            'formula': formula,
            'options': {'complexity': complexity, 'isotopic': isotopic},
            'orderBy': ORDERS.get(order),
            'orderDirection': DIRECTIONS.get(direction)
        }
        if molecular_weight is not None and molecular_weight_range is not None:
            json['molecularWeight'] = {'mass': molecular_weight, 'range': molecular_weight_range}
        if nominal_mass is not None and nominal_mass_range is not None:
            json['nominalMass'] = {'mass': nominal_mass, 'range': nominal_mass_range}
        if average_mass is not None and average_mass_range is not None:
            json['averageMass'] = {'mass': average_mass, 'range': average_mass_range}
        if monoisotopic_mass is not None and monoisotopic_mass_range is not None:
            json['monoisotopicMass'] = {'mass': monoisotopic_mass, 'range': monoisotopic_mass_range}
        response = self.post(api='compounds', namespace='filter', endpoint='intrinsicproperty', json=json)
        return response['queryId']

    def filter_mass(self, mass, mass_range, datasources=None, order=None, direction=None):
        """Search compounds by mass.

        Filter to compounds within ``mass_range`` of the given ``mass``.

        Optionally filter the results by data source. Use :meth:`~chemspipy.api.ChemSpider.get_datasources` to get the
        available datasources.

        The accepted values for ``order`` are: :data:`~chemspipy.api.RECORD_ID`, :data:`~chemspipy.api.MASS_DEFECT`,
        :data:`~chemspipy.api.MOLECULAR_WEIGHT`, :data:`~chemspipy.api.REFERENCE_COUNT`,
        :data:`~chemspipy.api.DATASOURCE_COUNT`, :data:`~chemspipy.api.PUBMED_COUNT` or
        :data:`~chemspipy.api.RSC_COUNT`.

        :param float mass: Mass between 1 and 11000 Atomic Mass Units.
        :param float mass_range: Mass range between 0.0001 and 100 Atomic Mass Units.
        :param list[string] datasources: (Optional) List of datasources to restrict the results to.
        :param string order: (Optional) Field to sort the result by.
        :param string direction: (Optional) :data:`~chemspipy.api.ASCENDING` or :data:`~chemspipy.api.DESCENDING`.
        :return: Query ID that may be passed to ``filter_status`` and ``filter_results``.
        :rtype: string
        """
        json = {
            'mass': mass,
            'range': mass_range,
            'dataSources': datasources,
            'orderBy': ORDERS.get(order),
            'orderDirection': DIRECTIONS.get(direction)
        }
        response = self.post(api='compounds', namespace='filter', endpoint='mass', json=json)
        return response['queryId']

    def filter_mass_batch(self, masses, datasources=None, order=None, direction=None):
        """Search compounds with a list of masses and mass ranges.

        The ``masses`` parameter should be a list of tuples, each with two elements: A mass, and a mass range::

            qid = cs.filter_mass_batch(masses=[(12, 0.001), (24, 0.001)])

        Optionally filter the results by data source. Use :meth:`~chemspipy.api.ChemSpider.get_datasources` to get the
        available datasources.

        The accepted values for ``order`` are: :data:`~chemspipy.api.RECORD_ID`, :data:`~chemspipy.api.MASS_DEFECT`,
        :data:`~chemspipy.api.MOLECULAR_WEIGHT`, :data:`~chemspipy.api.REFERENCE_COUNT`,
        :data:`~chemspipy.api.DATASOURCE_COUNT`, :data:`~chemspipy.api.PUBMED_COUNT` or
        :data:`~chemspipy.api.RSC_COUNT`.

        :param list[tuple[float, float]] masses: List of (mass, range) tuples.
        :param list[string] datasources: (Optional) List of datasources to restrict the results to.
        :param string order: (Optional) Field to sort the result by.
        :param string direction: (Optional) :data:`~chemspipy.api.ASCENDING` or :data:`~chemspipy.api.DESCENDING`.
        :return: Query ID that may be passed to ``filter_formula_batch_status`` and ``filter_formula_batch_results``.
        :rtype: string
        """
        masses = [{'mass': m, 'range': r} for m, r in masses]
        json = {
            'masses': masses,
            'dataSources': datasources,
            'orderBy': ORDERS.get(order),
            'orderDirection': DIRECTIONS.get(direction)
        }
        response = self.post(api='compounds', namespace='filter', endpoint='mass/batch', json=json)
        return response['queryId']

    def filter_mass_batch_status(self, query_id):
        """Get formula batch filter status using a query ID that was returned by a previous filter request.

        :param string query_id: Query ID from a previous formula batch filter request.
        :return: Status dict with 'status', 'count', and 'message' fields.
        :rtype: dict
        """
        endpoint = 'mass/batch/{}/status'.format(query_id)
        response = self.get(api='compounds', namespace='filter', endpoint=endpoint)
        return response

    def filter_mass_batch_results(self, query_id):
        """Get formula batch filter results using a query ID that was returned by a previous filter request.

        Each result is a dict containing a ``formula`` key and a ``results`` key.

        :param string query_id: Query ID from a previous formula batch filter request.
        :return: List of results.
        :rtype: list[dict]
        """
        endpoint = 'mass/batch/{}/results'.format(query_id)
        response = self.get(api='compounds', namespace='filter', endpoint=endpoint)
        return response['batchResults']

    def filter_name(self, name, order=None, direction=None):
        """Search compounds by name.

        The accepted values for ``order`` are: :data:`~chemspipy.api.RECORD_ID`, :data:`~chemspipy.api.MASS_DEFECT`,
        :data:`~chemspipy.api.MOLECULAR_WEIGHT`, :data:`~chemspipy.api.REFERENCE_COUNT`,
        :data:`~chemspipy.api.DATASOURCE_COUNT`, :data:`~chemspipy.api.PUBMED_COUNT` or
        :data:`~chemspipy.api.RSC_COUNT`.

        :param string name: Compound name.
        :param string order: (Optional) Field to sort the result by.
        :param string direction: (Optional) :data:`~chemspipy.api.ASCENDING` or :data:`~chemspipy.api.DESCENDING`.
        :return: Query ID that may be passed to ``filter_status`` and ``filter_results``.
        :rtype: string
        """
        json = {'name': name, 'orderBy': ORDERS.get(order),  'orderDirection': DIRECTIONS.get(direction)}
        response = self.post(api='compounds', namespace='filter', endpoint='name', json=json)
        return response['queryId']

    def filter_smiles(self, smiles):
        """Search compounds by SMILES.

        :param string smiles: Compound SMILES.
        :return: Query ID that may be passed to ``filter_status`` and ``filter_results``.
        :rtype: string
        """
        json = {'smiles': smiles}
        response = self.post(api='compounds', namespace='filter', endpoint='smiles', json=json)
        return response['queryId']

    def filter_status(self, query_id):
        """Get filter status using a query ID that was returned by a previous filter request.

        :param string query_id: Query ID from a previous filter request.
        :return: Status dict with 'status', 'count', and 'message' fields.
        :rtype: dict
        """
        endpoint = '{}/status'.format(query_id)
        response = self.get(api='compounds', namespace='filter', endpoint=endpoint)
        return response

    def filter_results(self, query_id, start=None, count=None):
        """Get filter results using a query ID that was returned by a previous filter request.

        :param string query_id: Query ID from a previous filter request.
        :param int start: Zero-based results offset.
        :param int count: Number of results to return.
        :return: List of results.
        :rtype: list[int]
        """
        endpoint = '{}/results'.format(query_id)
        params = {'start': start, 'count': count}
        response = self.get(api='compounds', namespace='filter', endpoint=endpoint, params=params)
        return response['results']

    def filter_results_sdf(self, query_id):
        """Get filter results as SDF file using a query ID that was returned by a previous filter request.

        :param string query_id: Query ID from a previous filter request.
        :return: SDF file containing the results.
        :rtype: bytes
        """
        endpoint = '{}/results/sdf'.format(query_id)
        response = self.get(api='compounds', namespace='filter', endpoint=endpoint)
        return zlib.decompress(base64.b64decode(response['results']), 16 + zlib.MAX_WBITS)

    def convert(self, input, input_format, output_format):
        """Convert a chemical from one format to another.

        Format: ``SMILES``, ``InChI``, ``InChIKey`` or ``Mol``.

        Allowed conversions: from InChI to InChIKey, from InChI to Mol file, from InChI to SMILES, from InChIKey to
        InChI, from InChIKey to Mol file, from Mol file to InChI, from Mol file to InChIKey, from SMILES to InChI.

        :param string input: Input chemical.
        :param string input_format: Input format.
        :param string output_format: Output format.
        :return: Input chemical in output format.
        :rtype: string
        """
        json = {'input': input, 'inputFormat': input_format, 'outputFormat': output_format}
        response = self.post(api='compounds', namespace='tools', endpoint='convert', json=json)
        return response['output']

    def validate_inchikey(self, inchikey):
        """Return whether ``inchikey`` is valid.

        :param string inchikey: The InChIKey to validate.
        :return: Whether the InChIKey is valid.
        :rtype: bool
        """
        json = {'inchikey': inchikey}
        try:
            response = self.post(api='compounds', namespace='tools', endpoint='validate/inchikey', json=json)
            return response['valid']
        except errors.ChemSpiPyHTTPError:
            return False

    def get_databases(self):
        """Get the list of datasources in ChemSpider.

        .. deprecated:: 2.0.0
           Use :py:meth:`~chemspipy.api.ChemSpider.get_datasources` instead.

        """
        warnings.warn('Use get_datasources instead of get_databases.', DeprecationWarning)
        return self.get_datasources()

    def get_extended_compound_info(self, csid):
        """Get extended record details for a CSID.

        .. deprecated:: 2.0.0
           Use :py:meth:`~chemspipy.api.ChemSpider.get_details` instead.

        :param string|int csid: ChemSpider ID.
        """
        warnings.warn('Use get_details instead of get_extended_compound_info.', DeprecationWarning)
        return self.get_details(record_id=csid)

    def get_extended_compound_info_list(self, csids):
        """Get extended record details for a list of CSIDs.

        .. deprecated:: 2.0.0
           Use :py:meth:`~chemspipy.api.ChemSpider.get_details_batch` instead.

        :param list[string|int] csids: ChemSpider IDs.
        """
        warnings.warn('Use get_details_batch instead of get_extended_compound_info.', DeprecationWarning)
        return self.get_details_batch(record_ids=csids)

    def get_extended_mol_compound_info_list(self, csids, mol_type=MOL2D, include_reference_counts=False,
                                            include_external_references=False):
        """Get extended record details (including MOL) for a list of CSIDs.

        A maximum of 250 CSIDs can be fetched per request.

        .. deprecated:: 2.0.0
           Use :py:meth:`~chemspipy.api.ChemSpider.get_details_batch` instead.

        :param list[string|int] csids: ChemSpider IDs.
        :param string mol_type: :data:`~chemspipy.api.MOL2D`, :data:`~chemspipy.api.MOL3D` or
                                :data:`~chemspipy.api.BOTH`.
        :param bool include_reference_counts: Whether to include reference counts.
        :param bool include_external_references: Whether to include external references.
        """
        warnings.warn('Use get_details_batch instead of get_extended_mol_compound_info_list.', DeprecationWarning)
        return self.get_details_batch(record_ids=csids)

    def get_record_mol(self, csid, calc3d=False):
        """Get ChemSpider record in MOL format.

        .. deprecated:: 2.0.0
           Use :py:meth:`~chemspipy.api.ChemSpider.get_mol` instead.

        :param string|int csid: ChemSpider ID.
        :param bool calc3d: Whether 3D coordinates should be calculated before returning record data.
        """
        warnings.warn('Use get_mol instead of get_record_mol.', DeprecationWarning)
        if calc3d:
            warnings.warn('calc3d parameter for get_record_mol is no longer supported.', DeprecationWarning)
        return self.get_mol(record_id=csid)

    def async_simple_search(self, query):
        """Search ChemSpider with arbitrary query, returning results in order of the best match found.

        This method returns a transaction ID which can be used with other methods to get search status and results.

        .. deprecated:: 2.0.0
           Use :py:meth:`~chemspipy.api.ChemSpider.filter_name` instead.

        :param string query: Search query - a name, SMILES, InChI, InChIKey, CSID, etc.
        :return: Transaction ID.
        :rtype: string
        """
        warnings.warn('Use filter_name instead of async_simple_search.', DeprecationWarning)
        return self.filter_name(name=query)

    def async_simple_search_ordered(self, query, order=CSID, direction=ASCENDING):
        """Search ChemSpider with arbitrary query, returning results with a custom order.

        This method returns a transaction ID which can be used with other methods to get search status and results.

        .. deprecated:: 2.0.0
           Use :meth:`~chemspipy.api.ChemSpider.filter_name` instead.

        :param string query: Search query - a name, SMILES, InChI, InChIKey, CSID, etc.
        :param string order: (Optional) Field to sort the result by.
        :param string direction: (Optional) :data:`~chemspipy.api.ASCENDING` or :data:`~chemspipy.api.DESCENDING`.
        :return: Transaction ID.
        :rtype: string
        """
        warnings.warn('Use filter_name instead of async_simple_search.', DeprecationWarning)
        return self.filter_name(name=query, order=order, direction=direction)

    def get_async_search_status(self, rid):
        """Check the status of an asynchronous search operation.

        .. deprecated:: 2.0.0
           Use :meth:`~chemspipy.api.ChemSpider.filter_status` instead.

        :param string rid: A transaction ID, returned by an asynchronous search method.
        :return: Unknown, Created, Scheduled, Processing, Suspended, PartialResultReady, ResultReady, Failed,
                  TooManyRecords
        :rtype: string
        """
        warnings.warn('Use filter_status instead of get_async_search_status.', DeprecationWarning)
        return self.filter_status(query_id=rid)['status']

    def get_async_search_status_and_count(self, rid):
        """Check the status of an asynchronous search operation. If ready, a count and message are also returned.

        .. deprecated:: 2.0.0
           Use :meth:`~chemspipy.api.ChemSpider.filter_status` instead.

        :param string rid: A transaction ID, returned by an asynchronous search method.
        :rtype: dict
        """
        warnings.warn('Use filter_status instead of get_async_search_status_and_count.', DeprecationWarning)
        return self.filter_status(query_id=rid)

    def get_async_search_result(self, rid):
        """Get the results from a asynchronous search operation.

        .. deprecated:: 2.0.0
           Use :meth:`~chemspipy.api.ChemSpider.filter_results` instead.

        :param string rid: A transaction ID, returned by an asynchronous search method.
        :return: A list of Compounds.
        :rtype: list[:class:`~chemspipy.objects.Compound`]
        """
        warnings.warn('Use filter_results instead of get_async_search_result.', DeprecationWarning)
        results = self.filter_results(query_id=rid)
        return [Compound(self, record_id) for record_id in results]

    def get_async_search_result_part(self, rid, start=0, count=-1):
        """Get a slice of the results from a asynchronous search operation.

        .. deprecated:: 2.0.0
           Use :meth:`~chemspipy.api.ChemSpider.filter_results` instead.

        :param string rid: A transaction ID, returned by an asynchronous search method.
        :param int start: The number of results to skip.
        :param int count: The number of results to return. -1 returns all through to end.
        :return: A list of Compounds.
        :rtype: list[:class:`~chemspipy.objects.Compound`]
        """
        warnings.warn('Use filter_results instead of get_async_search_result_part.', DeprecationWarning)
        if count == -1:
            count = None
        results = self.filter_results(query_id=rid, start=start, count=count)
        return [Compound(self, record_id) for record_id in results]

    def get_compound_info(self, csid):
        """Get SMILES, StdInChI and StdInChIKey for a given CSID.

        .. deprecated:: 2.0.0
           Use :meth:`~chemspipy.api.ChemSpider.get_details` instead.

        :param string|int csid: ChemSpider ID.
        :rtype: dict
        """
        warnings.warn('Use get_details instead of get_compound_info.', DeprecationWarning)
        return self.get_details(record_id=csid)

    def get_compound_thumbnail(self, csid):
        """Get PNG image as binary data.

        .. deprecated:: 2.0.0
           Use :meth:`~chemspipy.api.ChemSpider.get_image` instead.

        :param string|int csid: ChemSpider ID.
        :rtype: bytes
        """
        warnings.warn('Use get_image instead of get_compound_thumbnail.', DeprecationWarning)
        return self.get_image(record_id=csid)

    def simple_search(self, query):
        """Search ChemSpider with arbitrary query.

        .. deprecated:: 2.0.0
           Use :meth:`~chemspipy.api.ChemSpider.search` instead.

        :param string query: Search query - a chemical name.
        :return: Search Results list.
        :rtype: :class:`~chemspipy.search.Results`
        """
        warnings.warn('Use search instead of simple_search.', DeprecationWarning)
        return self.search(query=query)
