# -*- coding: utf-8 -*-
"""
chemspipy.search
~~~~~~~~~~~~~~~~

A wrapper for asynchronous search requests.

"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import datetime
import logging
import threading
import time

from six.moves import range

from . import errors, objects, utils


log = logging.getLogger(__name__)


# TODO: Use Sequence abc metaclass?
class Results(object):
    """Container class to perform a search on a background thread and hold the results when ready."""

    def __init__(self, cs, searchfunc, searchargs, raise_errors=False, max_requests=40):
        """Generally shouldn't be instantiated directly. See :meth:`~chemspipy.api.ChemSpider.search` instead.

        :param ChemSpider cs: ``ChemSpider`` session.
        :param function searchfunc: Search function that returns a transaction ID.
        :param tuple searchargs: Arguments for the search function.
        :param bool raise_errors: If True, raise exceptions. If False, store on ``exception`` property.
        :param int max_requests: Maximum number of times to check if search results are ready.
        """
        log.debug('Results init')
        self._cs = cs
        self._raise_errors = raise_errors
        self._max_requests = max_requests
        self._status = 'Created'
        self._exception = None
        self._qid = None
        self._message = None
        self._start = None
        self._end = None
        self._results = []
        self._searchthread = threading.Thread(name='SearchThread', target=self._search, args=(cs, searchfunc, searchargs))
        self._searchthread.start()

    def _search(self, cs, searchfunc, searchargs):
        """Perform the search and retrieve the results."""
        log.debug('Searching in background thread')
        self._start = datetime.datetime.utcnow()
        try:
            self._qid = searchfunc(*searchargs)
            log.debug('Setting qid: %s' % self._qid)
            for _ in range(self._max_requests):
                log.debug('Checking status: %s' % self._qid)
                status = cs.filter_status(self._qid)
                self._status = status['status']
                self._message = status.get('message', '')
                log.debug(status)
                time.sleep(0.2)
                if status['status'] == 'Complete':
                    break
                elif status['status'] in {'Failed', 'Unknown', 'Suspended', 'Not Found'}:
                    raise errors.ChemSpiPyServerError('Search Failed: %s' % status.get('message', ''))
            else:
                raise errors.ChemSpiPyTimeoutError('Search took too long')
            log.debug('Search success!')
            self._end = datetime.datetime.utcnow()
            if status['count'] > 0:
                self._results = [objects.Compound(cs, csid) for csid in cs.filter_results(self._qid)]
                log.debug('Results: %s', self._results)
            elif not self._message:
                self._message = 'No results found'
        except Exception as e:
            # Catch and store exception so we can raise it in the main thread
            self._exception = e
            self._end = datetime.datetime.utcnow()
            if self._status == 'Created':
                self._status = 'Failed'

    def ready(self):
        """Return True if the search finished.

        :rtype: bool
        """
        return not self._searchthread.is_alive()

    def success(self):
        """Return True if the search finished with no errors.

        :rtype: bool
        """
        return self.ready() and not self._exception

    def wait(self):
        """Block until the search has completed and optionally raise any resulting exception."""
        log.debug('Waiting for search to finish')
        self._searchthread.join()
        if self._exception and self._raise_errors:
            raise self._exception

    @property
    def status(self):
        """Current status string returned by ChemSpider.

        :return: 'Unknown', 'Created', 'Scheduled', 'Processing', 'Suspended', 'PartialResultReady', 'ResultReady'
        :rtype: string
        """
        return self._status

    @property
    def exception(self):
        """Any Exception raised during the search. Blocks until the search is finished."""
        self.wait()  # TODO: If raise_errors=True this will raise the exception when trying to access it?
        return self._exception

    @property
    def qid(self):
        """Search query ID.

        :rtype: string
        """
        return self._qid

    @property
    def message(self):
        """A contextual message about the search. Blocks until the search is finished.

        :rtype: string
        """
        self.wait()
        return self._message

    @property
    def count(self):
        """The number of search results. Blocks until the search is finished.

        :rtype: int
        """
        return len(self)

    @property
    def duration(self):
        """The time taken to perform the search. Blocks until the search is finished.

        :rtype: :py:class:`datetime.timedelta`
        """
        self.wait()
        return self._end - self._start

    @utils.memoized_property
    def sdf(self):
        """Get an SDF containing all the search results.

        :return: SDF containing the search results.
        :rtype: bytes
        """
        self.wait()
        return self._cs.filter_results_sdf(self._qid)

    def __getitem__(self, index):
        """Get a single result or a slice of results. Blocks until the search is finished.

        This means a Results instance can be treated like a normal Python list. For example::

            cs.search('glucose')[2]
            cs.search('glucose')[0:2]

        An IndexError will be raised if the index is greater than the total number of results.
        """
        self.wait()
        return self._results.__getitem__(index)

    def __len__(self):
        self.wait()
        return self._results.__len__()

    def __iter__(self):
        self.wait()
        return iter(self._results)

    def __repr__(self):
        if self.success():
            return 'Results(%s)' % self._results
        else:
            return 'Results(%s)' % self.status
