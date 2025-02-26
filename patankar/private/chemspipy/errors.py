# -*- coding: utf-8 -*-
"""
chemspipy.errors
~~~~~~~~~~~~~~~~

Exceptions raised by ChemSpiPy.

"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division


class ChemSpiPyError(Exception):
    """Root ChemSpiPy Exception."""
    pass


class ChemSpiPyHTTPError(ChemSpiPyError):
    """Base exception to handle HTTP errors."""

    #: Default message if none supplied. Override in subclasses.
    MESSAGE = 'ChemSpiPy Error'
    HTTP_CODE = None

    def __init__(self, message=None, http_code=None, *args, **kwargs):
        """

        :param string|bytes message: Error message.
        :param http_code: HTTP code.
        """

        # Decode message to unicode if necessary
        if isinstance(message, bytes):
            try:
                message = message.decode('utf-8')
            except UnicodeDecodeError:
                message = message.decode('iso-8859-1')

        self.message = message if message is not None else self.MESSAGE
        self.http_code = http_code if http_code is not None else self.HTTP_CODE
        super(ChemSpiPyHTTPError, self).__init__(*args, **kwargs)

    def __repr__(self):
        args = 'message={!r}'.format(self.message)
        if self.http_code is not None:
            args += ', http_code={!r}'.format(self.http_code)
        return '{}({})'.format(self.__class__.__name__, args)

    def __str__(self):
        return self.message


class ChemSpiPyBadRequestError(ChemSpiPyHTTPError):
    """Raised for a bad request."""
    MESSAGE = 'Bad request.'
    HTTP_CODE = 400


class ChemSpiPyAuthError(ChemSpiPyHTTPError):
    """Raised when API key authorization fails."""
    MESSAGE = 'Unauthorized.'
    HTTP_CODE = 401


class ChemSpiPyNotFoundError(ChemSpiPyHTTPError):
    """Raised when the requested resource was not found."""
    MESSAGE = 'Not found.'
    HTTP_CODE = 404


class ChemSpiPyMethodError(ChemSpiPyHTTPError):
    """Raised when an invalid HTTP method is used."""
    MESSAGE = 'Method Not Allowed.'
    HTTP_CODE = 405


class ChemSpiPyPayloadError(ChemSpiPyHTTPError):
    """Raised when a request payload is too large."""
    MESSAGE = 'Payload Too Large.'
    HTTP_CODE = 413


class ChemSpiPyRateError(ChemSpiPyHTTPError):
    """Raised when too many requests are sent in a given amount of time."""
    MESSAGE = 'Too Many Requests.'
    HTTP_CODE = 429


class ChemSpiPyServerError(ChemSpiPyHTTPError):
    """Raised when an internal server error occurs."""
    MESSAGE = 'Internal Server Error.'
    HTTP_CODE = 500


class ChemSpiPyUnavailableError(ChemSpiPyHTTPError):
    """Raised when the service is temporarily unavailable."""
    MESSAGE = 'Service Unavailable.'
    HTTP_CODE = 503


class ChemSpiPyTimeoutError(ChemSpiPyError):
    """Raised when an asynchronous request times out."""
    pass
