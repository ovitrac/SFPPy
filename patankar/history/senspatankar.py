#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SENSPATANKAR forks of SENSPATANKAR on Windows


@author: olivier.vitrac@agroparistech.fr
"""

# -- History ----
# Created on Mon Jan 17 19:20:53 2022
# RC not yet

oneline = "Simulate mass transfer from a n-multilayer to food"

docstr = """
Packaging layers are numbered 1..n
Food is numbered 0
Layer 1 is in contact with food.

Example of caller:
    from patankar.senspatankar import sens
    S=sens(nlayer=2,Bi=20)
    S

"""

# History
# 2022-01-17 first coding


# Dependencies
#import types
#from math import *
from uuid import uuid1 as autoid
import numpy as np

# Class definition
class sens:
    """ mass transfer solver from patankar package (formerly senspantakar)"""
    # --------------------------------------------------------------------
    # PRIVATE PROPERTIES (cannot be changed by the user if they start by __, _=protected)
    # --------------------------------------------------------------------
    __description = "SENSPATANKAR object"
    __version = 0.001             # version of sens
    __contact = "olivier.vitrac@agroparistech.fr" # contact person
    __sysid = str(autoid())       # internal id to check inheritances
    
    _printformat = '{0:.4g}';    # format to display D, k, l values
    CF0 = 0                      # initial concentration in food (static)

    # --------------------------------------------------------------------
    # CONSTRUCTOR OF INSTANCE PROPERTIES
    # --------------------------------------------------------------------
    def __init__(self,          # constructor (note that self is the object itself)
                 myid=__sysid,   # default if name
                 nlayer=3,      # number of layers (n)
                 k0=1,          # ref Henry like for food (a.u.)
                 k=1,           # ref Henry like for layers 1..n (a.u.)
                 l=50e-6,       # default layer thichness
                 D=1e-14,       # ref diffussivity for layers 1..n (m^2/s)
                 C0=1000,       # ref concentration (a.u.)
                 Bi=1e3,        # default mass Biot number (relative to l1)
                 L=200/1800,    # default dilution ratio (l1+l2..+ln)/l0
                 ):
        """ 

        Implemented properties
        ----------------------
         constructor
        myid  : TYPE, optional
                DESCRIPTION. User id
        nlayer : TYPE, optional
                 DESCRIPTION. number of layers. The default is 3.
        k0 : TYPE, optional
            DESCRIPTION. ref Henry like for food (a.u.). The default is 1.
        k : TYPE, optional
            DESCRIPTION. ref Henry like for layers 1..n (a.u.). The default is 1.
        l : TYPE, optional
            DESCRIPTION. default layer thichness. The default is 50e-6.
        D : TYPE, optional
            DESCRIPTION. ref diffussivity for layers 1..n (m^2/s). The default is 1e-14.
        C0 : TYPE, optional
            DESCRIPTION. ref concentration (a.u.). The default is 1000.
        Bi : TYPE, optional
            DESCRIPTION. default mass Biot number (relative to l1). The default is 1e3.
        L : TYPE, optional
            DESCRIPTION. default dilution ratio (l1+l2..+ln)/l0. The default is 200/1800.

        Returns
        -------
        None.

        """
        # object id
        self.myid = myid
        # scalar definitions
        self.nlayer = nlayer            # number of layers
        self.Bi = Bi                    # mass Biot number
        self.k0 = k0                    # liquid = 0 (reference value = 1)
        self.L  = L                     # dilution factor (respectively to iref) = Vf/Vp
        # vector definitions
        ones = np.ones(nlayer)
        self.k  = ones*k     # Henry like constants for each layer [l1,l2,l3,l4] with layer1=layer in contact with food
        self.D  = ones*D     # diffusion coeff in each layer (m^2/s)
        self.l  = ones*l     # in m, equivalent length l=V/Sref
        self.C0 = ones*C0    # concentration (a.u.)
        
    # --------------------------------------------------------------------
    # DISP method
    # --------------------------------------------------------------------
    def __repr__(self):
        """ disp method """
        print("[%s version=%0.4g, contact=%s]" % (self.__description,self.__version,self.__contact))
        print('l = %s' % (np.array2string(self.l, formatter={'float_kind':self._printformat.format})) )
        print('D = %s' % (np.array2string(self.D, formatter={'float_kind':self._printformat.format})) )
        print('k = %s' % (np.array2string(self.k, formatter={'float_kind':self._printformat.format})) )
        print('C0 = %s' % (np.array2string(self.C0, formatter={'float_kind':self._printformat.format})) )
        print('Bi = %0.4g, k0 = %0.4g, CF0 = %0.4g' % (self.Bi,self.k0,self.CF0))
        ret=('nlayer=%d %s with id="%s"' % (self.nlayer,self.__description,self.myid))
        return ret
    
    # --------------------------------------------------------------------
    # STRUCT method - returns the equivalent dictionary from an object
    # source: https://stackoverflow.com/questions/61517/python-dictionary-from-an-objects-fields
    # --------------------------------------------------------------------        
    def struct(self):
        """             returns the equivalent dictionary from an object """
        return dict((key, getattr(self, key)) for key in dir(self) if key not in dir(self.__class__))
  # --------------------------------------------------------------------