#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 19:48:16 2022

@author: olivier.vitrac@agroparistech.fr
"""

# caller layer
from patankar.layer import layer
A=layer(D=1e-14,l=50e-6)
A

# caller senspatankar
from patankar.senspatankar import sens
S=sens(Bi=5,l=1,myid="mine")
print(S)
S.struct()
