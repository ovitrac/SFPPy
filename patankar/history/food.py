#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Foodlayer builder for patankar package
# general object is foodlayer()
# derived ones yogurt(),ethanol(),ethanol50()...water()

Created on Tue Feb 15 10:26:44 2022

@author: INRAE\olivier.vitrac@agroparistech.fr

"""

# Revision history
# 2022-02-15 alpha version, check attributes in foodlayer.__repr__

# %% ALL ROOT CLASSES

# Top parent class FOODLAYER (for all food types and food simulants)
class foodlayer():
    """ parent foodlayer class """
    description = "general food class"
    name = "unknown"
    volume = 1e-3
    surfacearea = 0.06
    
    # display/representation method
    def __repr__(self):
        print('Food object "%s" (%s) with properties:' % (self.name,self.description))
        print("%15s: %0.8g" % ("volume",self.volume))
        print("%15s: %0.8g" % ("surface area",self.surfacearea))
        if hasattr(self,'k0'): print("%15s: %0.8g" % ("k0",self.k0))
        if hasattr(self,'h'): print("%15s: %0.8g" % ("h",self.h))
        
        return "%s (%s)" % (self.name,self.description)

# Top parent class TEXTURE
class texture():
    """ parent food texture class """
    description = "default class texture"
    name = "undefined"
    h = 1e-3

# Top parent class chemicalaffinity
class chemicalaffinity():
    """ parent chemical affinity class """
    description = "default chemical affinity"
    name = "undefined"
    k0 = 1    # p = k0*C with k0 Henry-like coefficient, p=partial pressure, C=concentration
    
# %% SECOND LEVEL CLASSES

# Second level of foodlayer class: real food vs simulant
class realfood(foodlayer):
    """ core real food class (second level)"""
    description = "real food class"
    
class simulant(foodlayer):
    """ core food simulant class (second level)"""
    name = "generic food simulant"
    description = "food simulant"

# Second level of class texture
class solid(texture):
    """ solid food texture """
    name = "solid food"
    description = "solid food products"
    h  = 1e-8 # m/s (convective mass transfer coefficient)

class semisolid(texture):
    """ semi-solid food texture """
    name = "solid food"
    description = "solid food products"
    h  = 1e-7 # m/s (convective mass transfer coefficient)

class liquid(texture):
    """ liquid food texture """
    name = "liquid food"
    description = "liquid food products"
    h  = 1e-6 # m/s (convective mass transfer coefficient)
    
class perfectlymixed(texture):
    """ perfectly mixed liquid (texture) """
    name = "perfectly mixed liquid"
    description = "maximize mixing, minimize the mass transfer boundary layer" 
    h  = 1e-4 # m/s (convective mass transfer coefficient)
    
# Second level of chemical affinity
class fat(chemicalaffinity):
    """ fat contact """
    name = "fat contact"
    description = "maximize mass transfer"
    k0 = 1
    
class aqueous(chemicalaffinity):
    """ aqueous food contact """
    name = "aqueous contact"
    description = "minimize mass transfer"
    k0 = 1000
    
class intermediate(chemicalaffinity):
    """ intermediate chemical affinity """
    name = "intermediate"
    description = "intermediate chemical affinity"
    k0 = 10

# %% THIRD LEVEL CLASSES

# Third level: list of food simulants
class ethanol(simulant,perfectlymixed,fat):
    """ ethanol food simulant """
    name = "ethanol"
    description = "ethanol = from pure ethanol down to ethanol 95%"
    
class ethanol50(simulant,perfectlymixed,intermediate):
    """ ethanol 50 food simulant """
    name = "ethanol 50"
    description = "ethanol 50, food simulant of dairy products"
    
class water(simulant,perfectlymixed,aqueous):
    """ water food simulant """
    name = "water"
    description = "water food simulant"

class tenax(simulant,solid,fat):
    """ Tenax(r) food simulant """
    name = "Tenax"
    description = "simulant of dry food products"
    

# Third level: list of food (first classes = higher precedence)
class yogurt(realfood,semisolid,ethanol50):
    """ yogurt as an example of real food """
    # static properties
    description = "yogurt"
    volume = 0.125e-3
    
    # constructor
    def __init__(self,name="no brand",volume=None):
        self.name = name
        if volume!=None: self.volume = volume
        

# ===================================================   
# main()
# ===================================================   
# for debugging purposes (code called as a script)
# the code is called from here
# ===================================================
if __name__ == '__main__':
    
    F = realfood()
    E95 = ethanol()
    Y = yogurt()
    YF = yogurt(name="danone")
    YF.description = "yogurt with fruits"
    
    # how to set a new food
    class sandwich(realfood,solid,fat): name = "sandwich"
    S = sandwich()