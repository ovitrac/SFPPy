#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Matlab-like Structure class

Created on Sun Jan 23 14:19:03 2022

@author: olivier.vitrac@agroparistech.fr

"""

class struct():
    """ mini class behaving as a Matlab Stucture """
    def __init__(self,**kwargs):
        self.Set(**kwargs)
    def Set(self,**kwargs):
        self.__dict__.update(kwargs)
    def SetAttr(self,key,val):
        self.__dict__[key] = val
    def GetAttr(self,key):
        return self.__dict__[key]
    def __repr__(self):
        keylengths = [len(key) for key in self.__dict__]
        width = max(10,max(keylengths)+2)
        fmt = "%%%ss:" % width
        line = ( fmt % ('-'*(width-2)) ) + ( '-'*(min(40,width*5)) )
        print(line)
        for key,value in self.__dict__.items(): print(fmt % key,value)
        print(line)
        return "structure object with %d fields" % len(self.__dict__)