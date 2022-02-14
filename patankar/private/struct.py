#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Matlab-like Structure class

Created on Sun Jan 23 14:19:03 2022

@author: olivier.vitrac@agroparistech.fr

"""

# revision history
# 2022-02-12 fix disp method for empty structures
# 2022-02-12 add type, format

class struct():
    
    type = "struct"
    ftype = "field"
    
    """ mini class behaving as a Matlab Stucture """
    def __init__(self,**kwargs):
        """ constructor """
        self.set(**kwargs)
        
    def set(self,**kwargs):
        """ initialization """
        self.__dict__.update(kwargs)
        
    def setattr(self,key,val):
        """ set field and value """
        self.__dict__[key] = val
        
    def getattr(self,key):
        """ get value """
        return self.__dict__[key]
         
    def __len__(self):
        """ return the number of fields """
        return len(self.__dict__)
       
    def __repr__(self):
        """ display method """
        if self.__dict__=={}:
            print("empty structure with no fields")
            return "empty structure"
        else:
            keylengths = [len(key) for key in self.__dict__]
            width = max(10,max(keylengths)+2)
            fmt = "%%%ss:" % width
            line = ( fmt % ('-'*(width-2)) ) + ( '-'*(min(40,width*5)) )
            print(line)
            for key,value in self.__dict__.items(): print(fmt % key,value)
            print(line)
            return "%s object with %d %ss" % (self.type,len(self.__dict__),self.ftype)
        
    def format(self,s):
        """ format a string with field (use {field} as placeholders) """
        return s.replace("${","{").format(**self.__dict__)
