#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Layer builder for patankar package
# general object is layer()
# derived ones LDPE(),HDPE(),PP()...air()

@author: olivier.vitrac@agroparistech.fr

"""

# -- History ----
# Created on Tue Jan 18 09:14:34 2022
# 2022-01-19 RC
# 2022-01-20 full indexing and simplification
# 2022-01-21 add split()
# 2022-01-22 add child classes for common polymers
# 2022-01-23 full implementation of units 


oneline = "Build multilayer objects"

docstr = """
Build layer(s) for SENSPATANKAR

Example of caller:
    from patankar.layer import layer
    A=layer(D=1e-14,l=50e-6)
    A

"""


# Package Dependencies 
# ====================
# <--  generic packages  -->
import numpy as np
from copy import deepcopy as duplicate
# <--  Internal to patankar package (note they are local)  -->
from private.struct import struct
if 'SIbase' not in dir():
    from private.pint import UnitRegistry as SIbase
    from private.pint import set_application_registry as fixSIbase
    
# Initialize unit conversion (intensive initialization)
# NB: degC and kelvin must be used for temperature
# conversion as obj,valueSI,unitSI = toSI(qSI(numvalue,"unit"))
# conversion as obj,valueSI,unitSI = toSI(qSI("value unit"))
def toSI(q): q=q.to_base_units(); return q,q.m,str(q.u)
NoUnits = 'a.u.'    # value for arbitrary unit
if ("SI" not in locals()) or ("qSI" not in locals()):
    SI = SIbase()      # unit engine
    fixSIbase(SI)      # keep the same instance between calls
    qSI = SI.Quantity  # main unit consersion method from string 
    # constants (usable in layer object methods)
    # define R,T0K,R*T0K,1/(R*T0K) with there SI units
    constants = struct()
    R,constants.R,constants.Runit = toSI(qSI(1,'avogadro_number*boltzmann_constant'))
    T0K,constants.T0K,constants.T0Kunit = toSI(qSI(0,'degC'))
    RT0K,constants.RT0K,constants.RT0Kunit = toSI(R*T0K)
    iRT0K,constants.iRT0K,constants.iRT0Kunit = toSI(1/RT0K)

# Concise data validator with unit convertor to SI
def check_units(value,ProvidedUnits,ExpectedUnits):
    """ check numeric inputs and convert them to SI units """
    if (ProvidedUnits==ExpectedUnits) or (ProvidedUnits==NoUnits) or (ExpectedUnits==None):
        conversion =1               # no conversion needed
        units = ExpectedUnits
    else:
        q0,conversion,units = toSI(qSI(1,ProvidedUnits))
    return np.array([value*conversion]),units
    

# default values (usable in layer object methods)
default = struct(l=1e-6,
                 D=1e-14,
                 k=1,
                 C0=1000,
                 T=25,                   # default temperature (°C)
                 lunit='m',              # do not change them
                 Dunit='m**2/s',         # do not change them
                 kunit=NoUnits,
                 Cunit=NoUnits,
                 layername="my layer",
                 layertype="unknown type",
                 layermaterial="unknown material"
                        )

# Main class definition
# =======================
class layer:
    """ 
        layer class from patankar package
        ...
        strings properties: layername, layertype, layermaterial
        scalar properties: D,k,l,C0
        ...
        Example:
            A = layer(D=1e-14,l=50e-6,layername="layer A",layertype="polymer",layermaterial="PP")
    """
    # --------------------------------------------------------------------
    # PRIVATE PROPERTIES (cannot be changed by the user)
    # __ read only properties
    #  _ private properties (not public)
    # --------------------------------------------------------------------
    __description = "LAYER object"                # description
    __version = 0.2                               # version of sens
    __contact = "olivier.vitrac@agroparistech.fr" # contact person
    _printformat = "%0.4g"   # format to display D, k, l values
    
    # --------------------------------------------------------------------
    # CONSTRUCTOR OF INSTANCE PROPERTIES
    # None = missing numeric value (managed by default)
    # --------------------------------------------------------------------
    def __init__(self,
                 l=None, D=None, k=None, C0=None, # default values will be assigned
                 lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername=default.layername,
                 layertype=default.layertype,
                 layermaterial=default.layermaterial):
        """
        
        Parameters
        ----------
        
        layername : TYPE, optional, string
                    DESCRIPTION. Layer Name. The default is "my layer".
        layertype : TYPE, optional, string
                    DESCRIPTION. Layer Type. The default is "unknown type".
        layermaterial : TYPE, optional, string
                        DESCRIPTION. Material identification . The default is "unknown material". 
        PHYSICAL QUANTITIES
        l : TYPE, optional, scalar
            DESCRIPTION. Thickness. The default is 50e-6 (m).
        D : TYPE, optional, scalar
            DESCRIPTION. Diffusivity. The default is 1e-14 (m^2/s).
        k : TYPE, optional, scalar
            DESCRIPTION. Henry-like coefficient. The default is 1 (a.u.).
        C0 : TYPE, optional, scalar
            DESCRIPTION. Initial concentration. The default is 1000 (a.u.).
        PHYSICAL UNITS
        lunit : TYPE, optional, string
                DESCRIPTION. Length units. The default unit is "m.
        Dunit : TYPE, optional, string
                DESCRIPTION. Diffusivity units. The default unit is 1e-14 "m^2/s".
        kunit : TYPE, optional, string
                DESCRIPTION. Henry-like coefficient. The default unit is "a.u.".
        Cunit : TYPE, optional, string
                DESCRIPTION. Initial concentration. The default unit is "a.u.".
        Returns
        -------
        a monolayer object which can be assembled into a multilayer structure

        """
        if l==None: l=default.l
        if D==None: D=default.D
        if k==None: k=default.k
        if C0==None: C0=default.C0
        if lunit==None: lunit=default.lunit
        if Dunit==None: Dunit=default.Dunit
        if kunit==None: kunit=default.kunit
        if Cunit==None: Cunit=default.Cunit
        l,lunit = check_units(l,lunit,default.lunit)
        D,Dunit = check_units(D,Dunit,default.Dunit)
        k,kunit = check_units(k,kunit,default.kunit)
        C0,Cunit = check_units(C0,Cunit,default.Cunit)
        self._name = [layername]
        self._type = [layertype]
        self._material = [layermaterial]
        self._nlayer = 1
        self._l = l[:1]
        self._D = D[:1]
        self._k = k[:1]
        self._C0 = C0[:1]
        self._lunit = lunit
        self._Dunit = Dunit
        self._kunit = kunit
        self._Cunit = Cunit
        
    # --------------------------------------------------------------------
    # overloading binary addition (note that the output is of type layer)
    # --------------------------------------------------------------------
    def __add__(self,other):
        """ C=A+B | overload + operator """
        if isinstance(other, layer):
            res = duplicate(self)
            for p in ["_name","_type","_material","_nlayer"]:
                setattr(res,p,getattr(self,p)+getattr(other,p))
            for p in ["_l","_D","_k","_C0"]:
                setattr(res,p,np.concatenate((getattr(self,p),getattr(other,p))))
            return res
        else: raise ValueError("invalid layer object")


    # --------------------------------------------------------------------
    # overloading binary multiplication (note that the output is of type layer)
    # --------------------------------------------------------------------
    def __mul__(self,ntimes):
        """ nA = A*n | overload * operator """
        if isinstance(ntimes, int) and ntimes>0:
            res = duplicate(self)
            if ntimes>1:
                for n in range(1,ntimes): res += self
            return res
        else: raise ValueError("multiplicator should be a strictly positive integer")


    # --------------------------------------------------------------------
    # len method
    # --------------------------------------------------------------------   
    def __len__(self):
        """ length method """
        return self._nlayer
    
    # --------------------------------------------------------------------
    # object indexing (get,set) method
    # --------------------------------------------------------------------
    def __getitem__(self,i):
        """ get indexing method """
        res = duplicate(self)
        # check indices
        isscalar = isinstance(i,int)
        if isinstance(i,slice):
            if i.step==None: j = list(range(i.start,i.stop))
            else: j = list(range(i.start,i.stop,i.step))
            res._nlayer = len(j)
        if isinstance(i,int): res._nlayer = 1
        # pick indices for each property
        for p in ["_name","_type","_material","_l","_D","_k","_C0"]:
            content = getattr(self,p)
            try:
                if isscalar: setattr(res,p,content[i:i+1])
                else: setattr(res,p,content[i])
            except IndexError as err:
                print("bad layer object indexing: ",err)
        return res
    
    def __setitem__(self,i,other):
        """ set indexing method """
        # check indices
        if isinstance(i,slice):
            if i.step==None: j = list(range(i.start,i.stop))
            else: j = list(range(i.start,i.stop,i.step))
        elif isinstance(i,int): j = [i]
        else:raise IndexError("invalid index")        
        islayer = isinstance(other,layer)
        isempty = not islayer and isinstance(other,list) and len(other)<1
        if isempty:         # empty right hand side
            for p in ["_name","_type","_material","_l","_D","_k","_C0"]:
                content = getattr(self,p)
                try:
                    newcontent = [content[k] for k in range(self._nlayer) if k not in j]
                except IndexError as err:
                    print("bad layer object indexing: ",err)
                if isinstance(content,np.ndarray) and not isinstance(newcontent,np.ndarray):
                    newcontent = np.array(newcontent)
                setattr(self,p,newcontent)
            self._nlayer = len(newcontent)
        elif islayer:        # islayer right hand side
            nk1 = len(j)
            nk2 = other._nlayer
            if nk1 != nk2:
                raise IndexError("the number of elements does not match the number of indices")
            for p in ["_name","_type","_material","_l","_D","_k","_C0"]:
                content1 = getattr(self,p)
                content2 = getattr(other,p)
                for k in range(nk1):
                    try:
                        content1[j[k]] = content2[k]
                    except IndexError as err:
                        print("bad layer object indexing: ",err)
                setattr(self,p,content1)            
        else:
            raise ValueError("only [] or layer object are accepted")
    
    
    # --------------------------------------------------------------------
    # Getter methods (show private/hidden properties and meta-properties)
    # --------------------------------------------------------------------
    @property
    def name(self): return self._name
    @property
    def type(self): return self._type
    @property
    def material(self): return self._material
    @property
    def l(self): return self._l
    @property
    def D(self): return self._D
    @property
    def k(self): return self._k
    @property
    def C0(self): return self._C0
    @property
    def lunit(self): return self._lunit
    @property
    def Dunit(self): return self._Dunit
    @property
    def kunit(self): return self._kunit
    @property
    def Cunit(self): return self._Cunit
    @property
    def n(self): return self._nlayer
    @property
    def resistance(self): return self.l*self.k/self.D
    @property
    def permeability(self): return self.D/(self.l*self.k)
    @property
    def lag(self): return self.l**2/(6*self.D)
    @property
    def pressure(self): return self.k*self.C0
    @property
    def thickness(self): return sum(self.l)
    @property
    def concentration(self): return sum(self.l*self.C0)/self.thickness
    @property
    def relative_thickness(self): return self.l/self.thickness
    @property
    def relative_resistance(self): return self.resistance/sum(self.resistance)
    @property
    def rank(self): return np.flip(np.argsort(np.array(self.resistance))+1).tolist()
    @property
    def referencelayer(self): return np.argmax(self.resistance)
    @property
    def Foscale(self): return self.D[self.referencelayer]/self.l[self.referencelayer]**2
    
    # --------------------------------------------------------------------
    # comparators based resistance
    # --------------------------------------------------------------------
    def __eq__(self, o):
        value1 = self.resistance if self._nlayer>1 else self.resistance[0]
        if isinstance(o,layer):
            value2 = o.resistance if o._nlayer>1 else o.resistance[0]
        else:
            value2 = o
        return value1==value2
    
    def __ne__(self, o):
        value1 = self.resistance if self._nlayer>1 else self.resistance[0]
        if isinstance(o,layer):
            value2 = o.resistance if o._nlayer>1 else o.resistance[0]
        else:
            value2 = o
        return value1!=value2
    
    def __lt__(self, o):
        value1 = self.resistance if self._nlayer>1 else self.resistance[0]
        if isinstance(o,layer):
            value2 = o.resistance if o._nlayer>1 else o.resistance[0]
        else:
            value2 = o
        return value1<value2

    def __gt__(self, o):
        value1 = self.resistance if self._nlayer>1 else self.resistance[0]
        if isinstance(o,layer):
            value2 = o.resistance if o._nlayer>1 else o.resistance[0]
        else:
            value2 = o
        return value1>value2

    def __le__(self, o):
        value1 = self.resistance if self._nlayer>1 else self.resistance[0]
        if isinstance(o,layer):
            value2 = o.resistance if o._nlayer>1 else o.resistance[0]
        else:
            value2 = o
        return value1<=value2

    def __ge__(self, o):
        value1 = self.resistance if self._nlayer>1 else self.resistance[0]
        if isinstance(o,layer):
            value2 = o.resistance if o._nlayer>1 else o.resistance[0]
        else:
            value2 = o
        return value1>=value2

    
    # --------------------------------------------------------------------
    # Getter methods and tools to validate inputs checknumvalue and checktextvalue
    # --------------------------------------------------------------------
    def checknumvalue(self,value):
        """ returns a validate value to set properties """
        if isinstance(value,int): value = float(value)
        if isinstance(value,float): value = np.array([value])
        if isinstance(value,list): value = np.array(value)
        if len(value)>self._nlayer:
            value = value[:self._nlayer]
            print('dimension mismatch, the extra value(s) has been removed')
        elif len(value)<self._nlayer:
            value = np.concatenate((value,value[-1:]*np.ones(self._nlayer-len(value))))
            print('dimension mismatch, the last value has been repeated')
        return value

    def checktextvalue(self,value):
        """ returns a validate value to set properties """
        if not isinstance(value,list): value = [value]
        if len(value)>self._nlayer:
            value = value[:self._nlayer]
            print('dimension mismatch, the extra entry(ies) has been removed')
        elif len(value)<self._nlayer:
            value = value + value[-1:]*(self._nlayer-len(value))
            print('dimension mismatch, the last entry has been repeated')
        return value        
        
    @l.setter
    def l(self,value): self._l =self.checknumvalue(value)
    @D.setter
    def D(self,value): self._D=self.checknumvalue(value)
    @k.setter
    def k(self,value): self._k =self.checknumvalue(value)
    @C0.setter
    def C0(self,value): self._C0 =self.checknumvalue(value)
    @name.setter
    def name(self,value): self._name =self.checktextvalue(value)
    @type.setter
    def type(self,value): self._type =self.checktextvalue(value)
    @material.setter
    def material(self,value): self._material =self.checktextvalue(value)
    
        
    # --------------------------------------------------------------------
    # hash methods (assembly and layer-by-layer)
    # note that list needs to be converted into tuples to be hashed
    # --------------------------------------------------------------------
    def __hash__(self):
        """ hash layer-object (assembly) method """
        return hash((tuple(self._name),
                     tuple(self._type),
                     tuple(self._material),
                     tuple(self._l),
                     tuple(self._D),
                     tuple(self.k),
                     tuple(self._C0)))
    
    # layer-by-layer @property = decoration to consider it 
    # as a property instead of a method/attribute
    # comprehension for n in range(self._nlayer) applies it to all layers
    @property
    def hashlayer(self):
        """ hash layer (layer-by-layer) method """
        return [hash((self._name[n],
                      self._type[n],
                      self._material[n],
                      self._l[n],
                      self._D[n],
                      self.k[n],
                      self._C0[n]))
                for n in range(self._nlayer)
                ]
 
    
    # --------------------------------------------------------------------
    # disp method (since the getter are defined, the '_' is dropped)
    # --------------------------------------------------------------------
    def __repr__(self):
        """ disp method """
        print("[%s version=%0.4g, contact=%s]" % (self.__description,self.__version,self.__contact))
        if self._nlayer==0:
            ret="empty %s" % (self.__description)
        else:
            for n in range(1,self._nlayer+1):
                print('-- [ layer %d of %d ] ---------- barrier rank=%d --------------'
                      % (n,self._nlayer,self.rank[n-1]))
                for p in ["name","type","material"]:
                    v = getattr(self,p)
                    print('%10s: "%s"' % (p,v[n-1]),flush=True)
                for p in ["l","D","k","C0"]:
                    v = getattr(self,p)                 # value
                    vunit = getattr(self,p[0]+"unit")   # value unit
                    print(('%10s: '+self._printformat+" [%s]")
                          % (p,v[n-1],vunit),flush=True)
            if self._nlayer==1:
                ret=('monolayer of %s' % self.__description)
            else:
                ret=('%d-multilayer of %s' % (self._nlayer,self.__description))
        return ret        

    # --------------------------------------------------------------------
    # STRUCT method - returns the equivalent dictionary from an object
    # source: https://stackoverflow.com/questions/61517/python-dictionary-from-an-objects-fields
    # --------------------------------------------------------------------        
    def struct(self):
        """ returns the equivalent dictionary from an object """
        return dict((key, getattr(self, key)) for key in dir(self) if key not in dir(self.__class__))
    # --------------------------------------------------------------------
    
    # --------------------------------------------------------------------
    # SIMPLIFY method ny collecting similar layers
    # --------------------------------------------------------------------        
    def simplify(self):
        """ merge continuous layers of the same type """
        nlayer = self._nlayer
        if nlayer>1:
           res = self[0]
           ires = 0
           ireshash = res.hashlayer[0]
           for i in range(1,nlayer):
               if self.hashlayer[i]==ireshash:
                   res.l[ires] = res.l[ires]+self.l[i]
               else:
                   res = res + self[i]
                   ires = ires+1
                   ireshash = self.hashlayer[i]
        else:
             res = duplicate(self)
        return res
    
    # --------------------------------------------------------------------
    # Split layers into a tuple
    # --------------------------------------------------------------------        
    def split(self):
        """ split layers """
        out = ()
        if self._nlayer>0:
            for i in range(self._nlayer):
                out = out + (self[i],) # (,) special syntax for tuple singleton
        return out


"""
=======================================================
Child classes derived from layer
this section can be extended to define specific layers
    * polymer
    * ink
    * air
    * paper and board
  
These classes are more flexible than the parent class layer.
They can include temperature dependence, refined tunning, etc.    
  
Properties taken from
    * linear thermal expansoin
    https://omnexus.specialchem.com/polymer-properties/properties/coefficient-of-linear-thermal-expansion  
  
Once the layers are incorporated in a multilayer structure,
they loose their original subclass and become only an object
layer. These subclasses are therefore useful to refine the
properties of each layer before standarizing them.

=========================================================
"""

# <<<<<<<<<<<<<<<<<<<<<<< P O L Y O L E F I N S >>>>>>>>>>>>>>>>>>>>>>

# <-- LDPE polymer ---------------------------------->
class LDPE(layer):
    """  extended pantankar.layer for low-density polyethylene LDPE  """
    def __init__(self,l=100e-6,D=1e-12,
                 k=None,C0=None,lunit=None,Dunit=None,kunit=None,Cunit=None,
                 layername="layer in LDPE"):
        """ LDPE layer constructor """
        layer.__init__(self,
                       l=l,D=D,k=k,C0=C0,
                       lunit=lunit,Dunit=Dunit,kunit=kunit,Cunit=Cunit,
                       layername=layername,
                       layertype="polymer", # set by default at inititialization
                       layermaterial="low-density polyethylene")
    def density(self,T=25):
        """ density of LDPE: density(T in K) """
        return 920 *(1-3*(T-default.T)*20e-5),"kg/m**3" # lowest temperature
    @property
    def Tg(self):
        """ glass transition temperature of LDPE """
        return -130,"C" # lowest temperature

        
# <-- HDPE polymer ---------------------------------->
class HDPE(layer):
    """  extended pantankar.layer for high-density polyethylene HDPE  """
    def __init__(self,l=500e-6,D=1e-13,
                 k=None,C0=None,lunit=None,Dunit=None,kunit=None,Cunit=None,
                 layername="layer in HDPE"):
        """ HDPE layer constructor """
        layer.__init__(self,
                       l=l,D=D,k=k,C0=C0,
                       lunit=lunit,Dunit=Dunit,kunit=kunit,Cunit=Cunit,
                       layername=layername,
                       layertype="polymer", # set by default at inititialization
                       layermaterial="high-density polyethylene")
    def density(self,T=25):
        """ density of HDPE: density(T in K) """
        return 940 *(1-3*(T-default.T)*11e-5),"kg/m**3" # lowest temperature
    @property
    def Tg(self):
        """ glass transition temperature of HDPE """
        return -100,"C" # highest temperature
        
# <-- PP polymer ---------------------------------->
class PP(layer):
    """  extended pantankar.layer for isotactic polypropylene PP  """
    def __init__(self,l=300e-6,D=1e-14,
                 k=None,C0=None,lunit=None,Dunit=None,kunit=None,Cunit=None,
                 layername="layer in PP"):
        """ PP layer constructor """
        layer.__init__(self,
                       l=l,D=D,k=k,C0=C0,
                       lunit=lunit,Dunit=Dunit,kunit=kunit,Cunit=Cunit,
                       layername=layername,
                       layertype="polymer", # set by default at inititialization
                       layermaterial="isotatic polypropylene")
    def density(self,T=25):
        """ density of PP: density(T in K) """
        return 910 *(1-3*(T-default.T)*7e-5),"kg/m**3" # lowest temperature
    @property
    def Tg(self):
        """ glass transition temperature of PP """
        return 0,"C" # highest temperature        


# <<<<<<<<<<<<<<<<<<<<<<< G A S E S  >>>>>>>>>>>>>>>>>>>>>>

# <-- air | ideal gas layer ---------------------------------->
class air(layer):
    """  extended pantankar.layer for ideal gases such as air """
    def __init__(self,l=1e-2,D=1e-6,T=25,
                 lunit=None,Dunit=None,Cunit=None,
                 layername="air layer"):
        """ air layer constructor """
        kair = 1/(constants.R *(T+constants.T0K))
        kairunit = constants.iRT0Kunit
        layer.__init__(self,
                       l=l,D=D,k=kair,C0=0,
                       lunit=lunit,Dunit=Dunit,kunit=kairunit,Cunit=Cunit,
                       layername=layername,
                       layertype="air", # set by default at inititialization
                       layermaterial="ideal gas")
  
    
# ===================================================   
# main()
# ===================================================   
# for debugging purposes (code called as a script)
# the code is called from here
# ===================================================
if __name__ == '__main__':
    G = air(T=60)
    P = LDPE(D=1e-8,Dunit='cm**2/s')
    A = LDPE()
    A=layer(D=1e-14,l=50e-6)
    print(A)
    A
    B=A*3
    D = B[1:2]
    B=A+A
    C=B[1]
    B.l = [1,2]
    A.struct()
    E = B*4
    #E[1:4]=[]
    E
    # ======
    A = layer(layername = "layer A")
    B = layer(layername = "layer B")
    C = layer(layername = "layer C")
    D = layer(layername = "layer D")
    # test = A+B+C+D
    # test[2] = test[0]
    # test[3] = []
    # test
    test = A+A+B+B+B+C
    print(test)
    testsimple = test.simplify()
    print(testsimple)