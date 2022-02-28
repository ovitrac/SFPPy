#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Matlab-like Structure class
Matlab MS.alias like class (with MS = Molecular Studio, author Olivier V.)

Created on Sun Jan 23 14:19:03 2022

@author: olivier.vitrac@agroparistech.fr

"""

# revision history
# 2022-02-12 fix disp method for empty structures
# 2022-02-12 add type, format
# 2022-02-19 integration of the derived class param()
# 2022-02-20 code optimization, iterable class- major update
# 2022-02-26 clarify in the help the precedence s=s1+s2

# Dependencies
from math import * # import math to authorize all math expressions in parameters
import types       # to check types
import re          # regular expression
from copy import copy as duplicate # to duplicate objects

# core struct cal
class struct():
    """ 
    mini class behaving as a Matlab Stucture 
    Use param() if you need evaluation as MS.alias()
    
    s=struct(a=1,b=2,c='${a}+${b} # evaluate me if you can')
    generates
      --------:----------------------------------------
             a: 1
             b: 2
             c: ${a}+${b} # evaluate me if you can
      --------:----------------------------------------
    Out: structure (struct object) with 3 fields
    
    
    s=param(a=1,b=2,c='${a}+${b} # evaluate me if you can')
    generates instead:
      --------:----------------------------------------
             a: 1
             b: 2
             c: ${a}+${b} # evaluate me if you can   (= 3 )
      --------:----------------------------------------
    Out: parameter list (param object) with 3 definitions    
        
    
    Data can be appended or changed
        s.a = 10
        s.d = 11
    
    Data can be indexed
        for x in s: print(x)
        s[0]
        s[-1]
    
    To delete a field, use either
        delattr(s,'d') # Python's way
    or
        s.d = []       # à la Matlab

    To check whether a field exist, use either
        hasattr(s,'b') # Python' way
    or
        'b' in s       # other Python's way'
        
    
    Concatenate or remove structures
    Note: the definitions in the most right structure will overwrite existing values
        a=struct(a=1,b=2)
        b=struct(c=3,d="d",e="e")
        c=a+b
        d=a+struct(c=30)
        e=c-a
        f=b-struct(e="")
        b-=struct(e="")
        
    Practical shorthands
    
    > construct a structure from keys
        s = struct.fromkeys(["a","b","c","d"])
          --------:----------------------------------------
                 a: None
                 b: None
                 c: None
                 d: None
          --------:----------------------------------------
        Out: structure (struct object) with 4 fields
        
    > build a structure from variables in a string
        s = struct.scan("${a}+${b}*${c}/${d} --- ${ee}")
        s.a = 1
        s.b = "test"
        s.c = [1,"a",2]
        s.generator()
    leads to
        X = struct(
                 a= 1,
                 b= "test",
                 c= [1, 'a', 2],
                 d= None,
                ee= None
                 )
    """
    
    # attributes to be overdefined
    _type = "struct"        # object type
    _fulltype = "structure" # full name
    _ftype = "field"        # field name
    _evalfeature = False    # true if eval() is available
    
    # attributes for the iterator method
    # Please keep it static, duplicate the object before changing _iter_
    _iter_ = 0              
    
    # excluded attributes (keep the , in the Tupple if it is singleton)
    _excludedattr = ('_iter_','__class__') # used by keys() and len()
    
    
    # Methods
    def __init__(self,**kwargs):
        """ constructor """
        self.set(**kwargs)
        
    def set(self,**kwargs):
        """ initialization """
        self.__dict__.update(kwargs)
        
    def setattr(self,key,value):
        """ set field and value """
        if isinstance(value,list) and len(value)==0 and key in self:
            delattr(self, key)
        else:
            self.__dict__[key] = value
        
    def getattr(self,key):
        """ get value """
        if key in self.keys():
            return self.__dict__[key]
        raise KeyError(f'the {self._ftype} "{key}" does not exist')
    
    def hasattr(self,key):
        """ return true if the field exist """
        return key in self.keys()
    
    def __getstate__(self):
        """ getstate for cooperative inheritance / duplication """
        return self.__dict__.copy()
    
    def __setstate__(self,state):
        """ setstate for cooperative inheritance / duplication """
        self.__dict__.update(state)
    
    def __getattr__(self,key):
        """ get attribure override """
        return self.getattr(key)

    def __setattr__(self,key,value):
        """ set attribute override """
        self.setattr(key,value)
        
    def __contains__(self,item):
        """ in override """
        return self.hasattr(item)
    
    def keys(self):
        """ return the fields """
        # keys() is used by struct() and its iterator
        return [key for key in self.__dict__.keys() if key not in self._excludedattr]
    
    def __getitem__(self,idx):
        """ return the ith element of the object  """
        if idx<len(self):
            return self.getattr(self.keys()[idx])
        raise IndexError(f"the {self._ftype} index should be comprised between 0 and {len(self)-1}")
         
    def __len__(self):
        """ return the number of fields """
        # len() is used by struct() and its iterator
        return len(self.keys())
    
    def __iter__(self):
        """ struct iterator """
        # note that in the original object _iter_ is a static property not in dup
        dup = duplicate(self)
        dup._iter_ = 0
        return dup
    
    def __next__(self):
        """ increment iterator """
        self._iter_ += 1
        if self._iter_<=len(self):
            return self[self._iter_-1]
        self._iter_ = 0
        raise StopIteration(f"Maximum {self._ftype} iteration reached {len(self)}")
        
    def __add__(self,s):
        """ add a structure """
        if not isinstance(s,struct):
            raise TypeError(f"the second operand must be {self._type}")
        dup = duplicate(self)
        dup.__dict__.update(s.__dict__)
        return dup
    
    def __iadd__(self,s):
        """ iadd a structure """
        if not isinstance(s,struct):
            raise TypeError(f"the second operand must be {self._type}")
        self.__dict__.update(s.__dict__)
        return self

    def __sub__(self,s):
        """ sub a structure """
        if not isinstance(s,struct):
            raise TypeError(f"the second operand must be {self._type}")
        dup = duplicate(self)
        listofkeys = dup.keys()
        for k in s.keys():
            if k in listofkeys:
                delattr(dup,k)
        return dup
    
    def __isub__(self,s):
        """ isub a structure """
        if not isinstance(s,struct):
            raise TypeError(f"the second operand must be {self._type}")
        listofkeys = self.keys()
        for k in s.keys():
            if k in listofkeys:
                delattr(self,k)
        return self
    

    def __repr__(self):
        """ display method """
        if self.__dict__=={}:
            print(f"empty {self._fulltype} ({self._type} object) with no {self._type}s")
            return f"empty {self._fulltype}"
        else:
            tmp = self.eval() if self._evalfeature else []
            keylengths = [len(key) for key in self.__dict__]
            width = max(10,max(keylengths)+2)
            fmt = "%%%ss:" % width
            fmteval = fmt[:-1]+"="
            line = ( fmt % ('-'*(width-2)) ) + ( '-'*(min(40,width*5)) )
            print(line)
            for key,value in self.__dict__.items():
                if key not in self._excludedattr:
                    print(fmt % key,value)
                    if self._evalfeature and isinstance(value,str):
                        print(fmteval % "",tmp.getattr(key))
            print(line)
            return f"{self._fulltype} ({self._type} object) with {len(self.__dict__)} {self._ftype}s"
        
    def format(self,s):
        """ format a string with field (use {field} as placeholders) """
        return s.replace("${","{").format(**self.__dict__)

    def fromkeys(self,keys):
        """ returns a structure from keys """
        return self+struct(**dict.fromkeys(keys,None))
    
    @staticmethod
    def scan(s):
        """ scan a string for variables """
        if not isinstance(s,str):
            raise TypeError("scan() requires a string")
        tmp = struct()
        #return tmp.fromkeys(set(re.findall(r"\$\{(.*?)\}",s)))
        found = re.findall(r"\$\{(.*?)\}",s);
        uniq = []
        for x in found:
            if x not in uniq: uniq.append(x)
        return tmp.fromkeys(uniq)

    def generator(self):
        """ generate Python code of the equivalent structure """
        nk = len(self)
        if nk==0:
            print("X = struct()")
        else:
            ik = 0
            fmt = "%%%ss=" % max(10,max([len(k) for k in self.keys()])+2)
            print("\nX = struct(")
            for k in self.keys():
                ik += 1
                end = ",\n" if ik<nk else "\n"+(fmt[:-1] % ")")+"\n"
                v = getattr(self,k)
                if isinstance(v,(int,float)) or v == None:
                    print(fmt % k,v,end=end)
                elif isinstance(v,str):
                    print(fmt % k,f'"{v}"',end=end)
                elif isinstance(v,(list,tuple)):
                    print(fmt % k,v,end=end)
                else:
                    print(fmt % k,"/* unsupported type */",end=end)
                    
    

# meta struct class param for scripting with evaluation capability
class param(struct):
    """ 
    class parameters derived from struct() with dynamic evaluation
    container obj.param = value
    
    Example:
        s=param(a=1,b=2,c='${a}+${b} # evaluate me if you can',
                d="$this is a string",e="1000 # this is my number")
        
    returns
          --------:----------------------------------------
                 a: 1
                 b: 2
                 c: ${a}+${b} # evaluate me if you can
                  = 3
                 d: $this is a string
                  = this is a string
                 e: 1000 # this is my number
                  = 1000
          --------:----------------------------------------
        Out: parameter list (param object) with 5 definitions
        
        s.a=10
        
    produces
          --------:----------------------------------------
                 a: 10
                 b: 2
                 c: ${a}+${b} # evaluate me if you can
                  = 12
                 d: $this is a string
                  = this is a string
                 e: 1000 # this is my number
                  = 1000
          --------:----------------------------------------
      Out: parameter list (param object) with 5 definitions
      
      
     Other example with text parameters
      s = param()
      s.mypath = "$/this/folder" 
      s.myfile = "$file"
      s.myext = "$ext"
      s.fullfile = "$${mypath}/${myfile}.${myext}"
      
    generates
          --------:----------------------------------------
            mypath: $/this/folder
                  = /this/folder
            myfile: $file
                  = file
             myext: $ext
                  = ext
          fullfile: $${mypath}/${myfile}.${myext}
                  = /this/folder/file.ext
          --------:----------------------------------------
    Out: parameter list (param object) with 4 definitions
    
    
    """
    
    # override
    _type = "param"
    _fulltype = "parameter list"
    _ftype = "definition"
    _evalfeature = True    # This class can be evaluated with .eval()    
    
    # magic constructor
    def __init__(self,**kwargs):
        """ constructor """
        super().__init__(**kwargs)
        
    # lines starting with # (hash) are interpreted as comments
    # ${variable} or {variable} are substituted by variable.value
    # any line starting with $ is assumed to be a string (no interpretation)
    # ^ is accepted in formula(replaced by **))
    def eval(self,s=""):
        """ eval method for structure such as MS.alias """
        # Evaluate all DEFINITIONS
        # the argument s is only used by formateval() for error management
        tmp = struct()
        for key,value in self.__dict__.items():
            # strings are assumed to be expressions on one single line
            if isinstance(value,str):
                # replace ${variable} (Bash, Lammps syntax) by {variable} (Python syntax)
                valuesafe = value.strip().replace("${","{")
                # replace "^" (Matlab, Lammps exponent) by "**" (Python syntax)
                valuesafe = valuesafe.replace("^","**")
                # Remove all content after # (they are considered comments)
                poscomment = valuesafe.find("#")
                if poscomment>=0: valuesafe = valuesafe[0:poscomment].strip()
                # Literal string starts with $
                if valuesafe.startswith("$"):
                    tmp.setattr(key,tmp.format(valuesafe[1:].lstrip())) # discard $
                else: # string empty or which can be evaluated
                    if valuesafe=="":
                        tmp.setattr(key,valuesafe) # empty content
                    else:
                        tmp.setattr(key, eval(tmp.format(valuesafe)))
            elif isinstance(value,(int,float)): # already a number
                tmp.setattr(key, value) # store the value with the key
            else: # unsupported types
                if s.find("{"+key+"}")>=0:
                    print(f'*** WARNING ***\n\tIn the {self._ftype}:"\n{s}\n"')
                else:
                    print(f'unable to interpret the "{key}" of type {type(value)}')
        return tmp

    # formateval obeys to following rules
    # lines starting with # (hash) are interpreted as comments
    def formateval(self,s):
        """ format method with evaluation feature """
        tmp = self.eval(s)
        # Do all replacements in s (keep comments)
        if len(tmp)==0:
            return s
        else:
            slines = s.split("\n")
            for i in range(len(slines)):
                poscomment = slines[i].find("#")               
                if poscomment>=0:
                    while (poscomment>0) and (slines[i][poscomment-1]==" "):
                        poscomment -= 1 
                    comment = slines[i][poscomment:len(slines[i])]
                    slines[i]  = slines[i][0:poscomment]
                else:
                    comment = ""
                slines[i] = tmp.format(slines[i])+comment
            return "\n".join(slines)