# Accessing Prediction Models of $D$ and $k$🍏⏩🍎

This document details the different methods to access $D$ and $k$ models used internally by `SFPPy`. This document is useful for the end-user and the developer that would like to add his/her own models.

> Sections 1 and 2 are accessible to a large audience. Sections 3-5 have much more specialized content and discuss how new models or rules can be added to `SFPPy`.
>
> The computational details of these models are not discussed in this document.



[TOC]

## 1 | Overview of models available in `SFPPy`

| Model       | Target | Interpretation                                               |
| ----------- | ------ | ------------------------------------------------------------ |
| `Dpiringer` | $D$    | Default overestimate of $D$ values (available for all classes of polymers, adhesives and as an ersatz for paper and cardboard).<br />This model is empirical. |
| `DFV`       | $D$    | Reference diffusivity model for evaluating recycled materials and recycling processes (available for polyolefins, polyesters, polyamides, polyvinyls ) - currently available for toluene.<br />This model relies on extensive theories. |
| `Dwelle`    | $D$    | Better $D$ model for “dry” (non-plasticized) PET and PS.<br />This model is empirical with similitudes with DFV for scaling. |
| `kFH`       | $k$    | any food and polymer associated with a molecule (solvent) or pattern (monomer) |

where $D$ refers to diffusivites (in $m^2\cdot s^{-1}$), which characterize the capability of a substance to spread in the material due to thermal agitation and $k$ to Henry-like coefficients (relative units), which assess the chemical affinity of the substance (higher $k$ values indicate lower chemical affinity). 



---



## 2 | Principles

### 2.1 | Without molecular information

> $D$ and $k$ models are defined in `patankar.property`. They can be accessed from the class method `evaluate()` by mentioning the name of the model: `Dpiringer`, `DFV`,`Dwelle`,`kFH`…
>
> > Low-level evaluations use a functional form.

```python
from patankar.property import Dpiringer
Dpiringer()
```

yields

```python
   property: Diffusivity
   notation: D
description: Piringer's overestimate of diffusion coefficients
       name: Piringer
parameters: {'polymer': {'polymer': 'polymer code/name', 'units': 'N/A'}, 'M': {'description': 'molecular mass', 'units': 'g/mol'}, 'T': {'description': 'temperature', 'units': 'degC'}}
    SIunits: m**2/s
    license: MIT
    version: 1.3
Out: <Dpiringer: Diffusivity:D>
```

Then the diffusivity (overestimate) in `LDPE` of a migrant with molecular weight $M=200\,\, g\cdot mol^{-1}$ at $T=25^\circ C$ is given by:

```python
D = Dpiringer.evaluate("LDPE",M=200,T=25)
print(D)
```

yields

```python
1.052531538821117e-12
```

### 2.2 | With substance information

> Relevant substance information (molecular weight, molar volume, molecular volume, polarity index…) are automatically retrieved or calculated with `migrant` from `patankar.loadpubchem`. The  `migrant` objects are instantiated with populated `Dtemplate` and `ktemplate` facilitating predictions.
>
> > In this case, the default models (`Dpiringer` and `kFH`) are managed and therefore imported by the class `migrant`. There is no need to import them directly.

In the example below, `BHT` is chosen as an example: 

```python
from patankar.loadpubchem import migrant
m = migrant("BHT")
m
```

yields

```python
<migrant object>
    Compound: BHT
        Name: ['Vianol', 'p-Cresol [...] yl 4-methyl phenol']
         cid: 31404
         CAS: ['128-37-0']
     M (min): 220.35
     M_array: [220.35]
     formula: C15H24O
      smiles: CC1=CC(=C(C(=C1)C(C)(C)C)O)C(C)(C)C
    InChiKey: NLZUEZXRPGMBCV-UHFFFAOYSA-N
        logP: [5.3]
   P' (calc): [0.]
Out: <migrant: Vianol - M=220.35 g/mol>
```

The object `m` includes two templates: `Dtemplate` (for $D$) and `ktemplate`($k$).

```python
m.Dtemplate
```

returns the default $D$ parameters as a `dict` (Python dictionary)

```python
{'polymer': 'LLDPE',
 'M': 220.35,
 'Vvdw': 217.12431614336657,
 'T': 40.0,
 'Tg': 76.0,
 'logP': array([5.3])}
```


```python
m.ktemplate
```

returns the default $k$ parameters as a `dict` (Python dictionary)

```python
{'Pi': array([0.]),
 'Pk': 3.97,
 'Vi': 258.2864390647352,
 'Vk': 30.9,
 'ispolymer': True,
 'alpha': 0.14,
 'lngmin': 0.0,
 'Psat': 1.0,
 'crystallinity': 0,
 'porosity': 0}
```

> These templates can be used in any of the corresponding prediction ($D$ or $k$) models installed SFPPy. 

The default evaluation models are accessible as Python `lambda` functions stored in the attribute `Deval` and `keval`:

```python
m.Deval
```
`Out: <function patankar.loadpubchem.migrant.Deval.locals>.func(**kwargs)`

```python
m.keval
```
`Out: <function patankar.loadpubchem.migrant.keval.<locals>.func(**kwargs)>`

> One possible approach for point-wise estimates is to copy the desired template (with `myD=m.Dtemplate.copy()` or `myk=m.ktemplate.copy()`) and to update the entries inside (with `myD.update(T=...)` or `myk.update(param...)`

We can do that in one single operation:

```python
# here the standard approach
#myD = m.Dtemplate.copy()
#myD.update(polymer="LDPE", T=20)

# Here the shortest syntax
# When a dictionary contains duplicate keys, the last occurrence of a key in the definition takes precedence. 
myD = {**m.Dtemplate, "polymer": "LDPE", "T": 20} # set LDPE @20°C // **means unpacking
DLDPE20 = m.Deval(**myD)
print("D in LDPE at 20°C=",DLDPE20)
```
`D in LDPE at 20°C= 4.5208446252136935e-13`

Now we can reuse myD for a new polymer (e.g., PP) and temperature (e.g., 35°C):

```python
# note that update() is returning None, its output cannot be used directly
myD.update(polymer="PP",T=35) 
DPP35 = m.Deval(**myD)
print("D in PP at 35°C=",DPP35)
```

`D in PP at 35°C= 7.610541702987448e-14`

## 2.3 | With substance and material/medium information

### 2.3.1 | Fully automatic approach (preferred)

When the substance/migrant is incorporated to a `layer` (from `patankar.layer`) or `foodlayer` (from `food.foodlayer`) instance, the system becomes fully informed and no action is expected from the user.

```python
from patankar.layer import LDPE, PP
from patankar.food import ethanol
A = LDPE(T=20,migrant=m)
B = PP(T=35,migrant=m)
F = ethanol(temperaturecontact=20,migrant=m)
# D is a dynamic property which is overriding any default value (_D)
print("D in LDPE at 20°C=",A.D)  
print("D in PP at 35°C=",B.D)
# k0 is a dynamic property which is overriding any default value (k)
print("k0 in ethanol=",F.k0)
# properties A.k and B.k are also available
```

yields

```python
D in LDPE at 20°C= [4.52084463e-13] # same value as DLDPE20
D in PP at 35°C= [7.6105417e-14] # same value as DPP35
k0 in ethanol= [0.25828644] # 0 refers to the index of tehe food layer in our publications
```

### 2.3.2 | Low-level control of models (not recommended)

Each `layer` and `foodlayer` instance offers methods to create lambda functions (*i.e.*, anonymous functions) that can accept `Dtemplate` and `ktemplate`:

```python
A._compute_Dmodel
```

`Out: <function patankar.layer.layer._compute_Dmodel.<locals>.func(**kwargs)>`

```python
A._compute_kmodel
```

`Out: <function patankar.layer.layer._compute_kmodel.<locals>.func(**kwargs)>`

### 2.3.3 | User overrides (advanced users)

`layer` instances offer two attributes `Dmodel` and `kmodel` to attach user functions accepting `Dtemplate` and `ktemplate` dictionaries. These attributes are initialized to `"default"`(str) meaning that the internal `_compute_Dmodel` or `_compute_kmodel` must be used instead.

When the user model returns `None`, the default value $D$ or $k$ values (set with `A.D=...` or `A.k=` ) is used instead:

```python
# Principles
# A.Dmodel=lambda *args, **kwargs: None # to force values set with A.D=...
# A.kmodel=lambda *args, **kwargs: None # to force values set with A.k=...
A.D = (1.234e-8,"cm**2/s")
A.Dmodel = A.Dmodel=lambda *args, **kwargs: None
print("new A.D =",A.D)
```

`new A.D = [1.234e-12]`



---



## 3 | Adding new models to `patankar.property`

The module `patankar.property`offers documentation and examples to add new models of the following types.

|               Base Class | Property Code | Description                                            |
| -----------------------: | :-----------: | ------------------------------------------------------ |
|          `Diffusivities` |      $D$      | Mathematical model to estimate diffusivities           |
|  `HenryLikeCoefficients` |      $k$      | Mathematical model to estimate Henri-like coefficients |
|   `ActivityCoefficients` |      $g$      | Mathematical model to estimate activity coefficients   |
| `PartitionCoeffcicients` |      $K$      | Mathematical model to estimate partition coefficients  |

For each property and model, it is imperative that it include a **class method** `evaluate` that can be called with a `dict` used as a template.



- A global variable `MigrationPropertyModels` indexes all available models

```python
MigrationPropertyModels = {
    "D":{
        "Piringer": Dpiringer,
        "FV": DFV,
        "Welle": Dwelle,
        # add other diffusivity models here
        },
    "k":{
        "FHP": kFHP
        # add other Henry-like models here
        },
    "g":{
        "FHP": gFHP
        # add other activity coefficients models here
        },
    "K":{
        },
    }
```



- Two helper functions are provided to check that the **right model** is used for the **right purpose**:

  - function`MigrationPropertyModel_validator(model)` which is returning True if all conditions are met.

  - function `PropertyModelSelector(rules,...)` that select the best model according to `rules` 



- A model not valid for production can be hidden to other classes by setting the class parameter `_available_to_import`

 ```python
 _available_to_import = False # this model can be directly imported if True
 ```



***



## 4 | Changing the precedence of models recognized by `patankar.loadpubchem`

### 4.1 | Precedence rules

The module `patankar.loadpubchem` uses chained rules (global variables `Dmodel_extensions` and `kmodel_extensions`) to define which alternative model must be used. It is possible to define rules based on the material id, substance id, a specific temperature range… 

```python
# Alternative Dmodels
Dmodel_extensions = {
    "DFV": {
      "description": "hole Free-Volume theory model for toluene in many polymers",
          "objects": ["material","migrant"], # it requires material and migrant (material always first)
            "rules": [ # we assume AND between all conditions
                        # Condition on the material (always first)
                        {"list": [
                            # medium must be a polymer (ispolymer == True)
                            {"attribute": "ispolymer",
                             "op": "istrue",
                            },
                            # the polymer with index must be of these types
                            {"attribute": "layerclass_history",
                             "index":0,
                             "op": "in",
                             "value": ("gPET","wPET","PMMA","PS","PVAc","LDPE")
                            }, # next condition
                                ] # close list of rules for rules[0]
                        }, # close rules[0]
                        # Condition on migrant
                        {"list": [
                            # migrant must be Toluene (based on its InChiKey)
                            {"attribute": "InChiKey",
                                   "op": "==",
                                   "value": "YXFVVABEGXRONW-UHFFFAOYSA-N"
                            }
                                ], # next condition
                        }, # next rule
                    ] # close rules
            }, # next model

    "Dwelle": {
    "description": "Frank-Welle diffusivity model based on VdW volumes",
        "objects": ["material"],
          "rules": [ # we assume AND between all conditions
                    # Condition on the material (always first)
                    {"list": [
                        # medium must be a polymer (ispolymer == True)
                        {"attribute": "ispolymer",
                         "op": "istrue",
                        },
                        # the polymer with index must be of these types
                        {"attribute": "layerclass_history",
                         "index":0,
                         "op": "in",
                         "value": ("gPET","PS","rPS","HIPS","rHIPS")
                        }, # next condition
                            ] # close list of rules for rules[0]
                    }, # close rules[0]
                    ] # close rules
            }
    }

```



### 4.2 Checking rules for a given `layer` index

The rules are managed by two helper methods of `migrant` instances at the level of each layer.

- `suggest_alt_Dmodel(material, index,...)`returns the name of the model as a string
- `suggest_alt_Dclass(material, index...)` returns a type object coding for the Python class of the model

> These methods are used internally by all `layer` instances to check that the model is valid for the requested conditions. If the alternative model does not exists or is returning None, the next available model is used automatically. The default model will be used as fall-back if none of the alternative models can operate.



## 5 | Changing the default models and templates used by `migrant`from `patankar.loadpubchem`

The constructor of the class `migrant`of the module`patankar.loadpubchem` defines the default models and templates via inputs:

- `Dmodel`, `kmodel`: are string names which should match the keys of `MigrationPropertyModels["D"]` and `MigrationPropertyModels["k"]` defined in `patankar.property`

- `Dtemplate` and `ktemplate`: they are `dict`s

  

```python
def __init__(self, name=None,   # substance identified by name
             M=None, logP=None, # substance identified by M, logP(less reliable)

             Dmodel = "Piringer", # <--- default D model

             Dtemplate = {"polymer":"LLDPE",
                          "M":50.0,  # used by Dpiringer (molecular mass in g/mol)
                       "Vvdw":100.0, # used by Dwelle (molecular volume)
                          "T":40.0,  # used by Dpringer, DFV
                         "Tg":76.0,  # used by DFV
                          }, # do not use None

             kmodel = "FHP",    # <--- default k model

             ktemplate = {"Pi":1.41,  # P'i (polarity index)
                          "Pk":3.97,  # P'k (polarity index)
                          "Vi":124.1, # molar volume of i
                          "Vk":30.9,  # molar volume of k
                          "ispolymer":True, # True if FH theory is applicable
                          "alpha":0.14, # \alpha \times (P'_i-P'k)^2
                          "lngmin":0.0, # min of log(\gamma_i) -- see theory at inifinite dilution
                          "Psat":1.0,   # partial saturation pressure for i (usually not defined)
                          "crystallinity":0, # k are calculated respectively to the volume fraction
                          "porosity":0       # of amorphous phase (1-crystallinity)(1-porosity)
                          }, # do not use None

             db=dbdefault, # cache.PubChem database

             raiseerror=True # raise an error if the susbtance is not found

             ):
```



> It is useful to change them only if new attributes/properties are required. Note that they are shared with all models and are not specific to one model (choose carefully their names). 



A low-level helper method `_validate_and_set_model` of instances `migrant` is controlling that entries provided by the user follow the rules of the template.



***

<div style="border: 2px solid #4CAF50; border-radius: 8px; padding: 10px; background: linear-gradient(to right, #4CAF50, #FF4D4D); color: white; text-align: center; font-weight: bold;">
  <span style="font-size: 20px;">🍏⏩🍎 <strong>SFPPy for Food Contact Compliance and Risk Assessment</strong></span><br>
  Contact <a href="mailto:olivier.vitrac@gmail.com" style="color: #fff; text-decoration: underline;">Olivier Vitrac</a> for questions |
  <a href="https://github.com/ovitrac/SFPPy" style="color: #fff; text-decoration: underline;">Website</a> |
  <a href="https://ovitrac.github.io/SFPPy/" style="color: #fff; text-decoration: underline;">Documentation</a>
</div>

