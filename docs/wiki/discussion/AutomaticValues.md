# Discussion: How diffusivities $D$ and Henry-like coefficients $k$ are automatically set in `SFPPy`🍏⏩🍎

## Preamble
Diffusivities $D_{i,k}$ and $k_{i,k}$ are  properties that depend on the pair substance (indexed $i$) and on the medium where  sorption/diffusion/desorption occurs (indexed $k$). For diffusivities, the medium is limited to food layers indexed in our publications from $1$ to $n$ (in `SFFPy` for `layer`s indexed from `0` to `n-1`). For Henry-like coefficients, the index range of $k$ includes also the food/food simulant indexed $0$ (in `SFPPy` denoted `foodlayer`).


## Update the SFPPy path 


```python
import sys

# Check if SFPPy path is already in sys.path
SFPPy_path = '/replace/it/by/your/path' # update it if needed
if SFPPy_path in sys.path:
    print(f"SFPPy installation '{SFPPy_path}' is already in sys.path.")
else:
    print(f"Add SFPPy path '{SFPPy_path}' to sys.path.")
    sys.path.append(SFPPy_path)
```

    Add SFPPy path '/replace/it/by/your/path' to sys.path.


## Let's define four migrants: $i=0..3$ <small>(Python indexing)</small>


```python
from patankar.loadpubchem import migrant
list_of_migrants = ["toluene","limonene","BHT","Irganox 1076"]
m = [migrant(s) for s in list_of_migrants]
print(m)
```

    <migrant object>
        Compound: toluene
            Name: ['phenyl-methane', ' [...]  for HPLC, >=99.9%']
             cid: 1140
             CAS: ['25013-04-1', '108-88-3']
         M (min): 92.14
         M_array: [92.14]
         formula: C7H8
          smiles: CC1=CC=CC=C1
        InChiKey: YXFVVABEGXRONW-UHFFFAOYSA-N
            logP: [2.7]
       P' (calc): [1.23097093]
       
    <migrant object>
        Compound: limonene
            Name: ['1-Methyl-4-isoprop [...] pentene, p.a., 95%']
             cid: 22311
             CAS: ['138-86-3', '7705-14-8']
         M (min): 136.23
         M_array: [136.23]
         formula: C10H16
          smiles: CC1=CCC(CC1)C(=C)C
        InChiKey: XMGQYMWWDOXHJM-UHFFFAOYSA-N
            logP: [3.4]
       P' (calc): [0.]
       
    <migrant object>
        Compound: BHT
            Name: ['2,6-di-t butyl-4-m [...] ', 'Antioxidant 29']
             cid: 31404
             CAS: ['128-37-0']
         M (min): 220.35
         M_array: [220.35]
         formula: C15H24O
          smiles: CC1=CC(=C(C(=C1)C(C)(C)C)O)C(C)(C)C
        InChiKey: NLZUEZXRPGMBCV-UHFFFAOYSA-N
            logP: [5.3]
       P' (calc): [0.]
       
    <migrant object>
        Compound: Irganox 1076
            Name: ['Octadecyl 3(3,5dit [...] e', 'Tominokusu SS']
             cid: 16386
             CAS: ['2082-79-3']
         M (min): 530.9
         M_array: [530.9]
         formula: C35H62O3
          smiles: CCCCCCCCCCCCCCCCCCOC [...] )C(C)(C)C)O)C(C)(C)C
        InChiKey: SSDSCDGVMJFTEQ-UHFFFAOYSA-N
            logP: [13.8]
       P' (calc): [0.]
       
    [<migrant: phenyl-methane - M=92.14 g/mol>, <migrant: 1-Methyl [...] lohexene - M=136.23 g/mol>, <migrant: 2,6-di-t [...] ylphenol - M=220.35 g/mol>, <migrant: Octadecy [...] opionate - M=530.9 g/mol>]


## Let's define four polymers and a temperature $T=50°C$
- PET with two states: dry PET (`gPET` for glassy PET) and plasticized PET (`wPET` for *wet* PET)
- PS (general purpose PS, unplasticized)
- PP (atactic polypropylene)


```python
from patankar.layer import gPET, wPET, PS, PP
materials = gPET()+wPET()+PS()+PP()
```

## Now we incorporate $i=0$ in all layers `materials`

>Without inserting a substance in `multilayer`, all layers have the same diffusivity:
`D: 1e-14 [m**2/s]`
>
> The default $k=1$ value is also assigned by default to all layers. 


```python
print("D=",materials.D) # we print the default values assigned
materials
```

    D= [1.e-14 1.e-14 1.e-14 1.e-14]
    
    [LAYER object version=1, contact=olivier.vitrac@agroparistech.fr]
    4-multilayer of LAYER object:
    -- [ layer 1 of 4 ] ---------- barrier rank=2 --------------
          name: "layer in gPET"
          type: "polymer"
      material: "glassy PET"
          code: "PET"
             l: 0.0002 [m]
             D: 1e-14 [m**2/s]
             k: 1 [a.u.]
            C0: 1000 [a.u.]
    -- [ layer 2 of 4 ] ---------- barrier rank=4 --------------
          name: "layer in wPET"
          type: "polymer"
      material: "plasticized PET"
          code: "PET"
             l: 0.0002 [m]
             D: 1e-14 [m**2/s]
             k: 1 [a.u.]
            C0: 1000 [a.u.]
    -- [ layer 3 of 4 ] ---------- barrier rank=3 --------------
          name: "layer in PS"
          type: "polymer"
      material: "polystyrene"
          code: "PS"
             l: 0.0001 [m]
             D: 1e-14 [m**2/s]
             k: 1 [a.u.]
            C0: 1000 [a.u.]
    -- [ layer 4 of 4 ] ---------- barrier rank=1 --------------
          name: "layer in PP"
          type: "polymer"
      material: "isotactic polypropylene"
          code: "PP"
             l: 0.0003 [m]
             D: 1e-14 [m**2/s]
             k: 1 [a.u.]
            C0: 1000 [a.u.]





    <multilayer with 4 layers: ['layer in gPET', 'layer in wPET', 'layer in PS', 'layer in PP']>



### We insert `m[0]` in all layers

> All values of $D$ and $k$ are updated automatically;


```python
materials.update(migrant=m[0],T=50)
print(materials.D) # we print the default values automatically calculated at temperature T
materials
```

    dimension mismatch, the last value has been repeated
    [1.86299457e-17 1.24073529e-15 2.64579177e-20 2.78163880e-12]
    
    [LAYER object version=1, contact=olivier.vitrac@agroparistech.fr]
    4-multilayer of LAYER object:
    -- [ layer 1 of 4 ] ---------- barrier rank=1 --------------
          name: "layer in gPET"
          type: "polymer"
      material: "glassy PET"
          code: "PET"
       crystal: 0.35
            Tg: 76 [degC]
             l: 0.0002 [m]
             D: 1.863e-17 [m**2/s]
              = DFV(gPET,<migrant: phenyl-methane - M=92.14 g/mol>,T=50.0 degC)
             k: 0.565 [a.u.]
              = kFHP(<ethylene terephthalate>,<migrant: phenyl-methane - M=92.14 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]
    -- [ layer 2 of 4 ] ---------- barrier rank=3 --------------
          name: "layer in wPET"
          type: "polymer"
      material: "plasticized PET"
          code: "PET"
       crystal: 0.35
            Tg: 46 [degC]
             l: 0.0002 [m]
             D: 1.241e-15 [m**2/s]
              = DFV(wPET,<migrant: phenyl-methane - M=92.14 g/mol>,T=50.0 degC)
             k: 0.565 [a.u.]
              = kFHP(<ethylene terephthalate>,<migrant: phenyl-methane - M=92.14 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]
    -- [ layer 3 of 4 ] ---------- barrier rank=4 --------------
          name: "layer in PS"
          type: "polymer"
      material: "polystyrene"
          code: "PS"
       crystal: 0
            Tg: 100 [degC]
             l: 0.0001 [m]
             D: 2.646e-20 [m**2/s]
              = DFV(PS,<migrant: phenyl-methane - M=92.14 g/mol>,T=50.0 degC)
             k: 0.3672 [a.u.]
              = kFHP(<styrene>,<migrant: phenyl-methane - M=92.14 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]
    -- [ layer 4 of 4 ] ---------- barrier rank=2 --------------
          name: "layer in PP"
          type: "polymer"
      material: "isotactic polypropylene"
          code: "PP"
       crystal: 0.5
            Tg: 0 [degC]
             l: 0.0003 [m]
             D: 2.782e-12 [m**2/s]
              = Dpiringer(PP,<migrant: phenyl-methane - M=92.14 g/mol>,T=50.0 degC)
             k: 1.106 [a.u.]
              = kFHP(<propylene>,<migrant: phenyl-methane - M=92.14 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]



    <multilayer with 4 layers: ['layer in gPET', 'layer in wPET', 'layer in PS', 'layer in PP']>



### Interpretation
> We note that the diffusivities are different for all four layers:
> >`[1.86299457e-17 1.24073529e-15 2.64579177e-20 2.78163880e-12]`
>
> The output is giving details on the models, which were used
> | Reported model                                               | Model       | Interpretation                                               |
> | ------------------------------------------------------------ | ----------- | ------------------------------------------------------------ |
> | `DFV(gPET,<migrant: Toluene-13C - M=92.14 g/mol>,T=50.0 degC)` | `DFV`       | best model for toluene when Free-volume parameters are available for this polymer |
> | `DFV(wPET,<migrant: Toluene-13C - M=92.14 g/mol>,T=50.0 degC)` | `DFV`       | best  model for toluene when Free-volume parameters are available for this polymer |
> | `DFV(gPS,<migrant: Toluene-13C - M=92.14 g/mol>,T=50.0 degC)` | `DFV`       | best  model for toluene when Free-volume parameters are available for this polymer |
> | `Dpiringer(PP,<migrant: Toluene-13C - M=92.14 g/mol>,T=50.0 degC)` | `Dpiringer` | fall-back mechanism                                          |

### We insert `m[1]` in all layers


```python
materials.update(migrant=m[1],T=50)
print(materials.D) # we print the default values automatically calculated at temperature T
materials
```

    dimension mismatch, the last value has been repeated
    [2.70631048e-19 1.72090422e-15 3.74480502e-21 1.39807262e-12]
    
    [LAYER object version=1, contact=olivier.vitrac@agroparistech.fr]
    4-multilayer of LAYER object:
    -- [ layer 1 of 4 ] ---------- barrier rank=1 --------------
          name: "layer in gPET"
          type: "polymer"
      material: "glassy PET"
          code: "PET"
       crystal: 0.35
            Tg: 76 [degC]
             l: 0.0002 [m]
             D: 2.706e-19 [m**2/s]
              = Dwelle(gPET,<migrant: 1-Methyl [...] lohexene - M=136.23 g/mol>,T=50.0 degC)
             k: 0.8452 [a.u.]
              = kFHP(<ethylene terephthalate>,<migrant: 1-Methyl [...] lohexene - M=136.23 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]
    -- [ layer 2 of 4 ] ---------- barrier rank=3 --------------
          name: "layer in wPET"
          type: "polymer"
      material: "plasticized PET"
          code: "PET"
       crystal: 0.35
            Tg: 46 [degC]
             l: 0.0002 [m]
             D: 1.721e-15 [m**2/s]
              = Dpiringer(wPET,<migrant: 1-Methyl [...] lohexene - M=136.23 g/mol>,T=50.0 degC)
             k: 0.8452 [a.u.]
              = kFHP(<ethylene terephthalate>,<migrant: 1-Methyl [...] lohexene - M=136.23 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]
    -- [ layer 3 of 4 ] ---------- barrier rank=4 --------------
          name: "layer in PS"
          type: "polymer"
      material: "polystyrene"
          code: "PS"
       crystal: 0
            Tg: 100 [degC]
             l: 0.0001 [m]
             D: 3.745e-21 [m**2/s]
              = Dwelle(PS,<migrant: 1-Methyl [...] lohexene - M=136.23 g/mol>,T=50.0 degC)
             k: 0.5494 [a.u.]
              = kFHP(<styrene>,<migrant: 1-Methyl [...] lohexene - M=136.23 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]
    -- [ layer 4 of 4 ] ---------- barrier rank=2 --------------
          name: "layer in PP"
          type: "polymer"
      material: "isotactic polypropylene"
          code: "PP"
       crystal: 0.5
            Tg: 0 [degC]
             l: 0.0003 [m]
             D: 1.398e-12 [m**2/s]
              = Dpiringer(PP,<migrant: 1-Methyl [...] lohexene - M=136.23 g/mol>,T=50.0 degC)
             k: 4.322 [a.u.]
              = kFHP(<propylene>,<migrant: 1-Methyl [...] lohexene - M=136.23 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]



    <multilayer with 4 layers: ['layer in gPET', 'layer in wPET', 'layer in PS', 'layer in PP']>



### Interpretation

> We note that the diffusivities are different for all four layers:
>
> >`[2.70631048e-19 1.72090422e-15 3.74480502e-21 1.39807262e-12]`
>
> The output is giving details on the models, which were used
>
> | Reported model                                               | Model       | Interpretation                                        |
> | ------------------------------------------------------------ | ----------- | ----------------------------------------------------- |
> | `Dwelle(gPET,<migrant: Cyclohex [...] )-, (R)- - M=136.23 g/mol>,T=50.0 degC)` | `Dwelle`    | best model for glassy PET when `DFV` is not available |
> | `Dpiringer(wPET,<migrant: Cyclohex [...] )-, (R)- - M=136.23 g/mol>,T=50.0 degC)` | `Dpiringer` | fall-back mechanism                                   |
> | `Dwelle(PS,<migrant: Cyclohex [...] )-, (R)- - M=136.23 g/mol>,T=50.0 degC)` | `Dwelle`    | best model for glassy PS when `DFV` is not available  |
> | `Dpiringer(PP,<migrant: Cyclohex [...] )-, (R)- - M=136.23 g/mol>,T=50.0 degC)` | `Dpiringer` | fall-back mechanism                                   |



## We can repeat the mechanism for `m[2]` and `m[3]`


```python
print("D for",m[2].compound,"=",materials.update(migrant=m[2]).D)
print("D for",m[3].compound,"=",materials.update(migrant=m[3]).D)
```

    D for BHT = [2.13915138e-20 5.73730849e-16 8.90819931e-24 4.66102284e-13]
    D for Irganox 1076 = [1.88965546e-23 2.87174640e-17 4.78703852e-31 2.33302351e-14]



```python
materials.update(migrant=m[2])
```


    [LAYER object version=1, contact=olivier.vitrac@agroparistech.fr]
    4-multilayer of LAYER object:
    -- [ layer 1 of 4 ] ---------- barrier rank=1 --------------
          name: "layer in gPET"
          type: "polymer"
      material: "glassy PET"
          code: "PET"
       crystal: 0.35
            Tg: 76 [degC]
             l: 0.0002 [m]
             D: 2.139e-20 [m**2/s]
              = Dwelle(gPET,<migrant: 2,6-di-t [...] ylphenol - M=220.35 g/mol>,T=50.0 degC)
             k: 1.387 [a.u.]
              = kFHP(<ethylene terephthalate>,<migrant: 2,6-di-t [...] ylphenol - M=220.35 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]
    -- [ layer 2 of 4 ] ---------- barrier rank=3 --------------
          name: "layer in wPET"
          type: "polymer"
      material: "plasticized PET"
          code: "PET"
       crystal: 0.35
            Tg: 46 [degC]
             l: 0.0002 [m]
             D: 5.737e-16 [m**2/s]
              = Dpiringer(wPET,<migrant: 2,6-di-t [...] ylphenol - M=220.35 g/mol>,T=50.0 degC)
             k: 1.387 [a.u.]
              = kFHP(<ethylene terephthalate>,<migrant: 2,6-di-t [...] ylphenol - M=220.35 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]
    -- [ layer 3 of 4 ] ---------- barrier rank=4 --------------
          name: "layer in PS"
          type: "polymer"
      material: "polystyrene"
          code: "PS"
       crystal: 0
            Tg: 100 [degC]
             l: 0.0001 [m]
             D: 8.908e-24 [m**2/s]
              = Dwelle(PS,<migrant: 2,6-di-t [...] ylphenol - M=220.35 g/mol>,T=50.0 degC)
             k: 0.9015 [a.u.]
              = kFHP(<styrene>,<migrant: 2,6-di-t [...] ylphenol - M=220.35 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]
    -- [ layer 4 of 4 ] ---------- barrier rank=2 --------------
          name: "layer in PP"
          type: "polymer"
      material: "isotactic polypropylene"
          code: "PP"
       crystal: 0.5
            Tg: 0 [degC]
             l: 0.0003 [m]
             D: 4.661e-13 [m**2/s]
              = Dpiringer(PP,<migrant: 2,6-di-t [...] ylphenol - M=220.35 g/mol>,T=50.0 degC)
             k: 7.093 [a.u.]
              = kFHP(<propylene>,<migrant: 2,6-di-t [...] ylphenol - M=220.35 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]





    <multilayer with 4 layers: ['layer in gPET', 'layer in wPET', 'layer in PS', 'layer in PP']>




```python
materials.update(migrant=m[3])
```


    [LAYER object version=1, contact=olivier.vitrac@agroparistech.fr]
    4-multilayer of LAYER object:
    -- [ layer 1 of 4 ] ---------- barrier rank=1 --------------
          name: "layer in gPET"
          type: "polymer"
      material: "glassy PET"
          code: "PET"
       crystal: 0.35
            Tg: 76 [degC]
             l: 0.0002 [m]
             D: 1.89e-23 [m**2/s]
              = Dwelle(gPET,<migrant: Octadecy [...] opionate - M=530.9 g/mol>,T=50.0 degC)
             k: 3.431 [a.u.]
              = kFHP(<ethylene terephthalate>,<migrant: Octadecy [...] opionate - M=530.9 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]
    -- [ layer 2 of 4 ] ---------- barrier rank=3 --------------
          name: "layer in wPET"
          type: "polymer"
      material: "plasticized PET"
          code: "PET"
       crystal: 0.35
            Tg: 46 [degC]
             l: 0.0002 [m]
             D: 2.872e-17 [m**2/s]
              = Dpiringer(wPET,<migrant: Octadecy [...] opionate - M=530.9 g/mol>,T=50.0 degC)
             k: 3.431 [a.u.]
              = kFHP(<ethylene terephthalate>,<migrant: Octadecy [...] opionate - M=530.9 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]
    -- [ layer 3 of 4 ] ---------- barrier rank=4 --------------
          name: "layer in PS"
          type: "polymer"
      material: "polystyrene"
          code: "PS"
       crystal: 0
            Tg: 100 [degC]
             l: 0.0001 [m]
             D: 4.787e-31 [m**2/s]
              = Dwelle(PS,<migrant: Octadecy [...] opionate - M=530.9 g/mol>,T=50.0 degC)
             k: 2.23 [a.u.]
              = kFHP(<styrene>,<migrant: Octadecy [...] opionate - M=530.9 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]
    -- [ layer 4 of 4 ] ---------- barrier rank=2 --------------
          name: "layer in PP"
          type: "polymer"
      material: "isotactic polypropylene"
          code: "PP"
       crystal: 0.5
            Tg: 0 [degC]
             l: 0.0003 [m]
             D: 2.333e-14 [m**2/s]
              = Dpiringer(PP,<migrant: Octadecy [...] opionate - M=530.9 g/mol>,T=50.0 degC)
             k: 17.55 [a.u.]
              = kFHP(<propylene>,<migrant: Octadecy [...] opionate - M=530.9 g/mol>)
            C0: 1000 [a.u.]
             T: 50 [degC]



    <multilayer with 4 layers: ['layer in gPET', 'layer in wPET', 'layer in PS', 'layer in PP']>



## Adapt automatically $k$ value for food

> We choose here ethanol as food simulant


```python
from patankar.food import ethanol
contactmedium = ethanol(T=50)
repr(contactmedium)
```

    Food object "ethanol" (ethanol = from pure ethanol down to ethanol 95%) with properties:
            CF0: 0 [a.u.]
    contacttime: 864000 [s]
        density: 1000 [kilogram / meter ** 3]
              h: 0.0001 [meter / second]
    surfacearea: 0.06 [meter ** 2]
         volume: 0.001 [meter ** 3]





    '<ethanol: ethanol>'



### Add `m[0]` to `contactmedium`
> This step makes appear the field `k0`
> > `  k0: 0.17613565 [a.u.]`
> >
> > `     = kFHP(<ethanol>,<migrant: phenyl-methane - M=92.14 g/mol>)`


```python
contactmedium.update(migrant=m[0])
```

    Food object "ethanol" (ethanol = from pure ethanol down to ethanol 95%) with properties:
            CF0: 0 [a.u.]
    contacttime: 864000 [s]
        density: 1000 [kilogram / meter ** 3]
              h: 0.0001 [meter / second]
             k0: 0.17613565 [a.u.]
               = kFHP(<ethanol>,<migrant: phenyl-methane - M=92.14 g/mol>)
    surfacearea: 0.06 [meter ** 2]
         volume: 0.001 [meter ** 3]



    <ethanol: ethanol>



print('k0 for ',m[0].compound," = ",contactmedium.k0)


```python
print('k0 for ',m[0].compound," = ",contactmedium.k0)
```

    k0 for  toluene  =  [0.17613565]


### We repeat the process for `m[1]`, `m[2]` and `m[3]` 


```python
print('k0 for ',m[1].compound," = ",contactmedium.update(migrant=m[1]).k0)
print('k0 for ',m[2].compound," = ",contactmedium.update(migrant=m[2]).k0)
print('k0 for ',m[3].compound," = ",contactmedium.update(migrant=m[3]).k0)
```

    k0 for  limonene  =  [0.55908817]
    k0 for  BHT  =  [0.25828644]
    k0 for  Irganox 1076  =  [0.63893733]

## Last words

Having a computational engine that choose properties and assumptions instead of you can be really annoying. You have several options to circumvent the rules:

-  **best option**: use `layerLink` objects (from `patankar.layer`) to enforce any $D$ and $k$ values for a specific layer (they are more general and manage also applicable to $C_0$, $l$, $T$…).

  > this methodology offers fine tunning at the level of each layer (one layer can use a model not the other).

  ```python
  D = layerLink("D", indices=0, values=some value)  # Create link for D
  object.Dlink = D  # Attach links to object
  ```

- use a generic `layer` instance (from module`pantankar.layer`) instead of a recognized polymer (LDPE, PP, PS, PET…)

- use a generic `foodlayer` instance (from module `patankar.food`) instead of a recognized food simulant (ethanol, water…)

- do not incorporate a migrant/substance into the `layer` and `food layer` 

- **easy option**: define a user `Dmodel` (it is initialized to “default”) or `kmodel` (it is initialized to “default”). For example, a dummy model returning `None` for any input will force the computational engine to use the values of $D$ and $k$ set with `myobject.D=...` or `myobject.k=...`

  ```python
   myobject.Dmodel=lambda *args, **kwargs: None # to force values set with myobject.D=...
   myobject.kmodel=lambda *args, **kwargs: None # to force values set with myobject.k=...
  ```

- **less recommended**: modify the model precedence rules used by the class `migrant` from `patankar.loadpubchem`.

- a bit more complicated: add your own models in `patankar.property` (follows instructions and examples inside)



***

<div style="border: 2px solid #4CAF50; border-radius: 8px; padding: 10px; background: linear-gradient(to right, #4CAF50, #FF4D4D); color: white; text-align: center; font-weight: bold;">
  <span style="font-size: 20px;">🍏⏩🍎 <strong>SFPPy for Food Contact Compliance and Risk Assessment</strong></span><br>
  Contact <a href="mailto:olivier.vitrac@gmail.com" style="color: #fff; text-decoration: underline;">Olivier Vitrac</a> for questions |
  <a href="https://github.com/ovitrac/SFPPy" style="color: #fff; text-decoration: underline;">Website</a> |
  <a href="https://ovitrac.github.io/SFPPy/" style="color: #fff; text-decoration: underline;">Documentation</a>
</div>
