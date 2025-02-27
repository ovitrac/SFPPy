<div style="padding: 20px; background-color: white; color: #C6392E; font-family: 'Georgia', serif; 
            font-size: 32px; font-weight: bold; text-align: center; 
            border: 3px solid #312289; border-radius: 10px; 
            box-shadow: 4px 4px 16px rgba(0, 0, 0, 0.3); line-height: 1.2">
  Risk Assessment of Migration From Packaging Materials Into Food
</div>



Yan Zhu^a^, Phuong-Mai Nguyen^b^, Olivier Vitrac^a^&dagger;

*^a^INRA, French National Institute of Agricultural Research, UMR 1145 Ingénierie Procédés Aliments, AgroParisTech, INRA, Université Paris-Saclay, 91300, Massy, France*

*^b^LNE, French National Laboratory of Metrology and Testing, 29 Avenue Roger Hennequin, 78197 Trappes Cedex, France*

 &dagger;[corresponding author](mailto:olivier.vitrac@agroparistech.fr)



# Abstract

This chapter reviews the principles to assess the risk of migration from food packaging materials and related application of such materials. The point of view considered in this chapter is that most of the past crises involving materials could have been avoided with a correct appraisal of mass transfer and the use of numerical tools. Migration modeling is broadly available today and accepted in major global regulations. This chapter does not advertise any platform but shows how their basic or advanced use can prevent the risk of migration regardless of existing or preexisting regulations. Risk assessment is proposed as a generic tool to revise the formulation of materials (thermoplastics, varnishes, coatings, lacquers, inks…), to optimize packaging design and supply chains, to promote new molecular design of additives, to devise improved recycling procedures, to accelerate the acceptation of new polymers and materials (*e.g.*, biosourced, biodegradable, recycled, with repeated use).



<small>Article Access and Citation Details</small>

<table style="width:100%; border-collapse: collapse; border: 1px solid #ddd; font-family: Arial, sans-serif;">
  <thead>
    <tr style="background-color: #C6392E; color: white; text-align: left;">
      <th style="padding: 10px; border: 1px solid #ddd;">Resource</th>
      <th style="padding: 10px; border: 1px solid #ddd;">Link</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 10px; border: 1px solid #ddd;">Original Reference</td>
      <td style="padding: 10px; border: 1px solid #ddd;">
        <a href="https://www.sciencedirect.com/science/article/pii/B9780081005965225018" target="_blank">Elsevier Record</a>
      </td>
    </tr>
    <tr>
      <td style="padding: 10px; border: 1px solid #ddd;">Converted HTML (Comfortable Layout)</td>
      <td style="padding: 10px; border: 1px solid #ddd;">
        <a href="MigrationModeling.html" target="_blank">HTML with original layout</a>
      </td>
    </tr>
    <tr>
      <td style="padding: 10px; border: 1px solid #ddd;">PDF to HTML Conversion</td>
      <td style="padding: 10px; border: 1px solid #ddd;">
        <a href="MigrationModeling.pdf" target="_blank">PDF</a>
      </td>
    </tr>
  </tbody>
</table>

<p style="margin-top: 15px; font-family: Arial, sans-serif; font-size: 14px;">
  <strong>Citation:</strong><br>
  Yan Zhu, Phuong-Mai Nguyen, Olivier Vitrac, <i>Risk Assessment of Migration From Packaging Materials Into Food</i>, 
  Reference Module in Food Science, Elsevier, 2019, ISBN 9780081005965, 
  <a href="https://doi.org/10.1016/B978-0-08-100596-5.22501-8" target="_blank">https://doi.org/10.1016/B978-0-08-100596-5.22501-8</a>.
</p>


# 1     Content

[TOC]

# 2     Introduction

**What is migration?**

Migration is a general term for spontaneous mass transfer of chemical substances, and, in the context of food packaging, it indicates an extraction of packaging constituents and their transfer to the food. Industry and authorities recognized the fate of cross mass transfer between materials and the resulting contamination of the food evenly. The term migration was consequently preferred to the contamination one in the scientific literature and legal documents. During the last decades, the concern about the safety of food contact materials (FCM) raised with our appetite for transformed food and ready to eat meals, and with our always growing needs for disposable packaging. Nowadays, FCM is identified as the prevalent source of exogenous chemical contaminants in food, ahead of pesticides, veterinary drug residues and other environmental contaminants . The ubiquitous contamination issue was thought to be restrained initially to contact layers and materials, but it is far from being the case. Modern food packaging systems are, indeed, printed, coated, laminated and associated with other materials. The whole history of these materials must be considered, as they may have been subjected to repeated use, brought to a second life via a recycling process, stored in reels or stacks contacting internal and external surfaces; shipped with other materials during long and warm periods. The whole contamination problem can be envisioned as a cross mass transfer of several substances between Matryoshka or nesting dolls. The pre-weighted food feeds the smallest doll and is surrounded by many layers including the rigid walls of the primary packaging, a sleeve, the transport cardboard box containing several sale units, the treated wood pallet, a wrap cling film, etc. until the freight container. In ready to eat meals and convenience foods, additional components may be present internally such as a bag preventing direct contact with walls, individual packages or wraps for portion control, separators and specific holders, sachets for seasonings, active elements to increase food shelf-life, etc. The whole picture is not complete without citing the many tie layers, the glued labels, the printed and coated surfaces. The ultravacuum and aerospace industry would have regarded similar materials and combinations as a substantial reservoir of organic compounds without exception. As an illustration, the NASA compiled more than 35,000 outgassing data [1](#_ENREF_1),  for almost any material, which could enter in a spacecraft, including many commodity materials such as thermoplastics, coatings, and tapes. 

 

**Our approach to risk assessment**

This chapter provides a comprehensive description of molecular processes responsible for the migration of packaging constituents and their pathways to the food with or without direct contact. The recurring structure of contamination scenarios leading to the past major crises is discussed and analyzed regardless of the modalities of their regulations. The distance between the facts and practices is maintained all along this chapter, as the regulations and the good manufacturing practices follow the crises and rarely precede them. Major regulations (the US, European and Chinese) are sketched to highlight their convergence on the use of modeling to demonstrate compliance and to evaluate the safety of recycled materials. The miscalculation of the connection between chemical structure and physicochemical properties (volatility, solubility, diffusivity) has been the foremost cause of past crises. Computer and proper simulation procedures can assist efficiently small and intermediate industries in overcoming internal knowledge limitations on materials, mass transfer, and physical-chemistry. By comparing with acceptable thresholds, migration modeling can be extended at low cost to non-evaluated and non-intentionally added substances. In the foreseeable future, similar techniques might be used to tackle the diffuse risks raised by endocrine-disrupting chemicals [3](#_ENREF_3), [4](#_ENREF_4) alone or in mixing cocktails [5](#_ENREF_5). More globally, the extension of predictive tools and approaches will benefit not only the evaluation of the contribution of food contact materials to the global [exposome](https://en.wikipedia.org/wiki/Exposome) [6](#_ENREF_6), but it will also facilitate the adoption of preventive approaches all along the supply chain.  Safe-by-design approaches, including additive redesign, optimization of the formulations (choice of substances and amounts), new packaging design and good manufacturing practices, will reduce the risk of unintended food packaging-interactions [7](#_ENREF_7), [8](#_ENREF_8). Improving the way food ingredients are stored and processed will bring additional risk reduction, beyond the reduction of the migration in the finished product.

 

**The scope of migration modeling**

The scope of migration modeling has been underestimated in the past and limited to compliance testing under worst-case scenarios. The uncertainty and the pioneered methods were too coarse thirty years ago. We suggest that the scope can be broadly extended today as shown in [Figure 1](#figure1). Earlier models could cover only single materials, simple geometries without any dynamic change of conditions. The most advanced models can today incorporate information at molecular scale and cover an entire supply chain. The chapter focuses on the key details and features required to get robust modeling of migration at the scale of a material, component (label, cap…), of an entire food packaging, of industrial practices. The methodology to get estimates of consumer exposure are not covered in this chapter because they are directly related to the design of the packaging itself. 

![fig1](./src/fig1.png)   

<a name="figure1"></a>**Figure 1**. Evolution of the scope of migration modeling during the last decades



**Beyond migration modeling: a holistic science**

Once used, the packaging losses all its value, but its cost remains high for the environment. Revised versions of life-cycle-analysis propose to add chronic effects on health due to long-term impacts on the environment[9-11](#_ENREF_9). We promote this initiative for applications in food contact but also for cosmetics, pharmaceutical, medical, biotechnological ones and for any situation where the release of chemicals by materials is of concern (*e.g.*, occupational exposure). Innovations could come shortly at the cost of a full revisitation of our engineering practices. A global science needs to emerge beyond common mistaken beliefs: migration can be avoided, biosourced materials are safer, biodegradable as recycled materials can be used without a safety assessment. We do not pretend to cope with all problems, but we suggest thoughtful routes to evaluate the risks of contamination at any stage of the supply chain and consequently to measure the benefits of alternative solutions. From step to step, we could approach better inertia as it has been achieved to minimize outgassing in ultravacuum systems.



# 3     What did we learn from crises?



------
**Highlights on crises associated with food contact materials**

- The contamination of food by materials in contact is never fortuitous, but it may be not avoidable.
- Only the nature of the migrating substances, the extent of the migration can be controlled by our choices.
- Most of all crises if not all could have been predicted with relatively simple descriptions of mass transfer phenomena and thermodynamics.
- As a corollary, consumer exposure to substances from food contact materials could be reduced.

------



## 3.1    A short history

Crises tend to predate regulations and frequently them — we coin crisis a situation where food safety is seriously questioned due to systematic contamination by one or several FCM. The contamination was usually unexpected but not unforeseeable due to its anthropogenic nature. Comparatively to food infection and food intoxication, the possibility of crises by FCM substances has been recognized lately in Europe. The obligation of traceability of all packaging components to organize the recall of packaged food products were implemented only in 2004 through the framework regulation [2035/2004/EC](https://eur-lex.europa.eu/legal-content/EN/TXT/PDF/?uri=CELEX:32004R1935&from=IT). 

Thirty years ago, western governments endorsed enthusiastically an early idea hypothesized by Jerome Nriagu [12](#_ENREF_12) and popularized by Clair Patterson [13](#_ENREF_13), whereby Roman civilization collapsed as a result of lead poisoning.  “Although today lead is no longer seen as the prime culprit of Rome’s demise” (quoting [14](#_ENREF_14)), lead poisoning from leaded pottery and earthenware was known for centuries. All possible sources of lead were individually tracked by authorities in the nineteenth century, in particular, after the development of wrought-iron canisters. The French ordinance of March 21th, 1879 (see p 231 of [15](#_ENREF_15)) prohibited, for instance, the use alloys of tin and lead for all inner parts including welding.  The French regulation of 1908 was even more explicit “no food substance should contain any harmful product or chemical substance” [16](#_ENREF_16). In modern times, the editorial of the New England Journal of Medicine was heading in 1972 “[The invisible pollution](https://www.nejm.org/doi/full/10.1056/NEJM197211302872213)” [17](#_ENREF_17) after the discovery of plasticizers in human blood stored in polyvinyl chloride (PVC) bags [18](#_ENREF_18), [19](#_ENREF_19). Similar contaminations were associated with tubbing used for culture tissues [20](#_ENREF_20). The exact nature of the contamination mechanism was not fully established at the time, and an analogy with the corrosion of metallic materials was falsely suggested [21](#_ENREF_21). The first mechanistic review of what called “extractivity” or “migration propensity” appeared only in 1980 [22](#_ENREF_22), [23](#_ENREF_23). Three mechanisms of contamination were listed:

- contamination only from the extreme surface of the material in contact;

- the diffusion-controlled release of the material in contact;

- penetration of the polymer matrix by the contacting phase (liquid) and subsequent extraction of materials constituents.


Early descriptions were strongly influenced by the behavior and the dominance of PVC in the seventies and eighties, and by the lack of sensitivity of contemporaneous analytical techniques. In spite of erratic results and difficulties in getting reproducible kinetics, the corollary reasoning supported the confidence of stakeholders in the apparent inertia of thermoplastics and thermosets lastingly. It was falsely thought that:

- the absence of direct or permanent contact, 

- aqueous contact,

- high molecular weight additive or residues,

- low temperatures


would prevent any significant migration and do not need proper attention. At the time, only a study using radio-labeled additives [24-27](#_ENREF_24) carried out under contract for the Food and Drug Administration and, subsequently, interpreted in detail during the Ph.D. thesis of Thomas P. Gandek at MIT [28](#_ENREF_28) highlighted several abnormalities, which were anticipating future crises. In poor barrier polymers, the hydrodynamic conditions in the contacting liquid were showed to control the release; but neglecting it them was not underestimating migration, on the contrary. In aqueous food simulants and presumably in any aqueous-type food, the decomposition of additives displaces continuous the apparent thermodynamical equilibrium between the material and the liquid in contact. Contrary to previous descriptions, the contamination was found unbounded [29](#_ENREF_29), [30](#_ENREF_30).

Similar conclusions were found for large additives migrating to dry food simulants for short periods at temperatures elevated but sufficiently close to those met during transportation [27](#_ENREF_27). The most outstanding finding was that the migration rate could not be overestimated by experiments using corn oil. Accelerating testing using food-simulating liquids is a common practice to evaluate the risk of migration, but it was emphasized that more rationale was required for evaluating with sufficient confidence the risk of migration for new polymers. Simulating liquids should be chosen respectively to the nature of both the migrating substance and of the original polymer. Fatty food simulants offer worst-case migration and extraction capabilities only for hydrophobic substances in apolar polymers. Aqueous simulants are more aggressive for polar or charged substances, and polar polymers [31](#_ENREF_31). Since most of the food products are multicomponent and multiphasic (*e.g.*, emulsions, gels, cake with chocolate, etc.), they cannot be reduced easily to a single contact phase when different classes of migrants are involved.

 

## 3.2    Why thermodynamics counts

Thermodynamics has been praised by the whole packaging community, including the chemical, compounding, processing, recycling and food industries, as well as authorities and safety agencies. It has been regularly as the primary argument to justify the conditions of compliance testing (choice of simulant, test temperature and contact time, extrapolation rules and migration calculations) and to authorize recycling process of polymers, active packaging, etc. A naïve reasoning may, however, lead to severe consequences, which should be underestimated. A common mistake is to assume that mass transfer stops after some long time. The transferred amount is assumed to reach, indeed, a maximum controlled by the partition coefficient between the packaging and the food. Statistic mechanics teaches us that this macroscopic description is oversimplified and proceeds with an analogy between a mechanical equilibrium and a chemical equilibrium. At molecular scale, all the substances continue to move freely at equilibrium as they were moving before the whole packaging-food system reaches a macroscopic equilibrium. In a closed system, the equilibrium is associated to a zero net mass balance across the packaging-food interface: the number of molecules of type A entering in the food is exactly compensated by the number of molecules of type A leaving the food. We might think that because the food has a larger volume than the packaging, the return of molecules is unlikely. It is however not correct because the substances are transported faster in food than in dense polymer matrices. Only when the volume ratio food-to-packaging becomes very large, the probability of return approaches zero and a total extraction of A is expected regardless of its affinity for the food.

When a chemical reaction transforms the species A into the species B in the food, the previous balance is profoundly modified, and more substances A are invading the food than substances A returning to the packaging. Similarly, when the substance A cleaves into two breakdown products $a1$ and $a2$ in the packaging, parts of substances A that migrated in the food are partly reabsorbed by the packaging itself. When the reaction is almost complete, only traces can be identified in the food. On the opposite, substances $a1$ and $a2$ can migrate (faster than A) and may remain undetected if not specifically targeted.

Adding more packaging components (several layers, cap, label, etc.), variable temperatures and complex contact conditions complicate the mass transfer description, but thermodynamics always provides the relationships to encompass all possibilities of exchanges including in multiphasic foods. The thermodynamics needed is not equilibrium thermodynamics, as it ignores the time-course of the migration processes, but a local version, where equilibrium is reached only at the interface between each phase and component. The classical non-equilibrium thermodynamical concept of local thermodynamic equilibrium is robust enough to integrate non-linear sorption isotherms, coupled mass transfer and plasticizing effects. When the relaxation of polymer systems is longer than the timescale of mass transfer (*e.g.* , in glassy materials, hysteresis effects), constitutive equations need to be modified to integrate the mechanical behavior of the materials (swelling or densification).

No food packaging system is enough isolated from the rest of the world so that a true thermodynamical equilibrium cannot be finally observed literally in real-life systems. Preventing the loss with surroundings in tests and calculations maximizes the amount transferred to the food. Conversely, not considering the possibility of redistribution of migrants between materials during their lifetimes hampers the proper management of non-intentionally added substances. 

In shorts, most of the experimental evidence supporting regulation and rules were obtained on simple materials, preferably monolayer and apolar ones at rubber state. Equilibrium was reached rapidly in conditions accelerated comparatively to the real shelf-life of food products. In this case and only in this case, the concentration in the simulating liquid increases monotonously with time.

Thermodynamics offer a robust framework to address all previous issues at the scales of molecules, where the interactions and the macroscopic properties emerge from the vibration of atoms. Theories or robust inference rules have been developed to relate the chemical structure of the migrants and the polymers to the diffusion and partition coefficients, without requiring an explicit description of molecular interactions. The principles of quantitative structure-relationships detailed in this chapter are sketched in [Figure 2](#figure2). They should not be considered as definitive rules, but as an ongoing process, where updates are regularly obtained.

<img src="./src/fig2.png" alt="fig2" style="zoom:80%;" />

<a name="figure2"></a>**Figure 2**. Principles of quantitative structure relationships to predict common transport and thermodynamic properties needed for migration modeling.

 

## 3.3    Modern crises

 “*Early civilizations adopted laws that punished sellers of tainted food*” as quoted by Merill [32](#_ENREF_32), In 1958, the US Food, Drug and Cosmetic Act introduced the Delauney clause, which prohibits the use of any substance as food additive “if it is found to induce a cancer when ingestion by man or animal” [33](#_ENREF_33). In the US regulation system, the concept of additives is broad and comprises substances, which may become a component of food or otherwise may affect the food characteristics [34](#_ENREF_34). They include therefore any substance released by food packaging, regardless of the nature of the materials. As a result, it is the responsibility of the industry to submit a dossier to the Food and Drug Administration (FDA) to get a new substance, a new polymer and new application of packaging approved. [Similar rules were enforced successively in EU](https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=LEGISSUM:l21301&from=FR) via the directive [1990/128/EEC](https://eur-lex.europa.eu/legal-content/EN/TXT/PDF/?uri=CELEX:31990L0128&from=EN) [35](#_ENREF_35), the regulations [2002/72/EC](https://eur-lex.europa.eu/legal-content/EN/TXT/PDF/?uri=CELEX:32002L0072&from=EN) [36](#_ENREF_36) and [10/2011/EC](https://eur-lex.europa.eu/legal-content/EN/TXT/PDF/?uri=CELEX:32011R0010&from=fr) [37](#_ENREF_37). The dossiers apply, however, only to initial substances of plastic materials monomers and additives. China adopted recently legal requirements close to the European system for plastics with a premarket approval for both the plastic resins and the additives [38](#_ENREF_38). The inventories of substances are compared in [Table 1](#table1).

<a name="table1"></a>**Table 1**. Comparison of the Inventories of food contact substances in US, EU and China. Sources [37-39](#_ENREF_37) 

| **MATERIALS**                                          | US regulation[39](#_ENREF_39)                                | EU regulation and provision                                  | Chinese regulation[38](#_ENREF_38)                   |
| ------------------------------------------------------ | :----------------------------------------------------------- | ------------------------------------------------------------ | ---------------------------------------------------- |
| **Plastics,   resins, additives, polymerization aids** | 1340 food contact   notifications + GRAS and prior sanctioned substances         (Title 21 CFR Parts 175, 177, 178) | 843 in the EU positive list[37](#_ENREF_37) including   428 with SML or SML   group, and   587 additives | 339 in positive list   86 with SML/QM or   SML group |
| **Rubber**                                             | 176 substances in the French positive   list [40](#_ENREF_40) | 88 substances   23 with SML/QM                               |                                                      |
| **Printing   inks**                                    | 5104 substances   (IAS) in EuPIA guidance[41](#_ENREF_41) and Swiss Ordinance[42](#_ENREF_42). For NIAS, see[43](#_ENREF_43) | 97 substances   33 with SML/QM                               |                                                      |
| **Paper and board (including recycled)**               | 482 substances   (Title 21 CFR Part 176)                     | 1556 substances additives in Council of   Europe resolutions [44](#_ENREF_44) and good manufacturing practices [45](#_ENREF_45). For recycled materials, see [46](#_ENREF_46), [47](#_ENREF_47) | 277 substances   77 with SML/QM                      |
| **Compliance   testing by modeling**                   | Yes (plastics)                                               | Yes (restricted to plastics)                                 | No restriction on   applicable materials             |
| **Risk Assessment including migration   modeling**     | Broad range of applications in   petitions in relation with consumer exposure determination | Broad range of applications in   petitions in relation with an upper estimation of consumer exposure | Not applicable in petitions                          |

*GRAS = generally recognized as safe; SML: specific migration limit (maximum concentration in food); QM: quantity maximum (maximum amount in the material before contact); IAS: intentionally added substances; NIAS: non-intentionally added substances.*



 A crisis occurs when a substance not authorized is found in the packaging material or when the levels in food raise concern due to its ubiquitous distribution or due to significant exposure of specific or global populations. As an illustration, the substances from materials intended to be in contact with food and exceeding EU migration tolerances during the period 2002-2018 are listed in Figure 3 in decreasing order of occurrence. Half of the 1956 cases are associated with imported tableware and with contaminations by heavy metals. Packaging alone represent less than 6% of the total with lids involved in 45% of cases. The main crises associated with packaging materials and their likely causes and consequences are reviewed in [Table 2](#table2).

 

<a name="table2"></a>**Table 2**. Analysis of the major crisis involving food contact materials

| **Crisis   (class)**                                         | **Period   (key references)**                   | **Cause**                                                    | **Consequence**                                              |
| ------------------------------------------------------------ | ----------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| **Levels   of plasticizers from PVC cling films above 100 mg****⋅****kg-1** | 198x   [48-53](#_ENREF_48)                      | Widely   used as additives (plasticizer, solvent) with large amounts, not covalently   bonded to the material backbone | Restrictions or   ban of many phthalates in FCM applications |
| **High levels BADGE and of its reaction   products of BADGE from can coatings** | 1997-today [54](#_ENREF_54)                     | starting substance for the manufacture of epoxy   resins, used as an additive, functioning as a stabilizer and as a   plasticizer à reaction with medium in contact | TDI   0.15 mg/kg b.w (BADGE)   SML   1 mg/kg of food (reaction products) |
| **Primary   aromatic amines from agglomerated cork stoppers** | 2000-today   [55](#_ENREF_55), [56](#_ENREF_56) | Surface   treatment, adhesive, lubricant à   reaction products | SML 0.01 mg:kg of food applied to the sum of PAA   released (annex II of [37](#_ENREF_37)) |
| **Primary aromatic amines from laminates**                   | 2000-today [57-59](#_ENREF_57)                  | polyurethane adhesive à reaction due to thermal treatment    | SML   0.01 mg/kg of food applied to the sum of PAA released (annex II of [37](#_ENREF_37)) |
| **High   levels of epoxidized soybean oils**                 | 1988   – today [60](#_ENREF_60)                 | Large   amounts in PVC plasticized, high lipid solubility, long time storage | Lowering of SML   to 30 mg/kg for infant food                |
| **High levels Nonyl-phenols**                                | 2004-today [61](#_ENREF_61)                     | Breakdown product of tris(nonylphenyl)phosphite   (authorized in EU 10/2011) | strictly   limited by Directive 2003/53/EC [62](#_ENREF_62)   restriction   REACH annex XVII [63](#_ENREF_63) |
| **Semicarbazide   leached from the thermal decomposition of azodicarboamide   from gaskets  used in baby food jar closure technology (press twist   and twist-off lids)** | 2003-2005   [64](#_ENREF_64), [65](#_ENREF_65)  | Breakdown product of azodicarbonamide during the   heat treatment | ban [66](#_ENREF_66)                                         |
| **ITX**                                                      | 2005-today [67](#_ENREF_67)                     | Bad identification of transfer/contamination paths   from ink to foodstuffs   Bad information exchange between stakeholders in supply chain | No   regulated in plastic material regulation   REACH   study |
| **Benzophenones   from printing inks**                       | 1995-today   [68](#_ENREF_68)                   | Bad   identification of transfer/contamination paths   Bad   information exchange between stakeholders in supply   chain | SML: 0.6 mg/kg                                               |
| **Bisphenol A leached from Baby Bottles**                    | 2004-today [69](#_ENREF_69), [70](#_ENREF_70)   | Widely used in infant products                               | Lowering   of SML to 0.05 mg/kg   Replaced   by other substances   New   risk assessment ongoing by EFSA |
| **Ubiquitous   contamination by mineral oils from recycled papers and boards** | 1993-today   [71](#_ENREF_71), [72](#_ENREF_72) | Recycled   paperboard, difficult in analytical analysis      | Proposition of   SML of 0.5 mg/kg                            |

 

The inset of [Figure 3](#figure3) and [Table 2](#table2) would suggest that the number of contamination cases would increase with time. The better effort of authorities to identify, track and contain the contamination by FCM led to more frequent reporting. Regardless of the real risks, the successive crises can coatings, [bisphenol A](https://en.wikipedia.org/wiki/Bisphenol_A) and [mineral oils](https://en.wikipedia.org/wiki/Mineral_oil) increased the awareness of the general population and, in return, promoted the sake of different management strategies and better integration along the supply chain.

​                                                  

![fig3](./src/fig3.png)

<a name="figure3"></a>**Figure 3**. List of contaminants from food contact materials reported in the European Rapid Alert System for Food and Feed[73](#_ENREF_73) (extraction on Dec 31st, 2018). The inset shows the evolution with time of the number of alerts issued by all member states, corresponding during the period 2002-2018 to 955 border rejections, 459 alerts, 288 information for attention, 284 information and 216 information for follow-up.



# 4     Risk assessment for decision making



------

**Highlights on modeling, risk assessment, and decision making**

- Migration modeling is a cognitive process aiming at capturing the essential mass transfer phenomena responsible for the contamination of food by substances originating from materials.
- Calculations are carried out in a way that they guarantee they are more severe than real conditions.
- Only compliance can be demonstrated by calculations and simulations.
- Calculations are in essence different than accelerated tests and should be used to reproduce real but conservative contact conditions (time, temperature, interactions with food).
- Sophistications in calculations require to be introduced progressively and in a comprehensive manner starting from the most conservative (severe) conditions.
- Variability and uncertainty must be characterized and documented.
- The whole process can be integrated into preventive approaches (safe-by-design) or stochastic calculations (see section 7)

------




## 4.1    Overview of migration modeling

### 4.1.1     What is currently permitted?

Notwithstanding the introduction of accelerated tests for compliance testing, simplified alternatives to traditional migration testing were sought in the nineties. The common interests of the industry and the authorities were summarized by Begley [74](#_ENREF_74), “Traditionally, migration tests are performed by using food-simulating liquids such as water, edible oils, ethanol/water solutions and sometimes food. These tests are time-consuming in two ways; generally, the accelerated tests run for at least ten days, and the analysis of the migrants at low concentrations in the simulants or food is generally difficult. These analyses are also expensive and generate hazardous laboratory waste”. Migration modeling was proposed both in the US (earlier works [75](#_ENREF_75) ) and in EU (earlier works [76-79](#_ENREF_76)) to tailor the process of migration assessment. Migration modeling is nowadays mostly extended in EU via a specific task force TF-MATHMOD publishing updated guidance [80](#_ENREF_80) Migration modeling is broadly accepted for compliance testing [37](#_ENREF_37) and risk assessment (see [81](#_ENREF_81)) of food contact materials with the following strict restrictions: 

1.  estimated values – whatever the calculation procedure and underlying assumptions – must be at least as severe as the real test and overestimate the migration;
2. calculations cannot be used to demonstrate the non-compliance.

Similar methodologies are used by the industry to calculate the level of decontamination (cleaning) of materials in mechanical recycling processes. But in this case, the amount released by the material should be estimated accurately and not overestimated. Additionally, it is worth noticing that the state of polymers is very different between the conditions met by the consumer (moderated temperatures) and the industrial conditions of recycling (high temperature, solvent-swollen). Only mechanistic modeling can cover both cases. Modeling is not restricted to any material (thermoplastic, thermoset, paper, and board), but it has been tested chiefly for thermoplastics and with non-uniform coverage.

There is no limit to the scope of modeling in current regulations and practices. Compliance testing can be seen as the simplest usage of modeling and validating a closed loop of recycling as the most complex application (see Figure 1). Repeated use, composite materials and safe-by-design approaches can be envisioned of intermediate complexity. All applications (small and large) obey in fact to same global scheme, where the result obtained at the scale immediately below is used at the upper scale, as shown in [Figure 4](#figure4), The complexity of modeling relies on the number of scales considered and not on the system (film, bottle, etc.) subjected to the modeling activity. The keys to identifying the number of scales, components and steps are discussed all along the chapter.

<img src="./src/fig4.png" alt="fig4" style="zoom:80%;" />

<a name="figure4"></a>**Figure 4**. Overview of the nested modeling strategy to predict migration for all applications depicted in Figure 1.

 

### 4.1.2     Modeling using a tiered approach: from worst-case scenarios to detailed conservative ones

The exact value of the contamination of the food is never achievable because the conditions of contact are variable (time, temperature) between comparable products and because our knowledge of molecular mechanisms is perfectible. As a result, the practice seeks successive approximations of the migration in a tiered approach, as shown in [Figure 5](#figure5). At the first tier, the estimation is very coarse and connected with the highest overestimation factor. If the determined concentration at tier is higher than the threshold of concern, the next tier is triggered by introducing substantial refinements and details, and so on. The process stops when no additional information can be introduced (experiments need to be preferred) or when the calculated concentration is lower than the threshold of concern. The lowest tier within the threshold of concern defines the proper level of knowledge required to demonstrate compliance or to guarantee the safe use of a material, substance, or process. There is no systematic procedure to identify the minimum tier to reach the goal, and only the needed information can be listed.



<img src="./src/fig5.png" alt="fig5" style="zoom:80%;" />

<a name="figure5"></a>**Figure 5**. Principle of the tiered approach to demonstrate compliance for food contact materials. Compliance is demonstrated as soon as the estimated concentration is greater than the threshold of concern. Tier 1 is usually associated with total migration (see Eq. $\ref{eq1}$).



The possibilities and prerequisites for using modeling in compliance testing are reviewed in [Table 3](#table3). The mentioned tiers $R_1$, $R_2$ and $R_3$* are indicative and corresponded to a modeling level more sophisticated than the first tier. The first tier is usually coined to “total migration” and corresponds to the total transfer of substances into food (see Figure 13). The corresponding concentration in food, $C_{F}^{tier\ 1}$, is determined by a “dilution” model:

$$
C_{F}^{tier\ 1}={{L}_{P/F}}C_{P}^{0}\label{eq1}
$$

where $C_{P}^{0}$ is the initial concentration in the material (regardless its distribution) expressed in mass per volume (preferred in this work) or in mass per mass (industrial practice). ${{L}_{P/F}}$ is the material-to-food volume or mass ratio. When one-dimensional representations are used, ${{L}_{P/F}}$is also the ratio between the thickness of packaging walls, denoted ${{l}_{p}}$, and the characteristic dimension of the food, denoted ${{l}_{F}}=\frac{{{V}_{F}}}{A}$ , where $A$ is the effective surface area in contact, counting usually the total surface area in contact with the food and the headspace. By contrast, the food volume ${{V}_{F}}$ is restricted to the condensed part of the food (solid or liquid).



<a name="table3"></a>**Table 3**. Prerequisites to use calculations as an alternative to migration testing.

| **Prerequisites**                                            | **Type of estimate**                                         | **Examples of works**                                        | **tier†**                                         |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------- |
| **Migration modeling or related   calculations**             | lectures , reference text books (specific or general)        | Lectures on migration [82-84](#_ENREF_82),  text books on migration [85-90](#_ENREF_85) ; reviews and case studies on migration modeling [86](#_ENREF_86), [89](#_ENREF_89), [91-95](#_ENREF_91) [96-98](#_ENREF_96) ; reference text books on packaging[90](#_ENREF_90), [99](#_ENREF_99); reference text books on mass transfer[100-107](#_ENREF_100) | R1   R2       R2   R3   R3       R3               |
| **Identity   of material**                                   | technical   specifications,     recycling   code, measurements | supplier, regulation,  standard [108](#_ENREF_108), [109](#_ENREF_109), FTIR spectra | R1   R2   R3                                      |
| **Characteristics of the polymer**                           | density , glass transition temperature $T_g$, crystallinity  | supplier; handbooks[110](#_ENREF_110), [111](#_ENREF_111); measurements | R1   R2   R3                                      |
| **Identity** **of the substance**                            | real   substance , chemical   structure, chemical   descriptors | supplier; deformulation [112](#_ENREF_112) and/or spectroscopy   [113](#_ENREF_113), [114](#_ENREF_114);   analoguous   substance | R1   R2   R3                                      |
| **Packaging geometry**                                       | 1 kg packed in 6 dm^2^  , 1D approximation of real geometry, 3D real geometry | regulation; supplier, end-user; research work[111](#_ENREF_111), [115-118](#_ENREF_115) | R2   R2       R3                                  |
| **Contact conditions**   **(time, temperature,   phase in contact…)** | standard   test conditions ,   accelerated   conditions, real   conditions | end-user; regulation                                         | R1   R2                                           |
| **Initial concentration**                                    | real values, overestimates                                   | supplier; guidance, orientation formulabrute force deformulation [113](#_ENREF_113), [114](#_ENREF_114) | R1   R2   R3                                      |
| **Diffusion   coefficients (see §****5****)**                | real   values, overestimates,  molecular   theory, molecular   modeling | measurements,   literature [119-130](#_ENREF_119)  *Piringer*   model[131](#_ENREF_131) or equivalent [100](#_ENREF_100); free-volume   theories and their extensions[123](#_ENREF_123), [124](#_ENREF_124), [126](#_ENREF_126), [132-136](#_ENREF_132) ; Molecular Dynamics simulation [137-142](#_ENREF_137) | R1   R2   R3   R3       R3       R3               |
| **Partition coefficients or sorption isotherms**    **(see §****6****)** | real values, overestimates, molecular theory, molecular modeling (*Flory Approximation*),  molecular modeling | measurements, literature; guidance [143](#_ENREF_143); *Kirkwood-Buff* theory [144](#_ENREF_144) ; *Flory -Huggins* theory[145-148](#_ENREF_145); solubility parameters [110](#_ENREF_110), [112](#_ENREF_112), [149](#_ENREF_149), [150](#_ENREF_150) ; theory of interacting liquids[151](#_ENREF_151); thermodynamic integration, insertion method [152-157](#_ENREF_152) at atomistic scale; atomistic calculations within the *Flory  approximation*[112](#_ENREF_112), [113](#_ENREF_113), [158-164](#_ENREF_158) | R1   R2   R3   R3   R3           R3   R3       R3 |
| **Mass transfer   resistance in the contacting phase**   **(see §****4.2.3.4****)** | none, correlation   diffusion   coefficient, flow            | worst-case; liquid   or gas model [100](#_ENREF_100), [101](#_ENREF_101), [103-105](#_ENREF_103), [165-167](#_ENREF_165);  Graetz   type problem; full simulation coupling[168](#_ENREF_168) | R1   R2   R2   R3                                 |
| **Acceptable thresholds**   **(see §****8.4.1****)**         | $QM$[37-39](#_ENREF_37), $SML$[37](#_ENREF_37), [38](#_ENREF_38), [169-174](#_ENREF_169), *threshold of regulation*[175](#_ENREF_175) , $TTC$[176](#_ENREF_176) | Regulations; regulation policy [74](#_ENREF_74) ; recommendations (GMP), literature | R1   R2   R3                                      |

†Recommended for compliance testing at the first ($R_1$), second ($R_2$) and third ($R_3$) tier. If the value of the migration is larger at the second tier than at the first one, use the second tier. If the value of the migration is larger at the third tier than at the second one, use the third one. $SML$=**specific migration limit*; $QM$=*maximum amounts*; $TTC$= *threshold of toxicological concern*; $GMP$=*good manufacturing practices*. $QM$ and maximum concentration in food $C_{F}^{\max }$ are related as $C_{F}^{\max }={{\left( \frac{{{V}_{F}}}{{{V}_{P}}}+\frac{1}{{{K}_{F/P}}} \right)}^{-1}}{{\rho }_{P}}QM$, with ${{V}_{F}}$ and ${{V}_{P}}$ the food and packaging volume respectively, ${{\rho }_{P}}$the density of the polymer and ${{K}_{F/P}}$ the partition coefficient.

 

### 4.1.3     Key steps in migration modeling and risk assessment approaches

Any migration modeling for compliance testing, risk assessment, and safe-by-design approaches should be initiated by the review of five important sections with an intent of providing an inventory on:

1. the formulation of materials (intentionally-added or not substances),

2. the components included in the design;

3. the steps followed by the material, the finished packaging, and the packaged food;

4. the information obtained from suppliers, regulations, industrial recommendations;
5. the described mechanisms of contamination.

For each section, the items need to be ranked and prioritized according to their suspected or foreseen importance on the contamination of the packaged food, as shown in [Figure 6](#figure6). The principles and illustrations described in this chapter can be used to extend the systems, steps, substances,… under scrutiny beyond primary food packaging and contact layers. 

![fig6](./src/fig6.png)

<a name="figure6"></a>**Figure 6**. Generic steps to review in migration modeling and safe-by design approaches. The depicted example corresponds to the review for a new aseptic carton packaging for milk to be consumed by infants.

 

## 4.2    Migration modeling for compliance testing and beyond

In agreement with EU recommendations [80](#_ENREF_80), this section details the assumptions and conditions suitable for compliance testing at high tiers using a comprehensive description of mass transfer. The levels of description would correspond to tiers $R_2$ and $R_3$ in [Table 3](#table3).

### 4.2.1     Principles of migration modeling
------

**Highlights on the principles of migration modeling**.

- Migration obey simply speaking to the well-known laws of diffusion.

- The mass transfer from one material to another material, or the food requires specific treatment and attention as it is not implemented by default generic numerical software (commercial or not).

- One-dimensional mass transfer calculations are sufficient for most of the applications if mass balances are well preserved.

- The cost of modeling is dramatically reduced by simulating the contacting phase implicitly with proper boundary conditions. It is important to note that the substances within any eventual mass transfer boundary layer are not included in the food mass balance when implicit models are used.

- Beyond its obvious numerical advantages (abacuses, master curves, pre-tabulated results), dimensionless formulations based on similitude principles facilitate the review of model assumptions and results.

------

Substances non-covalently bound to the polymer are subjected to thermal agitation, which causes in return a random translation of their center-of-mass. Each additive, monomer, residue jumps or walks randomly in the polymer matrix from one accessible void to the next one. The substances reach the interface eventually with the food, where the same hopping process is repeated usually at a higher pace. When volatile substances meet a gas phase, their skew trajectories are governed by the collision with gas molecules. In all cases, random walks occur in three-dimensions, but a concentration gradient develops only at leaching interfaces, in its perpendicular direction. As the walls of the packaging are thin compared to the characteristic food dimension, migration can be approximated as a one-dimensional mass transfer problem as shown in Figure 7. The migrating substances are depicted either as individual molecules or solutes (*i.e.*, scales are not respected) showing microscopic concentration fluctuations or as smooth macroscopic concentration profiles (continuous lines). Different symbols are used whether the solutes are in the polymer (position $0\le {x}/{{{l}_{p}}}\;<1$ ); in a small boundary layer in the food and of thickness ${{x}_{BL}}$, where a concentration gradient exists (position $1<{x}/{{{l}_{p}}}\;\le {1+{{x}_{BL}}}/{{{l}_{p}}}\;$); or in the food bulk (${x}/{{{l}_{p}}}\;>{1+{{x}_{BL}}}/{{{l}_{p}}}\;$) without any significant concentration gradient (*e.g.*, assumption of a perfectly mixed medium). Choosing ${{{x}_{BL}}}/{{{l}_{p}}}\;$ as a free parameter enables to cover almost all contacting phases (gas, liquid food or simulant, solid and semi-solid food) with reasonable complexity.

Figure 7 plots simulation results using the concepts of statistical physics (*i.e.*, the molecules jump randomly vertically and horizontally without “knowing” where the interface is located) and by using the concept of continuum mechanics (*i.e.*, balance on populations and macroscopic fluxes on elementary volumes). The stochastic and deterministic point of views are equivalent and highlight that the observed macroscopic gradients are the consequence of the evolution of the distributions of solutes with time and not its cause. In the upward direction, the random displacements are compensated by the same and opposite microscopic flux in the downward direction. The net balance is zero, and no concentration gradient can develop. The substances translate at the same frequency in the horizontal direction (*i.e.*, isotropic diffusion), but since no solute comes to compensate the flux from left to right at the beginning of the contact, a net flux develops from left to right, resulting in the spreading of a concentration gradient from the source (the polymer: $P$) to the food (the food: $F$). Statistical physics and continuum mechanics counts molecules in a very similar fashion, via the concepts of probability density,  $\rho =\frac{N}{{{N}_{0}}V}$ and of volume concentration $C=\rho {{N}_{0}}$, respectively. ${{N}_{0}}$ is the total number of migrating substances in the whole system and $N$ the number of molecules in an elemental volume $V$. The concentration at the interface between F and P (denoted P-F) requires specific treatment and analysis. Since the principle zero of thermodynamics does not hold for mass transfer, both density and concentration are discontinuous at the P-F interface. If no reaction occurs at the interface, the mass balance is kept across the interface (*i.e.*, no substance is lost). Additionally, the principle of microscopic reversibility entails that the amount of substances crossing the interface per time unit from left to right (denoted P→F)  is exactly matched by the amount of substances crossing the interface in the opposite direction (denoted F→P). In other words, any substance located exactly at the interface  $x={{l}_{p}}$has the same probability to go in P and F, irrespectively its origin. This principle of microscopic mass balance reads:
$$
{{\rho }_{F}}{{f}_{P\to F}}={{\rho }_{P}}{{f}_{F\to P}}\label{eq2}
$$
where ${{f}_{A\to B}}$ is the frequency of translating from the compartment A to the compartment $B$. The concept of local chemical equilibrium developed among others by Henry Eyring provides a robust framework to express the frequencies of passage from one state to another one without justifying the details on how the change in conformations and local velocities affect the passage from A to B. The frequency of passage is written as:
$$
{{f}_{A\to B}}=\kappa \exp \left( -\frac{{{G}^{\ddagger }}-{{G}_{A}}}{RT} \right)\label{eq3}
$$
with ${{G}^{\ddagger }}$ the free energy associated to the transition state and ${{G}_{A}}$ the free energy of sorption of the
solute in A. The preexponential factor $\kappa
$ is independent of temperature; its expression depends on the statistical ensemble used to express probabilities. By combining Eqs  and , the molecules distribute across the interface with a ratio of probabilities equal to:
$$
{{K}_{F/P}}=\frac{{{\rho }_{F}}}{{{\rho }_{P}}}=\frac{C\left( x={{l}_{P}}+\epsilon ,t \right)}{C\left( x={{l}_{P}}-\epsilon ,t \right)}=\exp \left( -\frac{{{G}_{P}}-{{G}_{F}}}{RT} \right)\,\,\,\,  with \,\,\epsilon\rightarrow0\,\,and\,\,t>0\label{eq4}
$$
where $x$ and $t$ indicate the position and time; the food-packaging interface is located at $x=l_p$.
The concentration profiles are, hence, continuous across the P-F interface only when the free sorption energies are similar in both compartments. In the general case, the concentration profile is discontinuous across the interface. [Figure 7](#figure7) plots cases for an apparent partition coefficient ${{K}_{F/P}}$ of 0.5. In other words, the substance has a twice more chemical affinity for P than for F, as expected for plastic additives and monomers. For conservative migration modeling, values of ${{K}_{F/P}}$ higher than unity are preferred to maximize the gradient and consequently the amount transferred to the food. When the release of substance does not modify the properties of the polymer (polymer densification, a shift of the glass transition temperature ${{T}_{g}}$ ), the ratio of concentrations is likely to be constant at the interface.

<img src="./src/fig7.png" alt="fig7" style="zoom: 67%;" />

<a name="figure7"></a>**Figure 7**. One-dimensional description of solute diffusion (*e.g.*, additive, monomer) from the packaging wall (position: $0\leq x/l_p\leq1$, individual solutes identified as ×) to the contacting phase (individual solutes identified as ■) via the food boundary layer (individual solutes identified as ●): (a) random distribution of solutes and corresponding concentration profile at     and (b) after contact times up to  $Fo=2$ . The percentages in the top part represent the residual amount in each compartment.



#### 4.2.1.2 Migration modeling and similitude properties

Under the assumption of uniform and constant transport and thermodynamic properties in each compartment (polymer: P, boundary layer: BL, bulk contacting phase or food: F), the mass transfer problem is self-similar according to a small number of dimensionless parameters or ratios. According to the principle of similitude, a real problem can be compared to a theoretical case without dimensions if all independent dimensionless quantities are similar. The key dimensionless quantities are reviewed in [Table 4](#table4).

<a name="table4"></a>Table 4. Key dimensionless quantities of the migration from monomaterials. Contact is assumed to be initiated at $t=0$ .

| Dimensionless quantity                                       | Meaning                      | Justification                                                |
| ------------------------------------------------------------ | ---------------------------- | ------------------------------------------------------------ |
| $u\left( x,t   \right)=\frac{C\left( x,t \right)}{{{C}_{ref}}}$ | Dimensionless concentration  | ${{C}_{ref}}$ is a reference concentration, usually the initial   concentration in the polymer $C_{P}^{t=0}$ . |
| ${{K}_{F/P}}$   it is defined from Eq..  (4)                 | Partition coefficient        | At macroscopic equilibrium, it is also defined as ${{K}_{F/P}}=\frac{{{C}_{F}}\left(   t\to \infty  \right)}{{{C}_{P}}\left(   t\to \infty  \right)}$ with ${{C}_{F}}\left(   t \right)$ and ${{C}_{P}}\left( t   \right)$ the volume-averaged concentrations in F and P, respectively. |
| $q=\frac{x}{{{l}_{p}}}$                                      | Dimensionless position       | $l_p$ is characteristic food dimension, usually the   thickness or the ratio $\frac{V_p}{A}$ if the geometry is complex, with   ${{V}_{P}}$ the volume of the material and $A$ the surface area in contact with $F$. It is recommended to maximize $A$ by   also including the headspace in the   calculations of $A$. |
| $Fo=\frac{{{D}_{P}}t}{l_{p}^{2}}=\int\limits_{0}^{t}{\frac{{{D}_{p}}\left(   t \right)dt}{l_{p}^{2}}}$ | Dimensionless time           | ${{D}_{P}}$ is the diffusion coefficient in the polymer at the   temperature of contact. The integral form needs to be preferred if the   diffusion coefficient is variable with time (temperature change) |
| $Bi$ or $Sh$   $=\frac{{{D}_{F}}}{{{D}_{P}}}\frac{{{l}_{p}}}{{{x}_{BL}}}=\frac{h{{l}_{p}}}{{{D}_{P}}}$ | mass Biot or Sherwood number | where ${{x}_{BL}}$  is the   thickness of the boundary layer and ${{D}_{F}}$ the diffusion coefficient in   F, such that $h=\frac{{{D}_{F}}}{{{x}_{BL}}}$ is   the effective mass transfer coefficient across the boundary. |
| ${{L}_{P/F}}=\frac{{{l}_{p}}}{{{l}_{F}}}=\frac{{{V}_{P}}/A}{{{V}_{F}}/A}$ | Dilution ratio               | This number is defined as the   ratio of characteristic lengths and controls how the substances are diluted in the food, usually much larger   than the packaging. |

 

#### 4.2.1.3     Explicit vs implicit food representation

It would be logical that migration models describe, how the migrants distribute in the food explicitly. Solid foods, such as a chicken or a pizza, are not expected to have all parts contaminated similarly. For risk assessment, we consider that all parts are intended to be ingested, including the most contaminated sauce in contact with the packaging. As a result, only a global estimate of the food contamination is required, as measured with a liquid simulating the entire food. Replacing a solid by a liquid or vice-versa has a consequence on the rates of mass transfer. This part discusses the differences between explicit and implicit representations of the food and of the risk to underestimate the real migration.

[Figure 7](#figure7) represents explicitly mass transfer in the food, that is the concentration profiles in the food are also calculated. The depicted cases correspond to a characteristic food length ${{l}_{F}}=\frac{{{V}_{F}}}{A}$ of $12{{l}_{p}}$ (**i.e.**, thin food to make the boundary layer visible). The total domain has a length of $13{{l}_{p}}$. When the food is represented explicitly, the amount transferred to the food is defined as:
$$
{{C}_{F}}\left( t \right)=\frac{1}{{{l}_{F}}}\int\limits_{x={{l}_{P}}}^{x={{l}_{p}}+{{l}_{F}}}{C\left( x,t \right)dx}\label{eq5}
$$

Eq.(6)  accumulates substances in the boundary layer (round symbols) and in the food bulk (square symbols). Implicit food representation will describe mass transfer only in the packaging and apply a proper boundary condition between the food and the packaging at the position $x={{l}_{p}}$. With the help of Eq. , the flux at the interface, denoted $j$,  can be expressed only with concentrations taken inside the packaging or in the food far from the interface,  $C\left( x\to \infty ,t \right)$.This assumption open the way of an implicit representation of the food via a boundary condition relating the diffusion at the packaging-food interface with the flux entering into the food:

$$
\begin{align}
  & j\left( t \right)=-{{\left. {{D}_{P}}\frac{\partial C}{\partial x} \right|}_{x=lp-\epsilon ,t}}=h\left( C\left( x+\epsilon ,t \right)-C\left( x\to \infty ,t \right) \right)  =h\left( {{K}_{F/P}}C\left( x-\epsilon ,t \right)-C\left( x\to \infty ,t \right) \right)  
\,\,\,\,  with \,\,\epsilon\rightarrow0\,\,and\,\,t>0\end{align}\label{eq6}
$$

 Eq. (6) offers a good approximation of the explicit representation when the contact time is sufficient to reach a fully developed concentration profile (linear, so-called Prandtl approximation) inside the boundary layer. The critical Fourier number is given by $F{{o}_{critical}}={{{\left( {{{x}_{BL}}}/{{{l}_{p}}}\; \right)}^{2}}}/{\left( 6{{{D}_{F}}}/{{{D}_{P}}}\; \right)=}\;{\left( {{{x}_{BL}}}/{{{l}_{p}}}\; \right)}/{\left( 6Bi \right)}\;$. The profile plotted in f[Figure 7a](#figure7) for $Bi=1$  and $Fo=0.1$ deviates from the assumption above. The critical Fourier number is, indeed, of 1.67 and the value of $C\left( x\to \infty ,t \right)$ is close to zero, whereas the food is already contaminated via its boundary layer.
Implicit food representation approximates the concentration in the food, ${{C}_{F}}\left( t \right)$ , by its concentration far from the interface. By noticing that the flux $j\left( t \right)$ is taken after the boundary layer, one gets:
$$
{{C}_{F}}\left( t \right)=C\left( x\to \infty ,t \right)\approx {{C}_{F}}\left( t=0 \right)+\int\limits_{0}^{t}{j\left( t \right)dt}\,\,\,\,for\,\,t>F{{o}_{critical}}\frac{l_{p}^{2}}{{{D}_{P}}}\label{eq7}
$$
For $Bi>50$ and if the threshold of concern is not too low, the amount present in the boundary layer can be neglected (less than 1% in Figure 7a when $Bi=50$ ) and an implicit representation can be used. Its main advantage is the dramatic reduction of the problem complexity and of the computational time. In very thin or low barrier films and in solid foods, the implicit food representation may underestimate the contamination of the food severely. The amounts in the boundary layer reported in [Figure 7a](#figure7) reach 9%, 18% and 22% for $Bi$= 10, 5, 1 respectively. When an implicit representation is used, calculating the concentration in the food from the mass balance in the packaging between $t=0$  and $t$ does not solve the issue as the closure equalities are enforced at any time:
$$
\frac{d}{dt}\left( \int\limits_{0}^{{{l}_{p}}}{C\left( x,t \right)dx} \right)=j\left( t \right)=\frac{d{{C}_{F}}}{dt}=\frac{dC\left( x\to \infty ,t \right)}{dt}\,\,\,\,for\,\,t>0\label{eq8}
$$
The amount present in the boundary layer is always neglected in migration representing the food implicitly. Only by choosing artificially $Bi\to \infty $as a worst-case scenario makes this amount negligible at the price of migration much faster than the one expected in the real conditions. 

#### 4.2.1.4     Other assumptions

We describe in this section the constitutive equations to describe mass transfer from monolayer and multilayer materials, when the food is represented implicitly via the boundary condition  and the food mass balance (8). The total packaging thickness is denoted ${{l}_{p}}=\sum\limits_{i=1}^{m}{{{l}_{i}}}$ for a packaging (*e.g.* laminate) consisting of $m$ layers. At the position $x=0$ (usually the surface exposed to the ambiance), an impervious boundary layer is assumed so that the amount transferred to the food is maximized. All substances are assumed to be well dispersed where they have been incorporated (*e.g.*, no blooming effect) and below their concentration at saturation (*i.e*., no supersaturation). 

### 4.2.2     Governing equations for monolayer materials

#### 4.2.2.1     Overview

The full set of equations for monolayer materials including the initial condition (IC), the transport equation (TE), the boundary conditions (BC) and the mass balance on the food compartment (MB) reads:
$$
\left\{ \begin{align}
  & IC:{{C}_{\left( x,t=0 \right)}}=C_{P}^{t=0}\,\text{for}\ 0\le x\le {{l}_{p}} \label{eq9}\\ 
 & TE:\frac{\partial {{C}_{\left( x,t \right)}}}{\partial t}=\frac{1}{{{x}^{n}}}\frac{\partial }{\partial x}\left( {{D}_{P}}\left( {{C}_{\left( x,t \right)}} \right){{x}^{n}}\frac{\partial {{C}_{\left( x,t \right)}}}{\partial x} \right)\text{for}\ 0\le x\le {{l}_{p}},t>0\label{eq10} \\ 
 & BC:j\left( t \right)=-{{\left. \left( {{D}_{P}}\frac{\partial C}{\partial x} \right) \right|}_{x={{l}_{p}}-\epsilon ,t}}=h\left[ {{K}_{F/P}}{{C}_{\left( x={{l}_{p}}-\epsilon ,t \right)}}-{{C}_{F}}\left( t \right) \right]; \label{eq11}\\ 
 & {{\left. \,\,\frac{\partial C}{\partial x} \right|}_{x=0,t}}=0\ \text{for}\ t>0\ \text{and}\ \epsilon \to 0\, \label{eq12}\\ 
 & MB:{{C}_{F}}\left( t \right)={{C}_{F}}\left( t=0 \right)+\frac{A}{{{V}_{F}}}\int\limits_{0}^{t}{j\left( \tau  \right)d\tau }\ \text{for}\ \text{t}\ge 0 \label{eq13} \\
 \end{align} \right.
$$
where $n$ is an exponent controlling the system of coordinates ($n=0$ : Cartesian, $n=1$ : cylindrical; $n=2$ : spherical); ${{D}_{P}}\left( {{C}_{\left( x,t \right)}} \right)$ is the diffusion coefficient possibly variable with concentration (*e.g.*, case of plasticizers used at high concentrations). 

When the diffusion coefficient in the packaging is considered constant along with ${{C}_{F}}\left( t=0 \right)=0$, [Table 4](#table4) and system 9-13 yield the following dimensionless formulation for Cartesian coordinates:
$$
\left\{ \begin{align}
  & IC:{{u}_{\left( q,Fo=0 \right)}}=1\,\text{for}\ 0\le q\le 1 \label{eq14}\\ 
 & TE:\frac{\partial {{u}_{\left( q,Fo \right)}}}{\partial Fo}=\frac{{{\partial }^{2}}{{u}_{\left( q,Fo \right)}}}{\partial {{q}^{2}}}\text{for}\ 0\le q\le 1,Fo>0 \label{eq15}\\ 
 & BC:{{j}^{*}}\left( Fo \right)=-{{\left. \frac{\partial {{u}_{\left( q,Fo \right)}}}{\partial q} \right|}_{q=1-\epsilon }}=Bi\left[ {{K}_{F/P}}{{u}_{\left( q=1-\epsilon ,Fo \right)}}-\frac{{{C}_{F}}\left( t \right)}{C_{P}^{t=0}} \right]; \label{eq16}\\ 
 & {{\left. \,\,\frac{\partial {{u}_{\left( q,Fo \right)}}}{\partial q} \right|}_{q=0}}=0\ \text{for}\ t>0\ \text{and}\ \epsilon \to 0\, \label{eq17}\\ 
 & MB:\frac{{{C}_{F}}\left( t \right)}{C_{P}^{t=0}}={{L}_{P/F}}\int\limits_{0}^{Fo}{{{j}^{*}}\left( \tau  \right)d\tau }\ \text{for}\ \text{t}\ge 0 \label{eq18}\\ 
\end{align} \right.
$$
When $Bi\to \infty $,  the third type boundary condition at $q=1$ (Robin boundary condition) can be replaced by a simple coupling with the mass balance. By differentiating BM respectively to $Fo$, one gets:
$$
{{j}^{*}}\left( Fo \right)=-{{\left. \frac{\partial {{u}_{\left( q,Fo \right)}}}{\partial q} \right|}_{q=1-\epsilon }}=\frac{1}{{{L}_{P/F}}C_{P}^{t=0}}\frac{d{{C}_{F}}\left( t \right)}{dFo}=\frac{{{K}_{F/P}}}{{{L}_{P/F}}}\frac{\partial {{u}_{\left( q=1,Fo \right)}}}{\partial Fo}
$$
The worst-case scenario, $Bi\to \infty $and ${{K}_{F/P}}\to \infty $corresponds to the Dirichlet boundary condition ${{u}_{\left( q=1,Fo \right)}}=0$.

#### 4.2.2.2     Concentration in the contact phase at thermodynamical equilibrium

According to Eqs. , the maximum concentration in food is obtained at thermodynamical equilibrium (${{j}^{*}}\to 0$ ):


$$
\frac{C_{F}^{eq}}{C_{P}^{t=0}}=\frac{{{C}_{F}}\left( Fo\to \infty  \right)}{C_{P}^{t=0}}=\frac{1}{{1}/{{{L}_{P/F}}}\;+{1}/{{{K}_{F/P}}}\;}
$$
 In practice it is convenient to express the kinetics of desorption in food as a function of the fraction of the equilibrium value: 


$$
\begin{align}
  & {{C}_{F}}\left( Fo,Bi,{{K}_{F/P}},{{L}_{P/F}} \right)=\frac{C_{P}^{t=0}}{{1}/{{{L}_{P/F}}}\;+{1}/{{{K}_{F/P}}}\;}{{{\bar{v}}}^{*}}\left( Fo,Bi,{{K}_{F/P}},{{L}_{P/F}} \right)=C_{F}^{eq}{{{\bar{v}}}^{*}}\left( Fo,Bi,{{K}_{F/P}},{{L}_{P/F}} \right)  
\end{align}
$$
​                                                   

#### 4.2.2.3     Dimensionless migration kinetics and their analytical approximations

The dimensionless migration kinetics ${{\bar{v}}^{*}}\left( Fo,Bi,{{K}_{F/P}},{{L}_{P/F}} \right)$ are plotted in Figure 8. The analytical solution associated with the Dirichlet condition $C\left( x={{l}_{p}},t \right)=0$ is given by Eq. 4.18 in Crank’s book [102](#_ENREF_102) and reads :

$$
{{\bar{v}}^{*}}\left( Fo \right)=\frac{{{C}_{F}}\left( Fo \right)}{C_{F}^{eq}}=1-\frac{8}{{{\pi }^{2}}}\underset{n\to \infty }{\mathop{\lim }}\,{{S}_{n}}\left( Fo \right) \,\,\,\,with \,\,\,\,{{S}_{n}}\left( Fo \right)=\sum\limits_{i=0}^{n}{\frac{\exp \left( -\tfrac{{{\pi }^{2}}}{4}{{\left( 2i+1 \right)}^{2}}Fo \right)}{{{\left( 2i+1 \right)}^{2}}}}
$$

For small $Fo$ values, the approximation  requires $n$ to be very large (103 or 104) and a more efficient approximation can be obtained by combining an approximation of Eq 4.20 in Crank’s book [102](#_ENREF_102) with Eq.  for $i=0$ as:

$$
{{\bar{v}}^{*}}\left( Fo \right)=\min \left[ \frac{2\sqrt{Fo}}{\sqrt{pi}},1-\frac{2}{{{\pi }^{2}}}\exp \left( -\frac{{{\pi }^{2}}}{4}Fo \right) \right]\le \min \left[ \frac{2\sqrt{Fo}}{\sqrt{pi}},1 \right]
$$

Approximations  and  are plotted along with the results of numerical simulations in  [Figure 8](#figure8). The common assumption of the linearity of ${{\bar{v}}^{*}}$ with $\sqrt{Fo}$ is well verified while $Fo\le 0.8$ and when the equivalent length of the contacting phase ${{\ell }_{F}}={{{K}_{F/P}}}/{{{L}_{P/F}}}\;$ and $Bi$ are large. At intermediate$Bi$ values, the Dirichlet condition offers a conservative approximation (*i.e.*, ${{\bar{v}}^{*}}$ is overestimated) and Eq.  can be used safely for compliance testing. However, when ${{\ell }_{F}}\ll 100$ (*e.g.*, small food volume, low chemical affinity for the food), Eq.  must be avoided due to a significant risk of underestimation of low $Fo$ values. The equilibrium is, indeed, reached much faster due to a much smaller amount to transfer and because the concentration at the F-P interface never vanishes. An accurate estimation of migration requires ad-hoc numerical solutions or special analytical solutions. In mathematical terms, finite volume effects cause the propagation of shockwaves between the polymer and the contacting phase: the addition of substances to F modifies instantaneously the capacity of P to transfer additional substances, and so on for the next substances creating positive and destructive interference across the mass transfer boundary layer when it exists. One practical consequence is that analytical solutions are without closed-forms and therefore more complex than Eq. XXX.

- For $Bi\to \infty $ and ${{L}_{P/F}}<1$, the analytical solution verifying finite volume effects and boundary condition  is given by Eq. 4. 37 in Crank’s book [102](#_ENREF_102) and is very accurate at a reasonable cost for $Fo>{{10}^{-4}}$.

- For arbitrary $Bi$ and ${{\ell }_{F}}$ values, new solutions verifying the general boundary condition  (see Refs. [127](#_ENREF_127), [177-179](#_ENREF_177)). Short-time and long-time solutions were optimized for efficiency and to integrate more complex physics such as non-linear sorption isotherms or for boundary conditions variable in time. 


 The general solutions are not detailed here as their expressions exceed the scope of the article. When $Bi\to \infty $, the Eq. 4. 37 in Crank’s book [102](#_ENREF_102) reads:

$$
{{\bar{v}}^{*}}\left( Fo \right)=1-\sum\limits_{n=1}^{\infty }{2\frac{{{\ell }_{F}}\left( 1+{{\ell }_{F}} \right)}{1+{{\ell }_{F}}+\ell _{F}^{2}q_{n}^{2}}\exp \left( -q_{n}^{2}Fo \right)} where {{q}_{n}}$are zeros of $\tan {{q}_{n}}=-{{\ell }_{F}}{{q}_{n}}
$$

The zeros of the transcendental equation, ${{q}_{n}}$, increase with $n$ and also when ${{\ell }_{F}}$ decreases. This behavior demonstrates that ${{\bar{v}}^{*}}\left( Fo \right)$ converges exponentially and more rapidly to equilibrium when${{\ell }_{F}}$values are low. The linearization with $\sqrt{Fo}$ ceases to be correct earlier and is associated to slopes varying with ${{\ell }_{F}}$. As the values of ${{q}_{n}}$are usually not tabulated beyond $n=6$in reference text books (see Table 4.1, p. 379 in [102](#_ENREF_102), it is recommended to restrict the use of Eq. XXXX  to $Fo>{{10}^{-4}}$ . 

 

![fig8](./src/fig8.png)

<a name="figure8"></a>**Figure 8**. Dimensionless desorption kinetics ${{\bar{v}}^{*}}\left( Fo \right)=\frac{{{C}_{F}}\left( Fo \right)}{C_{F}^{eq}}$  for various values of ${{L}_{P/F}}$, ${{K}_{F/P}}$ and $Bi$ with $C_{F}^{eq}=\frac{C_{P}^{0}}{1/{{L}_{P/F}}+1/{{K}_{F/P}}}$ .  Approximations [1] and [2] are given by Eqs. XXX  and XXX , respectively.

#### 4.2.3 Governing equations for multilayers



------
Highlights for multilayers
- Multilayer can be seen as a generalization of monolayer materials but with additional features such as functional barriers (delayed migration) and reservoir layers (accumulation inside one or several internal virgin layers).
- When the partition coefficients between layers and their diffusion coefficients are constant with time, the total migration is the sum of the migrations associated with the contribution of individual layers.
- Uncertainty in partitioning and the initial distribution of migrants can be overcome by moving “artificially” the content on one layer to a layer closer to the contacting phase.
- The calculation procedures presented in this section have been introduced in the guidance to EU Regulation 10/2011/EC.

------


4.2.3.1	Thermodynamical assumptions
The case of materials consisting of $m$ materials or layers ($m>1$ ) can be seen as a generalization of the monolayer case  ($m=1$ ) at the expense of a few additional assumptions and conventions. Because monolayer systems were dominating in the 20th century, the reasoning supporting US and EU regulations was, hence, devised based on an assumption of a migration without delay and obeying to a scaling proportional to the square root of time. The conventional condition of ten days at 40°C was thus thought to be equivalent to a test of one hundred days with a factor comprised between unity (equilibrium) and $\sqrt{\tfrac{100\ \text{days}}{10\ \text{days}}}\approx 3.16$ ($Fo\ll 1)$. Moving the substances away from the P-F interface delays substantially the mass transfer to the contacting phase. This behavior is central to the concept of functional barrier 180; it was initially explored to promote the incorporation of recycled material – possibly contaminated due to post-consumer misuse or mixing with non-food grade materials – in co-extrudates. The recycled polymer is sandwiched within two layers of virgin polymers. Similar problems can be resolved only by adapting the initial condition $C\left( x,t=0 \right)$  in the equation system (9) to the need instead of using a uniform distribution.
When the materials are different in nature, a contact condition similar to Eq. (4) needs to be considered. The proposed description relies on an assumption of a linear and reversible sorption of substances in each layer so that a linear sorption isotherm is assumed in any layer, including in the food. By denoting $p\left( x,t \right)$ the partial pressure of the migrating substance, the isotherm associated to the Henri constant ${{\left\{ {{k}_{j}} \right\}}_{j=0..m}}$ reads for each layer:
$$
p\left( x,t \right)={{k}_{j}}C\left( x,t \right)\,\,\,\,for\,\,\,\,0\le {{L}_{j-1}}\le x<{{L}_{j}}\,\,\,\,and\,\,\,\,t>0\,\,\,\,;\,\,\,\,with\,\,\,\,{{L}_{j}}=\sum\limits_{i=1}^{j}{{{l}_{i}}}\,\,\,\,and\,\,\,\,0\le j\le m
$$
(17)
With a sake of generalization, the layers including the food are indexed from $j=0$ (food: F) to $j=m$ (layer the most distant to F), with $j=1$ being the contact layer, as shown in [Figure 9](#figure9). The residual concentration in each layer is ${{C}_{j}}\left( x,t \right)=\frac{1}{{{l}_{j}}}\int\limits_{{{L}_{j-1}}}^{{{L}_{j}}}{C\left( x,t \right)dt}$ for $1\le j\le m$and ${{C}_{0}}\left( t \right)={{C}_{F}}\left( t \right)$.

<img src="./src\fig9.png" alt="fig9" style="zoom:67%;" />

<a name="figure9"></a>**Figure 9**. Indexing rule of a material including m layers (total thickness $l_p$) in contact with a food product indexed 0. The left and right external boundaries are considered impervious (no mass loss). The concepts of “functional barrier” are “reservoir” assume that the layer j is the source (with non-zero initial concentration). It is used in §4.2.3.5 and §4.2.3.6.

#### 4.2.3.2	Concentration in the contact phase at thermodynamical equilibrium
The Henri isotherm defined in Eq. (17) offers a robust but simplified thermodynamic representation of the variation of the chemical potential with the local composition in the system. The validity of the model and its generalization are discussed in §6.2 . Partial pressure \[p\left( x,t \right)\] is a continuous potential, and thermodynamical equilibrium is achieved when its value is uniform across the structure. By neglecting mass losses, Eq. (17) and the mass balance between $t=0$ and $t\to \infty $ $\sum\limits_{j=0}^{m}{C_{j}^{t\to \infty }{{l}_{j}}}=\sum\limits_{j=0}^{m}{C_{j}^{t=0}{{l}_{j}}}$ enables to generalizes Eq. (12):
$$
C_{F}^{eq}=C_{0}^{t\to \infty }=\frac{\sum\limits_{j=0}^{m}{\frac{{{l}_{j}}}{{{l}_{0}}}C_{j}^{t=0}}}{1+\sum\limits_{j=1}^{m}{\frac{{{k}_{0}}}{{{k}_{J}}}\frac{{{l}_{j}}}{{{l}_{0}}}}}
$$
(18)
where ${{{k}_{0}}}/{{{k}_{j}}}\;={C_{j}^{t\to \infty }}/{C_{0}^{t\to \infty }}$ is the partition coefficient between the layer $j$ and the food.

#### 4.2.3.3	Transport equations
Transport equations are unchanged and are piecewise-defined:
$$
\frac{\partial C\left( x,t \right)}{\partial t}={{D}_{j}}\frac{{{\partial }^{2}}C\left( x,t \right)}{\partial {{x}^{2}}}\,\,\,\,for\,\,\,\,0\le {{L}_{j-1}}\le x<{{L}_{j}}\,\,\,\,and\,\,\,\,1\le j\le m
$$
(19)
and connected at internal interfaces via the double conditions of mass conservation and local thermodynamical equilibrium:
$$
\begin{align}
  & {{D}_{j}}{{\left. \frac{\partial C\left( x,t \right)}{\partial x} \right|}_{x={{L}_{j}}-\epsilon }}={{D}_{j+1}}{{\left. \frac{\partial C\left( x,t \right)}{\partial x} \right|}_{x={{L}_{j}}+\epsilon }}\ \text{for}\ 1\le j<m \\ 
 & {{k}_{j}}C\left( x={{L}_{j}}-\epsilon ,t \right)={{k}_{j+1}}C\left( x={{L}_{j}}-\epsilon ,t \right)\ \ \ \text{for}\ 0\le j<m 
\,\,\,\,and\,\,\,\,\epsilon \to 0\,\,\,\, , t>0\end{align}
$$
(20)

#### 4.2.3.4	Limiting mass transfer resistance
A dimensionless formulation of Eqs (6), (19) and (20) is achievable at the expense of choosing a reference layer, $1\le {{j}_{ref}}\le m$ , setting the reference time scale $\frac{l_{{{j}_{ref}}}^{2}}{{{D}_{{{j}_{ref}}}}}$. In numerical algorithms, where stability and convergence are very stringent, it is convenient to choose as a reference layer, the layer associated to the highest mass transfer resistance (lowest permeability) is the best choice:
$$
{{j}_{ref}}\ \text{such}\ \text{that}\ \frac{{{l}_{{{j}_{ref}}}}{{k}_{{{j}_{ref}}}}}{{{D}_{{{j}_{ref}}}}}=\max \left[ {{\left\{ \frac{{{l}_{j}}{{k}_{j}}}{{{D}_{j}}} \right\}}_{j=1..m}} \right]
$$
(21)
If several conditions need to be compared, a natural choice is to choose the contact layer as the reference layer (${{j}_{ref}}=1$).

#### 4.2.3.5     Typologies of migration behaviors

The main behaviors, which can be observed with multilayers, are illustrated in simple configurations corresponding to a bilayer structure (each layer has a thickness ${{{l}_{p}}}/{2}\;$ ) in contact with a liquid phase with a characteristic thickness $50{{l}_{p}}$ and associated to ${{k}_{0}}=1$. The five considered cases are summarized in [Table 5](#table5) and were associated with a similar initial amount and comparable final concentration in the contacting phase *ca*. $2.0\cdot {{10}^{-2}}\cdot {\left( {{C}_{1}}+{{C}_{2}} \right)}/{2}\;$. The other parameters were $Bi\to \infty $, ${{l}_{1}}={{l}_{2}}=\frac{{{l}_{p}}}{2}$ , ${{D}_{1}}={{D}_{2}}$  and the time scale $\frac{l_{p}^{2}}{4{{D}_{1}}}$. The calculated profiles and kinetics corresponding to the five scenarios are plotted in [Figure 10](#figure10).

<a name="table5"></a>**Table 5**. Illustration of the main behaviors associated with multilayer structures. The concepts of functional barrier and reservoir are illustrated in [Figure 9](#figure9).

| $C_{j}^{0}$ |       | ${{k}_{j}}$ |       | interpretation                                               | code        |
| ----------- | ----- | ----------- | ----- | ------------------------------------------------------------ | ----------- |
| $j=1$       | $j=2$ | $j=1$       | $j=2$ |                                                              |             |
| 1           | 1     | 1           | 1     | uniform   distribution (equivalent to a monolayer)           | [1,1]×[1,1] |
| 0           | 2     | 1           | 1     | functional   barrier (barrier to diffusion only)             | [0,2]×[1,1] |
| 2           | 0     | 1           | 1     | reservoir   layer (same capacity)                            | [2,0]×[1,1] |
| 0           | 2     | 2           | 1     | functional   barrier (barrier to diffusion and of solubility) | [0,2]×[2,1] |
| 2           | 0     | 1           | 2     | reservoir   layer (twice less capacity)                      | [2,0]×[1,2] |

Monolayers and functional barriers lead to uniformly decreasing concentration profiles. The corresponding desorption kinetics in F are respectively proportional to the square-root of time and proportional to time after some lag time equal to $\frac{l_{1}^{2}}{6{{D}_{1}}}$. For the same initial amount in the structure and after the lag time, the functional barrier ceases to operate and leads to a migration proportional to $\frac{{{D}_{1}}}{{{k}_{1}}{{l}_{1}}}$. Only a functional barrier combining a diffusion barrier ($\frac{{{D}_{1}}}{{{l}_{1}}}\ll \frac{{{D}_{2}}}{{{l}_{2}}}$ ) and solubility one ($\frac{{{k}_{1}}}{{{k}_{0}}}\gg \frac{{{k}_{2}}}{{{k}_{0}}}$ ) can slow down the desorption durably. Reservoirs behave very differently; they are associated with non-monotonous concentration profiles, accelerated desorption kinetics while converging to a very similar concentration at equilibrium. For the same initial content and when the barrier on the right is higher than the barrier on the left $\frac{{{k}_{2}}{{l}_{2}}}{{{D}_{2}}}\gg \frac{{{k}_{0}}}{Bi}$ , the “*reservoir*” configuration overestimates the migration kinetics associated with all other configurations.

 

<img src="./src/fig10.png" alt="fig10" style="zoom: 67%;" />

<a name="figure10"></a>**Figure 10**. Concentration profiles (top) and migration kinetics (bottom) for the bilayer structure and scenarios detailed in [Table 5](#table5).



#### 4.2.3.6     Superposition principles and conservative scenarios for multilayer and multicomponent systems

 

**Mathematical principles**

Multilayer structures offer a broad range of behaviors. In the simplest cases, as shown in [Figure 10](#figure10), desorption kinetics are monotonous with time. It may be not the case if the functional barrier and reservoir effects are combined. The calculations for complex multilayer are complicated by the difficulty to associate the uncertainty in diffusion ${{\left\{ {{D}_{j}} \right\}}_{j=1..m}}$  and sorption properties ${{\left\{ {{k}_{j}} \right\}}_{j=1..m}}$on the concentration in food. For monolayer materials, a conservative scenario is achieved by overestimating simultaneously ${{D}_{1}}$ ,$\frac{{{k}_{0}}}{{{k}_{1}}}$ and $Bi$. For multilayer materials, no similar rule holds. Intuitively based on the illustrated configurations in [Figure 10](#figure10), it can be stated that for the capacity of the layer $j$ to transfer its content, ${{M}_{j}}\left( t \right)$ (amount transferred to F at the time $t$), is maximized (denoted $\left\lceil {{M}_{j}} \right\rceil \left( t \right)$ ) if the following properties are fulfilled:

$$
\begin{align}

  & {{D}_{i}}\,\text{is}\ \text{replaced}\ \text{by}\ \left\lceil {{D}_{i}} \right\rceil \text{}{{\text{D}}_{i}}\ \ \text{for}\ 1\le i\le j\ \text{and}\ {{k}_{i}}\ \text{is replaced by }\left\lfloor {{k}_{i}} \right\rfloor \ <{{k}_{i}}\ \text{for}\ 0\le i\le j \\ 

 & {{D}_{i}}\,\text{is}\ \text{replaced}\ \text{by}\ \left\lfloor {{D}_{i}} \right\rfloor \text{}{{\text{D}}_{i}}\ \ \text{and}\ {{k}_{i}}\ \text{is replaced by }\left\lceil {{k}_{i}} \right\rceil \ >{{k}_{i}}\ \text{for}\ j<i\le m\  \\ 

\end{align}
$$

Conditions  have a mathematical justification in the linear properties of the equations -. The solution of any linear decomposition of the initial concentration profile is equal to the sum of the individual solutions:

$$
\begin{align}

  & C\left( t=0,x \right)=\sum\limits_{j=1}^{p}{{{C}_{j}}\left( t=0,x \right)}\ \text{where}\ p\ge 1\ \text{is}\ \text{the}\ \text{number}\ \text{of}\ \text{profiles, }0\le x\le {{l}_{p}} \\ 

 & {{C}_{F}}\left( t,C\left( t=0,x \right) \right)=\frac{M\left( t,C\left( t=0,x \right) \right)}{{{V}_{F}}}=\frac{\sum\limits_{j=1}^{p}{{{M}_{j}}\left( t,{{C}_{j}}\left( t=0,x \right) \right)}}{{{V}_{F}}} \\ 

 & \le \frac{\sum\limits_{j=1}^{p}{\left\lceil {{M}_{j}} \right\rceil \left( t,{{C}_{j}}\left( t=0,x \right) \right)}}{{{V}_{F}}}=\left\lceil {{C}_{F}} \right\rceil \left( t,C\left( t=0,x \right) \right)  

\end{align}
$$


**Example for a trilayer material ABC**

Eq. XXX is particularly significant as it is valid for any partitioning of the source terms in the material, respecting or not the positions of the layers. An application of the additivity of${{M}_{j}}$values (concentration profiles and kinetics) is shown in [Figure 11](#figure11) for a trilayer structure ABC detailed in [Table 6](#table6). 



<a name="table6"></a>**Table 6**. Parameters used to construct realistic and conservative migration scenarios depicted in [Figure 11](#figure11) and [Figure 12](#figure12). Quantities are expressed respectively to the likely values† for the first layer (the three layers ABC are indexed 1,2,3). They are scalar when the contribution of each layer as a source is considered in combination with others (the three sources are considered at once). The contributions of individual sources are indicated by 3×1 vectors mentioning the properties of all layers considered as a source or not.

| property                                                     | **case-study (likely)**   contribution of the *j*th layer | **conservative scenario (high tier)**   contribution of the *j*th layer | **worst-case scenario (low tier)**   contribution of the *j*th layer |               |            |           |       |       |       |
| ------------------------------------------------------------ | --------------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------- | ---------- | --------- | ----- | ----- | ----- |
|                                                              | $j=1$                                                     | $j=2$                                                        | $j=3$                                                        | $j=1$         | $j=2$      | $j=3$     | $j=1$ | $j=2$ | $j=3$ |
| $\frac{{{l}_{j}}}{l_{1}^{likely}+l_{2}^{likely}+l_{3}^{likely}}$ | 0.2                                                       | 0.5                                                          | 0.3                                                          | 0.2           | 0.5        | 0.3       | 0.2   | 0.5   | none  |
| ${C_{j}^{0}}/{C_{1}^{0,likely}}\;$                         | 1                                                         | 1                                                            | 1                                                            | [1,0,0]       | [0,1,0]    | [0,0,1]   | 3.5   | 0.6   | none  |
| ${{{D}_{j}}}/{D_{1}^{likely}}\;$                             | 1                                                         | 1                                                            | 1                                                            | [1,10-3,10-3] | [1,1,10-3] | [1,10,10] | 1     | 1     | none  |
| ${{{k}_{j}}}/{k_{0}^{likely}}\;$                           | 0.3                                                       | 0.8                                                          | 2.0                                                          | [10,10,10]    | [1,10,10]  | [1,1,10]  | 10    | 10    | none  |

†Likely value = true value or close to the true value in the considered scenario.



The simulation of each layer individually underlines the different mechanisms controlling the contribution of each layer: reservoir effect for source A (scaling of desorption kinetics with the square root of time) and functional barriers for B and C (desorption kinetics linear with time after significant lag-times). The depicted profiles are assumed to the “likely” or “true” ones. They are considered inaccessible to simulation and should be approximated at some tier in a way that the concentration in F is always overestimated (see [Figure 13](#figure13)). 



<img src="./src/fig11.png" alt="fig11" style="zoom:80%;" />

<a name="figure11"></a>**Figure 11**. Illustration of the additivity of the sources (see Eq. ) for a trilayer structure ABC associated with the case study detailed in Table 6: concentration profiles (top), kinetics (bottom). The case “sources ABC” is obtained by simulating the whole structure. The result A+B+C corresponds to the mathematical addition of the contributions of the three sources.

 

Eq. XXX provides a numerical procedure to devise conservative scenarios for multilayer structures. A similar procedure has been detailed in the section 4.2 of the European guidance document [143](#_ENREF_143). We repeat here the procedure for the sole overestimation of the chemical affinity effects and by keeping the diffusion coefficients to their “likely” values. The core idea is to prevent or hinder the diffusion of the substances in the *j*^th^ layer to the right (by assuming that the food is on the left) and to facilitate their desorption to the left, to bring the contaminants closer to the food. The “*conservative scenario*” of [Table 6](#table6) applies a factor ten to the Henry-like coefficient(s) ${{k}_{j}}$of the source and behind. The likely ${{k}_{j}}$ and  $D_j$ values are kept for the layers between the food and the *j*^th^ layer. To prevent back diffusion in the reservoir layers the diffusion coefficients were divided arbitrarily by a factor 103. The corresponding kinetics are shown in [Figure 12](#figure12), with their parameters listed in the section “conservative scenario” of [Table 6](#table6). The diffusion coefficients towards the food are overestimated by a factor 10. At intermediate Fourier numbers below 0.4, the kinetics is not significantly overestimated, but a significant conservatism is achieved at larger Fourier numbers when the overall migration is controlled by the second and the third layers as sources. This case study demonstrates how migration scenarios can be finely tuned to decrease the uncertainty related to the behavior of internal layers. The contact layer as a source can always be considered as a single monolayer.

 

<img src="./src/fig12.png" alt="fig12" style="zoom:80%;" />

<a name="figure12"></a>**Figure 12**. Illustration of the conservative scenario of Table 6 based on the overestimation of the contribution of each source. The reference corresponds to the initial case-study configuration also depicted in Figure 11.

 

In practice, any uncertainty on the internal partitioning between layer can be converted into a conservative scenario by forcing mass transfer to the food and by relocating “artificially” the content from the layer $j$ to layer $j-1$. The iterative procedure is illustrated in [Figure 13](#figure13) and can be applied iteratively to decrease a $m$-layer problem into a $m-1$ , $m-2$ , etc. layer problem, until to reach tier 1 (total migration). The first iteration applied to the case-study is denoted “worst-case” scenario in [Table 6](#table6). The procedure may overestimate dramatically the real migration, but it keeps the applicability of conservative calculations to demonstrate safety at low cost. For risk assessment, the procedure may be applied with caution as it may lead to unrealistic consumer exposure.

<img src="./src/fig13.png" alt="fig13" style="zoom:80%;" />

<a name="figure13"></a>**Figure 13**. Principles of the simplification of a $m$-layer problem (here $m=3$) into a problem with a lower number of layers and, therefore, easier parameterization. The represented distributions in the packaging corresponds to initial conditions at various tiers.



### 4.2.4     Strategies and equations to simulate multiple steps and conditions 

This part discusses the invariance of migration estimates, namely ${{C}_{F}}\left( t \right)$ , with the order of contact and thermodynamical conditions met by the components of the packaging before and during food contact.

------

**Highlights on the strategies to simulate multiple steps and variable conditions.**

- The initial dispersion of the substances between materials and the subsequent migration into the contact phase requires at least a two-step modeling and simulation.
- Temperature and relative humidity are variable during storage, transportation, etc.  and subjected to uncontrolled diurnal and seasonal variations. Their variations can be integrated in chained simulations (one condition = one step), where the outputs calculated at the previous step are used as inputs at the next steps.
- Chained scenarios can be factorized under certain conditions in a shorter series of conditions without respecting the original order of steps or variations.
- Factorization  keeps unchanged migration results only if the contact conditions are unchanged and if apparent partition coefficients are independent of temperature.
- Factorization should be avoided if the dispersion in food plays a significant role and in the presence of polar migrants.
- Factorized scenarios offer practical estimates for repeated uses.
- Probabilistic migration modeling (see §7) can be used to analyze the effects of known and unknown variations.

------

#### 4.2.4.1     Problem formulation

Mass transfer between components and materials occur insidiously along the supply chain. [Figure 14](#fig14) illustrates conditions triggering or altering migration from printed materials. Many uncontrolled factors may affect the extent of mass transfer: i) variable contact or exposure times, ii) random combinations of storage and transportation steps for intermediate, finished packaging materials and packaged foods, iii) changes in temperature and relative humidity (*e.g.*, seasonal, diurnal, international transportation), iv) modifications of boundary conditions during any stage of materials lifetime and product shelf-life. The redistribution of migrants between materials, layers, and components deserve special attention as it remains usually ignored by end-users. In practice, cross-contamination occurrences can be also considered indirectly (*i.e.*, without causaility) as impurities and non-intentionally added substances. Without being exhaustive, cross-contamination is highly likely from cured adhesives and printing inks, from recycled materials and any material with rich volatile organic compounds. Packaging and materials stored in stacks and reels ease cross contaminations by contacting internal and external layers, regardless the presence of a functional barrier (absolute such as a metallic layer or relative such as barrier polymer) in the structure. Due to periodic conditions, the inner layer can act as a reservoir of contaminants before the food is put in contact.

![fig14](./src/fig14.png)

<a name="figure14"></a>**Figure 14**. Illustration of the redistribution of the migrants from UV-cured printing ink and their subsequent migration in food for long shelf-life products. The depicted cases cover hot-filled products (*e.g.*, soups, pasteurized juices, sterilized dairy products), dry or ready-to-eat products stored in cardboard boxes.

 

From a mathematical point, the succession of steps and temperature variations can be seen as a sequence of constant conditions occurring in variable order. In the presence of ${{n}_{steps}}$ conditions occurring at the time $t=0,{{t}_{1,}}{{t}_{2}},...,{{t}_{{{n}_{steps}}}}$ , the composite solution is obtained by integrating the coupled system  via the Chasles’ relation:

$$
\begin{align}

  & C\left( x,{{t}_{{{n}_{steps}}}} \right)=\int\limits_{0}^{{{t}_{nsteps}}}{{{\left. \frac{\partial }{\partial t}C\left( x,t \right) \right|}_{C\left( x,t=0 \right)}}dt}=\sum\limits_{i=0}^{{{n}_{steps}}}{\int\limits_{{{t}_{i}}}^{{{t}_{i+1}}}{{{\left. \frac{\partial }{\partial t}C\left( x,t \right) \right|}_{C\left( x,t={{t}_{i}} \right)}}dt}} \\ 

 & {{C}_{F}}\left( {{t}_{{{n}_{steps}}}},C_{P}^{t=0},C_{F}^{t=0} \right)=\int\limits_{0}^{{{t}_{nsteps}}}{{{\left. \frac{d{{C}_{F}}\left( t,C(x,t) \right)}{dt} \right|}_{C_{P}^{t=0},C_{F}^{t=0}}}dt}=\sum\limits_{i=0}^{{{n}_{steps}}}{\int\limits_{{{t}_{i}}}^{{{t}_{i+1}}}{{{\left. \frac{d{{C}_{F}}\left( t,C\left( x,t \right) \right)}{dt} \right|}_{C_{P}^{t={{t}_{i}}},C_{F}^{t={{t}_{i}}}}}dt}}  

\end{align}
$$

The determinations of the duration of each step ${{\left\{ \Delta {{t}_{i}} \right\}}_{i={{t}_{i+1}}-{{t}_{i}}}}$ and their corresponding temperatures are critical. Representing all diurnal and seasonal temperature variations shown on the timeline of [Figure 14 ](#figure14) would require 2×450=900 successive simulations (one every 12 hours). On the one hand, a rigorous approach suggests that the congruence of all the steps should be strictly preserved to get reliable conclusions. In this case, how to identify the worst-case combination of conditions or steps? If the temperature variations are uncontrolled, how to build a conservative sequence? On the other hand, a naïve approach would suggest that the times series in Eq. XXX could be built from mutually independent steps assembled as the sum of a series of standard and well-controlled steps (*e.g.*, cold, moderate and warm days), and one series and of stochastic contributions, representing an extra safety margin. 

#### 4.2.4.2     A first intuitive approach

The variations of ${{C}_{F}}$ between steps are not factorizable but the dimensionless times are. Their effects are additive and commutative ${{C}_{F}}\left( \Delta F{{o}_{1}}+\Delta F{{o}_{2}} \right)={{C}_{F}}\left( \Delta F{{o}_{2}}+\Delta F{{o}_{1}} \right)$ (see Eq. XXX for the demonstration). If the diffusion coefficient is the only quantity varying with temperature and the plasticizer content (*i.e.*, no change in the partitioning, and the packaging dimensions), the effects of ${{n}_{steps}}$ is captured via the generalized Fourier number:

$$
Fo\left( {{t}_{nsteps}} \right)=\frac{\int\limits_{0}^{{{t}_{nsteps}}}{{{D}_{{{j}_{ref}}}}\left( T\left( t \right),\text{plasticizer}\left( t \right) \right)dt}}{l_{{{j}_{ref}}}^{2}}
$$

where ${{D}_{{{j}_{ref}}}}$ is the diffusion coefficient of the most limiting material, component or layer, and ${{l}_{{{j}_{ref}}}}$ the one-dimensional equivalent thickness (see Eq.  to identify ${{j}_{ref}}$ ).


Eq. XXX is equivalent to the low of the composition of velocities along a curvilinear coordinate system tangent to the trajectory going from A (${{\bar{v}}^{*}}\left( t=0 \right)=0$ ) to B (${{\bar{v}}^{*}}\left( t\to \infty  \right)=1$). The analogy between spatial translation along a winding road and the translation along the curve ${{\bar{v}}^{*}}$ vs. $Fo$ is illustrated in [Figure 15](#figure15) by comparing the travel via three modes of transport (each of duration ${{\left\{ \Delta {{t}_{i}} \right\}}_{i=1,2,3}}$at a speed: ${{\left\{ {{v}_{i}} \right\}}_{i=1,2,3}}$) with the cumulated contamination after three steps (of duration ${{\left\{ \Delta {{t}_{i}} \right\}}_{i=1,2,3}}$ at temperature ${{\left\{ {{T}_{i}} \right\}}_{i=1,2,3}}$. The total distance is ${{v}_{1}}\Delta {{t}_{1}}+{{v}_{2}}\Delta {{t}_{2}}+{{v}_{3}}\Delta {{t}_{3}}$ , and the cumulated migration is ${{\bar{v}}^{*}}\left( {\left( D\left( {{T}_{1}} \right)\Delta {{t}_{1}}+D\left( {{T}_{1}} \right)\Delta {{t}_{2}}+D\left( {{T}_{3}} \right)\Delta {{t}_{3}} \right)}/{{{l}^{2}}}\; \right)$ independently on the order of the steps. Replacing the physical time by a cumulant measuring the total arc-length of the curve or the road makes it possible to use endpoint estimates (one simulation) instead of chained simulations (${{n}_{steps}}$ simulations). In this new space, the residence times – represented by the density of markers on the curves – are not uniform. They are distributed more densely at departure and destination, but more sparsely in the middle region has the studied system follows different routes.


![fig15](./src/fig15.png)

<a name="figure15"></a>**Figure 15**. Illustration of the composition rules (a) for distances and (b) for the migration from a monolayer material, and of their invariance with the order of the steps (see Eqs.  and ).

 

Eq. XXX  suffers, however, a lack of generality as it applies only to the limiting mass transfer resistance and not to all layers. As a rule of thumb, It offers an acceptable if the function $\bar{v}*$ is monotone with $Fo$ (*i.e.*, ${{\bar{v}}^{*}}\left( F{{o}_{1}} \right)\le {{\bar{v}}^{*}}\left( F{{o}_{2}} \right)$ when $F{{o}_{1}}\le F{{o}_{2}}$. A simple counter-example can be, however, constructed by noticing that the concentration at equilibrium $C_{F}^{eq}$ depends only on the initial and final states, but not on intermediate steps. 

The conditions of exchangeability of steps (which is more generic than Eq. ) is discussed hereafter in more general terms. Two conditions are analyzed:  i) when the effect of the mass transfer resistance is considered (explicit representation) and when the number of molecules does not change, and ii) when an implicit food representation is used (*i.e.*, Eqs. XXX,XXX and XXX). The distinction between explicit and implicit food representation is relevant, as the boundary layer delays the effects of perturbations and may contain a significant amount of contaminants, which are ignored at low $Bi$ values in implicit representations (see §4.2.1.3).

#### 4.2.4.3     Strict conditions of exchangeability with explicit food models 

**Microscopic description of the random walk of molecules between P and F**

The visited distance by a deterministic system is $\ell =\Delta {{\ell }_{1}}+\Delta {{\ell }_{2}}+\Delta {{\ell }_{3}}+\cdots $ and is invariant with the order of visits. At microscopic scale (**i.e.**, at a scale where they can be separated), the trajectories of migrants in materials and the food verify this property (distances are additives), but with a different relationship with time (displacements are proportional with time). The random displacements of additives and residues shares, instead, notable features with random walks and continuous stochastic paths. In a loose sense, substances jump randomly at discrete times (random walks) or as continuous events (stochastic paths). Their skewed trajectories are nowhere differentiable, and velocities cannot be defined in a classical sense. Under a hypothesis of stationarity of the microscopic random process (the quantity $X\left( t \right)$ has the same statistics as $X\left( t+\epsilon  \right)$ for any $\epsilon $), we will justify a law of composition for the mean-square-distances, denoted $\left\langle {{\ell }^{2}} \right\rangle $, visited by the molecules: $\left\langle {{\ell }^{2}} \right\rangle $= ${{D}_{1}}\Delta {{t}_{1}}+{{D}_{2}}\Delta {{t}_{2}}+{{D}_{3}}\Delta {{t}_{3}}+\cdots $. $\left\langle {{\ell }^{2}} \right\rangle $ is mathematically defined as $\left\langle {{x}^{2}}\left( t \right) \right\rangle -{{\left\langle x\left( t \right) \right\rangle }^{2}}$, where $\left\langle x\left( t \right) \right\rangle $ is the average distance (first moment) and $\left\langle {{x}^{2}}\left( t \right) \right\rangle $the square distance to the food-packaging interface located at$x=0$:

$$
\begin{align}

  & \left\langle x\left( t \right) \right\rangle =\int\limits_{-\infty }^{+\infty }{x\rho \left( x,t \right)}dx \\ 

 & \left\langle {{x}^{2}}\left( t \right) \right\rangle =\int\limits_{-\infty }^{+\infty }{{{x}^{2}}\rho \left( x,t \right)}dx \\ 

\end{align}
$$

The one-dimensional space approximation with explicit food representation is defined on the domain $-{{l}_{F}}\le x\le {{l}_{p}}$ so that the probability density $\rho \left( x,t \right)$ verifies: $\int\limits_{-{{l}_{F}}}^{{{l}_{p}}}{\rho \left( x,t \right)dx}=1$ and $\rho \left( x,t \right)=0$ for $x<-{{l}_{F}}$ or $x>{{l}_{p}}$. A differential form of the growth of $\left\langle {{x}^{2}}\left( t \right) \right\rangle$ with time is inferred from a special case of the Fokker-Planck equation $\frac{\partial \rho \left( x,t \right)}{\partial t}=D\frac{{{\partial }^{2}}\rho \left( x,t \right)}{\partial {{x}^{2}}}$ (similar to Eq. XXX). By multiplying both sides by ${{x}^{2}}$ and by integrating over the entire domain, one gets:

$$
\int\limits_{-\infty }^{+\infty }{{{x}^{2}}\frac{\partial \rho \left( x,t \right)}{\partial t}dx}=\int\limits_{-\infty }^{+\infty }{{{x}^{2}}D\frac{{{\partial }^{2}}\rho \left( x,t \right)}{\partial {{x}^{2}}}dx}\ \ \text{with}\ D={{D}_{p}}\ \text{when}\ x\ge 0\ \text{and}\ D={{D}_{F}}\ \text{otherwise}
$$

with ${{D}_{p}}$ and ${{D}_{F}}$ being the diffusion coefficients in P and F, respectively. The left-hand side is equivalent to $\frac{\partial }{\partial t}\left\langle {{x}^{2}}\left( t \right) \right\rangle $. The right-hand side requires two successive integration by parts along with the impervious boundary conditions at $x=-{{l}_{F}}$ and  $x={{l}_{P}}$, and the conservation of the flux at $x=0$ . The simplifications associated with the compact support of $\rho \left( x,t \right)$ leads to:

$$
\begin{align}

  & \frac{\partial }{\partial t}\left\langle {{x}^{2}}\left( t \right) \right\rangle =\int\limits_{-\infty }^{+\infty }{D{{x}^{2}}\frac{{{\partial }^{2}}\rho \left( x,t \right)}{\partial {{x}^{2}}}dx}=\left[ {{D}_{F}}{{x}^{2}}\frac{\partial \rho \left( x,t \right)}{\partial x} \right]_{-\infty }^{0}+\left[ {{D}_{p}}{{x}^{2}}\frac{\partial \rho \left( x,t \right)}{\partial x} \right]_{0}^{+\infty } \\ 

 & -2{{D}_{F}}\int\limits_{-\infty }^{0}{x\frac{\partial \rho \left( x,t \right)}{\partial x}dx}-2{{D}_{P}}\int\limits_{0}^{+\infty }{x\frac{\partial \rho \left( x,t \right)}{\partial x}dx} \\ 

 & =0+0-2{{D}_{F}}\left( 0-\int\limits_{-\infty }^{0}{\rho \left( x,t \right)dx} \right)-2{{D}_{P}}\left( 0-\int\limits_{0}^{+\infty }{\rho \left( x,t \right)dx} \right) \\ 

 & =2{{D}_{F}}\int\limits_{-{{l}_{F}}}^{0}{\rho \left( x,t \right)dx}+2{{D}_{P}}\int\limits_{0}^{+{{l}_{p}}}{\rho \left( x,t \right)dx}=2{{D}_{P}}+2\left( {{D}_{F}}-{{D}_{P}} \right)\int\limits_{-{{l}_{F}}}^{0}{\rho \left( x,t \right)dx} \\ 

 & =2{{D}_{eff}}\left( t \right) \\ 

\end{align}
$$

For any initial distribution of substances $\rho \left( x,t=0 \right)$ (single molecule or collection of molecules), Eq.  describes the evolution of the mean-square distance $\left\langle {{x}^{2}}\left( t \right) \right\rangle$ to the FP interface. Since equations do not include any thermodynamical consideration (Eq.  was not enforced), Eq.  is valid while $\left\langle {{x}^{2}}\left( t \right) \right\rangle <{{\left( \frac{{{l}_{p}}+{{l}_{f}}}{2} \right)}^{2}}$ (beyond this length scale, finite size effects dominate and ${{D}_{eff}}\left( t \right)\to 0$). It shows that $\left\langle {{x}^{2}}\left( t \right) \right\rangle$ increases as $2{{D}_{p}}t$, when the concentration in F (${{C}_{F}}=\frac{{{l}_{p}}}{{{l}_{F}}}C_{P}^{0}\int\limits_{-{{l}_{F}}}^{0}{\rho \left( x,t \right)dx}$) is low. When the amount of substances in F becomes larger and because the molecules diffuse faster in F than in P, $\left\langle {{x}^{2}}\left( t \right) \right\rangle$ increases more rapidly. ${{D}_{eff}}\left( t \right)$ is the effective diffusion coefficient between P and F when both materials are replaced by an equivalent medium. A version similar to Eq.XX  for a homogeneous medium and known as Einstein equation is presented in §5.1.1 (see Eq. XXX).

 

**Condition of invariance of the dispersion of solutes with the properties of the contacting phase**

Eq. XXX shows that the composition  $\left\langle {{x}^{2}}\left( t \right) \right\rangle$ is independent of ${{D}_{F}}$ and of the amount already transferred to the food only if the inequality $\int\limits_{-\infty }^{0}{\rho \left( x,t \right)dx}$ $\ll$  $\frac{{{D}_{P}}}{{{D}_{F}}}\int\limits_{0}^{+\infty }{\rho \left( x,t \right)dx}$ is verified. The dispersion of contaminants is, hence, independent of the order of variations of ${{D}_{P}}$ with time when:

$$
{{C}_{F}}\left( t \right)\ll \frac{{{l}_{p}}}{{{l}_{F}}}\frac{{{D}_{P}}}{{{D}_{P}}+{{D}_{F}}}C_{P}^{0}=\frac{{{D}_{P}}}{{{D}_{P}}+{{D}_{F}}}{{L}_{P/F}}C_{P}^{0}\approx \frac{{{D}_{P}}}{{{D}_{F}}}{{L}_{P/F}}C_{P}^{0}
$$

By noticing that ${{L}_{P/F}}C_{P}^{0}$ is the total migration (see Eq. ) and that the diffusion in polymers is 1:100 or less lower than in the food one, we see that the invariance with the order of the steps hold only at the beginning of the migration process or when the chemical affinity for the food is very low. Eq.  could be justified with the example of a large food volume submitted to a cold condition during $\Delta {{t}_{0}}$, denoted $\left( \Delta {{t}_{0}},{{T}_{1}} \right)$ followed by a warm condition during the same time, denoted $\left( {{t}_{0}},{{T}_{2}} \right)$ with ${{T}_{2}}>{{T}_{1}}$ . The normal order $\left( \Delta {{t}_{0}},{{T}_{1}} \right)\times \left( \Delta {{t}_{0}},{{T}_{2}} \right)$  would lead to a small mass transfer during $\left( \Delta {{t}_{0}},{{T}_{1}} \right)$ and a very strong during $\left( \Delta {{t}_{0}},{{T}_{2}} \right)$. Performing the transfer in the reverse order $\left( \Delta {{t}_{0}},{{T}_{2}} \right)\times \left( \Delta {{t}_{0}},{{T}_{1}} \right)$ will cause even a higher mass transfer during the first step. If the food is large, the extra number of molecules transferred to the food during the first step $\left( \Delta {{t}_{0}},{{T}_{2}} \right)$ will not have enough time to be transferred back to the packaging. As the food-packaging contact is not symmetric ${{D}_{F}}\gg {{D}_{P}}$ and ${{l}_{F}}\gg {{l}_{p}}$, the two orders might not lead to similar irreversible behaviors: the food is contaminated in both cases but not in the same extent.

 

**Discussion on the limits introduced by implicit models**

The fundamental results exposed here rely on an explicit representation of the food, where molecules displace at a finite velocity. This subtle detail is not reproduced in implicit representations, which assume a perfect mixing outside the mass transfer boundary layer (the velocity of molecules). Only a delay is considered in the boundary assuming a linear profile instead of a parabolic one (see the distribution of molecules depicted in green symbols [Figure 7](#figure7)). The next paragraph reviews the conditions of commutativity of implicit models under variable conditions. The condition of commutativity is less severe as the back flux from the contacted phase is immediately compensated in the numerical scheme. But as shown here, the condition of commutativity is expected to be verified in real life, only when the condition  is met. The equivalence between time and temperature is acceptable only far from the equilibrium (${{\bar{v}}^{*}}\to 0$ ), but it is unacceptable closer to equilibrium (${{\bar{v}}^{*}}\to 1$ ), where the effective mass transfer is governed by effective partition coefficient across the P-F interface. From Eq. X, partitioning coefficients are independent of the temperature only if the difference of free sorption energies between P and F are kept constant (${{G}_{P}}={{G}_{F}}$) when the temperature is raised.

 

#### 4.2.4.4     Conditions of exchangeability in food implicit models

------

Overview of implicit numerical models and their solutions**

Food implicit models are by far the most used. They are more flexible to accommodate variable conditions and chained conditions. They have been implemented with various numerical models using different spatial discretization schemes. The finite difference method is the dominant approach in one-dimension problems, but it loses accuracy at interfaces when large jumps in concentrations and diffusion coefficients occur. The finite element method is the standard in the industry as it enables to integrate any partial differential equation system on arbitrary geometrical domains, using a grid approximation (consisting of triangles, quadrangles, curvilinear polygons). The finite volume method is in essence similar (values are calculated on a meshed geometry), but the equations are integrated on small, but not infinitesimal, volumes. By positioning the interface between volumes at the exact location, where thermodynamic constraints such as Eq. XXX needs to be strictly verified, the method enables to maintain exact mass balances and the continuity of chemical potentials between materials. The pros and cons of each method are discussed in Ref.[181](#_ENREF_181). The three methods can be put in a matrix form, coding for a system of ordinary differential equation:

$$
\frac{\partial \underline{\mathbf{C}}\left( t \right)}{\partial t}=\underline{\underline{\mathbf{M}}}\underline{\mathbf{C}}\left( t \right)
$$

The common practice is to include the concentration in P discretized ${{n}_{nodes}}$ and the concentration in F in $\underline{\mathbf{C}}\left( t \right)$. Since the food is represented implicitly, one node is sufficient for the food and $\underline{\mathbf{C}}\left( t \right)$is a $\left( {{n}_{nodes}}+1 \right)\times 1$ vector mapping its continuous version$C\left( x,t \right)$. $\underline{\underline{\mathbf{M}}}$ is a triband matrix $\left( {{n}_{nodes}}+1 \right)\times \left( {{n}_{nodes}}+1 \right)$ for discretization scheme at order 1 and pentaband matrix for quadratic finite elements techniques.

When the transport and thermodynamic properties are constant, the solution of Eq. XXX respectively to the initial condition $\underline{\mathbf{C}}\left( t \right)=0$ is $\exp \left( -\underline{\underline{\mathbf{M}}}t \right)\underline{\mathbf{C}}\left( t=0 \right)$, with $\exp \left( -\underline{\underline{\mathbf{M}}}t \right)=\sum\limits_{k=0}^{\infty }{\frac{{{\left( -t \right)}^{k}}}{k!}{{\underline{\underline{\mathbf{M}}}}^{k}}}$. 

 

**Composition rules when chained simulations are used (example with three steps)**

The solution of the mass transfer associated with three conditions: $\left( {{\underline{\underline{\mathbf{M}}}}_{1}},\Delta {{t}_{1}} \right)\times \left( {{\underline{\underline{\mathbf{M}}}}_{2}},\Delta {{t}_{2}} \right)\times \left( {{\underline{\underline{\mathbf{M}}}}_{3}},\Delta {{t}_{3}} \right)$, with ${{\left\{ {{\underline{\underline{\mathbf{M}}}}_{i}} \right\}}_{i=1..3}}$ the stiffness matrix for the *i*th step (*e.g.*, coding the effect of temperature on diffusion and partition coefficients) reads:

$$
\underline{\mathbf{C}}\left( \Delta {{t}_{1}}+\Delta {{t}_{2}}+\Delta {{t}_{3}} \right)=\exp \left( -{{\underline{\underline{\mathbf{M}}}}_{3}}\Delta {{t}_{3}} \right)\exp \left( -{{\underline{\underline{\mathbf{M}}}}_{2}}\Delta {{t}_{2}} \right)\exp \left( -{{\underline{\underline{\mathbf{M}}}}_{1}}\Delta {{t}_{1}} \right)\underline{\mathbf{C}}\left( t=0 \right)
$$

The steps are exchangeable if the equality  satisfies also:

$$
\underline{\mathbf{C}}\left( \Delta {{t}_{1}}+\Delta {{t}_{2}}+\Delta {{t}_{3}} \right)=\exp \left( -{{\underline{\underline{\mathbf{M}}}}_{1}}\Delta {{t}_{1}}-{{\underline{\underline{\mathbf{M}}}}_{2}}\Delta {{t}_{2}}-{{\underline{\underline{\mathbf{M}}}}_{3}}\Delta {{t}_{3}} \right)\underline{\mathbf{C}}\left( t=0 \right)
$$

which is verified only if ${{\left\{ {{\underline{\underline{\mathbf{M}}}}_{i}} \right\}}_{i=1..3,i\ne j}}$ and ${{\left\{ {{\underline{\underline{\mathbf{M}}}}_{j}} \right\}}_{j=1..3,i\ne j}}$ commutes (meaning that ${{\underline{\underline{\mathbf{M}}}}_{i}}{{\underline{\underline{\mathbf{M}}}}_{j}}={{\underline{\underline{\mathbf{M}}}}_{j}}{{\underline{\underline{\mathbf{M}}}}_{i}}$ for $i\ne j$).  

 

**Conditions of exchangeability imposed by the physical-chemistry**

Reference [181](#_ENREF_181) demonstrated that a necessary and sufficient condition to have the solution  applicable is that the ratios of Henry coefficients ${{\left\{ \frac{{{k}_{u}}\left( t \right)}{{{k}_{v}}\left( t \right)} \right\}}_{u,v=0...m,u\ne v}}$ (see their definitions in Eq. XXX) remain constant with time between all considered steps. As a result, adding or removing a material/layer/food (*i.e.*, changing  ${{k}_{u}}$ from $\infty $ to 0 or the reverse) breaks the condition of exchangeability of steps. The condition of exchangeability is also likely to be lost for polar solutes dispersed between polar and apolar phases. As shown in §6.2.2 (see Eq. XXX at infinite dilution when ${{\phi }_{i}}\to 0$ ), the excess enthalpies of mixing are negative for the firsts indeed and positive for the lasts. In this particular case, the apparent activation energy of the partition coefficient between $u$ and $v$ is maximal, and the effect of temperature needs to be simulated by respecting the order of the temperature variations (*i.e.*, by using Eq.  XXX of Eq. XXX). For apolar solutes, the activation of ${{k}_{u}}$ is, conversely, almost compensated by a symmetric variation of ${{k}_{v}}$ , when the temperature is changed.

 

### 4.2.5     Discussion on the choice of accelerated conditions and the identification of critical steps

------
**Highlights for multiple steps and variable conditions**

- Accelerated conditions and time-temperature relationships need to be verified according to the rules of factorization mentioned in §4.2.4.
- The contribution of any step or component by obtained by difference if the comparison includes all steps and components.
- Probabilistic migration modeling (see section 7) can be used to evaluate the safety margin associated with accelerated and equivalent conditions.

------

#### 4.2.5.1     Factors of acceleration and possible biases. 

The different theoretical developments presented in this part highlighted the complications to reach exact modeling of migration for a given product. The real contact conditions and couple time×temperature need to be considered and incorporated in a proper simulation scenario. Subtle variations of temperature and contact times may modify the estimates. We remember the reader that a tiered approach, as shown in Figure 5 is always preferable. A robust approach to assess the effects of various sources of variability (*e.g.*, temperature variations) and uncertainty (*e.g.*, effects of temperature of diffusion and partition coefficients) is presented in §7. Probabilistic modeling applied to Eq. XXX can be used to derive conservative estimates either with a prescribed risk or with a prescribed safety margin (see. [Figure 23](#figure23)).

Accelerated tests are used experimentally to reproduce the level of contamination after the end of the food product shelf-life while using much shorter contact times. By assuming that migration is controlled by one single mass transfer resistance (index ${{j}_{ref}}$ ), the principle relies on replacing the real condition $\left( t_{life}^{shelf},T_{life}^{shelf} \right)$ by the accelerated condition $\left( t_{life}^{shelf},T_{conditions}^{accelerated} \right)$. When amount transferred at equilibrium is not affected by temperature (see. §6.2.2  for discussion), the ratio of ${{C}_{F}}$ values between accelerated and real conditions read:

$$
\begin{align}

  & \frac{{{C}_{F}}\left( t_{conditions}^{accelerated},T_{conditions}^{accelerated} \right)}{{{C}_{F}}\left( t_{life}^{shelf},T_{life}^{shelf} \right)}={{{\bar{v}}}^{*}}\left( \frac{{{D}_{{{j}_{ref}}}}\left( T_{conditions}^{accelerated} \right)}{{{D}_{{{j}_{ref}}}}\left( T_{life}^{shelf} \right)}\frac{t_{conditions}^{accelerated}}{t_{life}^{shelf}} \right) \\ 

 & ={{{\bar{v}}}^{*}}\left( \exp \left[ -\frac{E_{a}^{{{D}_{{{j}_{ref}}}}}}{R}\left( \frac{1}{T_{conditions}^{accelerated}}-\frac{1}{T_{life}^{shelf}} \right) \right]\frac{t_{conditions}^{accelerated}}{t_{life}^{shelf}} \right) \\ 

 & \le \sqrt{\exp \left[ -\frac{E_{a}^{{{D}_{{{j}_{ref}}}}}}{R}\left( \frac{1}{T_{conditions}^{accelerated}}-\frac{1}{T_{life}^{shelf}} \right) \right]\frac{t_{conditions}^{accelerated}}{t_{life}^{shelf}}}\ \ \text{for}\ \text{monolayer}\ \text{materials} \\ 

\end{align}
$$

where $E_{a}^{{{D}_{{{j}_{ref}}}}}$ is the apparent activation energy of the diffusion. This parameter is critical and will be discussed in §5.4. For monolayer materials, an activation energy of 60 kJ⋅mol-1 yields an activation factor of ca. 5 when the temperature is raised to 40°C from 20°C. Testing food packaging at 40°C during ten days is equivalent to a test of 50 days at 20°C. If the shelf-life is of 100 days, the risk of underestimation of the true migration can reach up to factor ${\sqrt{2}}/{2}\;\approx 0.7$. In the presence of a functional barrier, the factor of underestimation can be even larger as shown in [Figure 10](#figure10). 

 

#### 4.2.5.2     Causality, critical steps, and crtical components. 

In the presence of multiple steps, it is recommended to perform a test, which mimics the most critical steps for the final value of ${{C}_{F}}$. Due to the lack of generality of Eqs. XXX  and XXX, the factorization strategies should be avoided to identify critical steps. The reference [181](#_ENREF_181) demonstrated that the contribution of one single step, denoted $\Delta {{C}_{F}}\left( \text{step}\ i \right)$ ,  or of one component (label, cap, film, bottle, etc.), denoted $\Delta {{C}_{F}}\left( \text{component}\ i \right)$ ,  could always be determined by difference as:

$$
\begin{align}

  & \Delta {{C}_{F}}\left( \text{step}\ i \right)=\max \left( {{\left. {{C}_{F}} \right|}_{1\to 2\to \cdots \to {{n}_{\text{steps}}}}}-{{\left. {{C}_{F}} \right|}_{{1\to 2\to \cdots \to {{n}_{\text{steps}}}}/{i}\;}},{{\left. {{C}_{F}} \right|}_{\text{step}\ i}} \right) \\ 

 & \Delta {{C}_{F}}\left( \text{component}\ i \right)=\max \left( \left. {{C}_{F}} \right|_{1+2+\cdots +{{n}_{components}}}^{assembly}-\left. {{C}_{F}} \right|_{1+2+\cdots +{{n}_{\text{components}}}-i}^{assembly},{{\left. {{C}_{F}} \right|}_{\text{component}\ i}} \right) \\ 

\end{align}
$$

where ${1\to 2\to \cdots \to {{n}_{\text{steps}}}}/{i}\;$ represents the sequence including all steps except the step $i$; $1+2+\cdots +{{n}_{\text{components}}}-i$ including all components except the component $i$ .

The contributions of the sum of all steps and components are larger than the real contamination. This behavior is expected as several steps can lead to similar effects. Corollary,  Eq. XXX enables the identification of steps, which occur without contact with food (see [Figure 14](#figure14)) or of any component, even if there is no direct or intended contact with food. When the contributions $\Delta {{C}_{F}}$  are normalized by a threshold of concern, $SML$, they can be used to seek the likely causes (industrial practices or components) which lead or which could lead ${{C}_{F}}$ to exceeding tolerances. The effects can be cumulated over several substances to get a global criticality of the step or the component of interest. This procedure has been used as the foundation of a preventive approach of mass transfer reusing the principles of the methodology FMECA (Failure Mode Effects and Criticality Analysis) [181](#_ENREF_181). The approach can be employed to analyze an entire supply chain or industrial practices, as illustrated in [Figure 13](#figure13). The assessment of the risk of cross-contamination of the constituents of UV-cured printing ink during the long-term storage of pots in stacks is shown in [Figure 16](#figure16). The stacking step is indexed 1.  The redistribution of contaminants triggers contamination immediately during and after filling the pot at 80°C. Migrants accumulates continuously in the food during storage until the final consumption. By difference, the contribution of the first step can be evaluated. After normalization by $SML$, it is demonstrated that the first step provides the largest contribution to the risk of exceeding the legal limit. 

![fig16](./src/fig16.png)

<a name="figure16"></a>**Figure 16**. Contribution of the stacking of pots before filling them at hot temperatures with a Chinese soup along the supply chain: a. steps of supply chain; b. assessment of a severity of a single step; c. comparison of a several step. After ref. [181](#_ENREF_181).



# 5     Diffusion properties in polymers

------
**Highlights on the diffusion coefficients in polymers (thermoplastics, thermosets, elastomers)**

- Diffusion coefficients are essential properties to evaluate the performances of functional barriers and to assess the kinetics of migration in food.
- Diffusion coefficients are strongly affected by the size of migrants: volume effects dominate for bulk and rigid solutes (*e.g.*, aromatic solutes), the number of rigid sub-units dominates in flexible solutes.
- The activation of diffusion follows an Arrhenius relationship only far from the glass transition temperature ( ). Near  , free-volume effects dominate, and apparent activation energies depend on both temperature and the size of diffusants.
- For compliance testing, it is recommended to use overestimates of diffusion coefficients to reduce the risk of underestimation of migration.
- When modeling is used to validate mechanical recycling processes of polymers, it is recommended to use either realistic estimates or probabilistic modeling (see section 7).

------

  

## 5.1    Definitions of diffusion coefficients

Diffusion coefficients in the polymer, ${{D}_{P}}$ , are essential properties to calculate migration according to Eqs.  and . Reference textbooks [see chapter 11 of Ref.[102](#_ENREF_102), see Ref.[107](#_ENREF_107) proposed several definitions of diffusion coefficients including overestimates (see chapters in Ref.[89](#_ENREF_89)). We review in this section the definitions, which are relevant in macroscopic migration models and which are capable of incorporating the effects of the molecular structure of the migrants and the polymer.

### 5.1.1     Self- and trace-diffusion coefficients

As a starting point, it is worth noticing that the random walk of migrants (without gradient) does not occur in an empty space but among other molecules (polymer segments and other solutes, such as plasticizers). By reusing the illustrations of [Figure 7](#figure7), the net flux density perpendicular to a cross section located at a position $x$ is proportional to the difference of the net velocities between the solutes (*i.e.*, molecules of the migrant of interest), $u\left( x,t \right)$ and the Stokesian velocity of the surrounding molecules ${{u}_{0}}\left( x,t \right)$ :

$$
[J\left( x,t \right)=C\left( x,t \right)\left[ u\left( x,t \right)-{{u}_{0}}\left( x,t \right) \right]=-{{D}_{P}}{{\left. \frac{\partial C\left( x,t \right)}{\partial x} \right|}_{x}}
$$

In the limit case where the migrant is of the same nature as the surrounding molecules (*e.g.*, pure liquid plasticizers, pure solvents), $u\left( x,t \right)={{u}_{0}}\left( x,t \right)$, there is no net flux ($J\left( x,t \right)=0$) and no macroscopic gradient ${{\left. \frac{\partial C\left( x,t \right)}{\partial x} \right|}_{x}}=0$. The diffusion coefficient is still defined, and it is called the self-diffusion coefficient but not correlated to any macroscopic gradients.

The opposite limiting case corresponding to infinite dilution ($C\left( x,t \right)\to 0$ ) and a solid behavior (${{u}_{0}}\left( x,t \right)\to 0$ ) also leads to an extremely low net flux ($J\left( x,t \right)\to 0$), but with a concentration gradient. The corresponding diffusion coefficient is called trace diffusion coefficient. Experimentally, the velocity $u\left( x,t \right)$can be determined from the lag-time, ${{t}_{lag}}$ , associated with the migration across a functional barrier of thickness ${{l}_{fb}}$. By noticing that ${{t}_{lag}}=\tfrac{l_{fb}^{2}}{6{{D}_{P}}}$ (see §4.2.3.5), one gets

$$
u\left( x=0,t={{t}_{lag}} \right)=\frac{J\left( x=0,t={{t}_{lag}} \right)}{C\left( x=0,t={{t}_{lag}} \right)}=\frac{{{l}_{fb}}}{{{t}_{lag}}}=6\frac{{{l}_{fb}}}{l_{fb}^{2}}{{D}_{P}}=6\frac{{{D}_{P}}}{{{l}_{fb}}}
$$

Eq. XXX is correct but suffer a lack of generality for an arbitrary initial distribution of solutes. The approach presented in §4.2.4.3 is more general. Eq  associated to${{D}_{F}}={{D}_{P}}$ shows that:

$$
{{D}_{P}}=\frac{1}{6{{N}_{migrants}}}\underset{t\to \infty }{\mathop{\lim }}\,\frac{d}{dt}\left\langle \sum\limits_{i=1}^{{{N}_{migrants}}}{{{\left\| \mathbf{r}_{i}^{CM}\left( t \right)-\mathbf{r}_{i}^{CM}\left( 0 \right) \right\|}^{2}}} \right\rangle =\frac{1}{6}\underset{t\to \infty }{\mathop{\lim }}\,\frac{d}{dt}{{g}_{CM}}\left( t \right)\approx \frac{{{g}_{CM}}\left( t \right)}{6t}
$$

Factor 6 appears in Eqs. XXX-XXX instead of factor 2 shown in Eq , because the random walks are considered in three dimensions and not anymore in one dimension. The mean-square-displacement ${{g}_{CM}}\left( t \right)={\left\langle \sum\limits_{i=1}^{{{N}_{migrants}}}{{{\left\| \mathbf{r}_{i}^{CM}\left( t \right)-\mathbf{r}_{i}^{CM}\left( t=0 \right) \right\|}^{2}}} \right\rangle }/{{{N}_{migrants}}}\;$ is calculated accordingly from the 3D positions of the center-of-mas of all migrants, ${{\left\{ \mathbf{r}_{i}^{CM}\left( t \right) \right\}}_{i=1..{{N}_{migrants}}}}$ and averaged over all possible initial positions ${{\left\{ \mathbf{r}_{i}^{CM}\left( t=0 \right) \right\}}_{i=1..{{N}_{migrants}}}}$. The definition  is used to calculate diffusion coefficients by molecular dynamics simulations in polymers [140](#_ENREF_140), [182](#_ENREF_182). In practice and to remove any uncontrolled drift (**i.e.*,* $\left| {{u}_{0}} \right|>0$ ) due to the displacement of surrounding molecules, the positions are not absolute but taken respectively to the center-of-mass of the whole system.

 

### 5.1.2     Mutual diffusion coefficients

Section 5.1.1 describes diffusion either in a stream or when there is only one diffusing species. Diffusion in a medium with variable composition generates fluctuations in ${{u}_{0}}$ (see Eq. ) and change the nature of the interactions with the surrounding molecules. The picture is complete if we also consider the point of view of the surrounding molecules, whose diffusion is also affected. An example could be the diffusion of antioxidant (or any large molecule) in a heterogeneously plasticized polymer (*i.e.*, see Ref. [183](#_ENREF_183)), where the local mobility of the antioxidant is strongly enhanced by the local amount of plasticizer in the ternary system (migrant+plasticizer+polymer). Similarly, the diffusion of the plasticizer increases the mobility of polymer segments. Strictly speaking, all these effects cannot be described by the sole shift of the glass transition temperature and the increase of free volumes presented in §5.3 (see also discussion in chapters 5 and 6 of the reference textbook [107](#_ENREF_107)). As the diffusion of one species affects the diffusion of all other species and reciprocally, self- and trace diffusion coefficients cannot be used. The correct description requires the use of the Onsager’s theory and generalized forces. The entire system (polymer+solutes) reorganizes to minimize its total free-energy. A cationic surfactant will, for example, accumulate at the surface of the material to minimize its interaction with the polymer. A poorly soluble colorant or pigment will do the same at high concentration. The resulting concentration profiles evolve therefore spontaneously from an initial uniform distribution after processing to a highly heterogeneous one. Such evolutions with matter moving from low to high concentrated regions cannot be predicted with concentration gradients as effective driving forces. The gradients of chemical potentials need to be used as an effective driving force. In the sake of simplicity and to avoid a tensor definition of the diffusion coefficient, we consider only the effect of the chemical potential of the migrant, denoted $\mu \left( x,t \right)$ , as driving force. Rigorously, a linear relationship should be considered between the flux of the migrating substance and the driving forces associated with all species in the mixture. The generalized driving force, $f$ , induced by the local variation of the chemical potential of the migrant in the mixture reads: 

$$
f=-{{\left. \frac{\partial \mu \left( x,t \right)}{\partial x} \right|}_{t}}={{\zeta }_{mutual}}\left[ u\left( x,t \right)-{{u}_{0}}\left( x,t \right) \right]={{\zeta }_{mutual}}\frac{J\left( x,t \right)}{C\left( x,t \right)}
$$

where ${{\zeta }_{0}}$ is a friction coefficient; $\mu \left( x,t \right)={{\mu }_{0}}+RT\ln \left( {{\gamma }^{v}}\phi  \right)$ is the chemical potential of the migrant defined respectively to its volume fraction $\phi $ and activity coefficient ${{\gamma }^{v}}$. The relationship between the gradient of chemical potential and the concentration gradient is given by:

$$
\begin{align}

  & {{\left. \frac{\partial \mu \left( x,t \right)}{\partial x} \right|}_{t}}={{\left. \frac{\partial \mu \left( x,t \right)}{\partial \phi } \right|}_{t}}{{\left. \frac{\partial \phi \left( x,t \right)}{\partial x} \right|}_{t}}=RT{{\left. \frac{\partial \ln \left( {{\gamma }^{v}}\phi \left( x,t \right) \right)}{\partial \phi } \right|}_{t}}{{\left. \frac{\partial \phi \left( x,t \right)}{\partial x} \right|}_{t}} \\ 

 & =\frac{RT}{C\left( x,t \right)}\left[ {{\left. \frac{\partial \ln \left( {{\gamma }^{v}} \right)}{\partial \ln \left( \phi \left( x,t \right) \right)} \right|}_{t}}+1 \right]{{\left. \frac{\partial C\left( x,t \right)}{\partial x} \right|}_{t}}=\Gamma \frac{RT}{C\left( x,t \right)}{{\left. \frac{\partial C\left( x,t \right)}{\partial x} \right|}_{t}}  

\end{align}
$$

with $\Gamma \left( \phi  \right)=1+{{\left. \frac{\partial \ln \left( {{\gamma }^{v}} \right)}{\partial \ln \left( \phi \left( x,t \right) \right)} \right|}_{t}}$ being the thermodynamical factor depending on the composition, whose dependence with concentration can be determined from the expression of the sorption isotherm see §6.2). An effective diffusion coefficient, $D_{P}^{mutual}$ ,  related to the flux $J\left( x,t \right)=-D_{P}^{mutual}{{\left. \frac{\partial C\left( x,t \right)}{\partial x} \right|}_{t}}$ is identified from Eqs. XXX and XXX as:

$$
D_{P}^{mutual}\left( \phi  \right)=\Gamma \left( \phi  \right)\frac{RT}{{{\zeta }_{mutual}}\left( \phi  \right)}=\Gamma \left( \phi  \right)\frac{{{\zeta }_{trace}}}{{{\zeta }_{mutual}}\left( \phi  \right)}D_{P}^{trace}
$$

where subscripts *trace* and *mutual* refer to the value of the property at infinite dilution and in mixture, respectively. For binary solvent-polymer mixtures, Vrentas and Vrentas [184](#_ENREF_184) proposed the following evolution of the friction coefficient:

$$
\frac{{{\zeta }_{mutual}}\left( \phi  \right)}{{{\zeta }_{trace}}}=\alpha {{\phi }^{2}}+\left( 1-\phi  \right)\left( 1+2\phi  \right)
$$

with $\alpha =\frac{{{\overline{V}}_{\begin{smallmatrix} pure solvent\end{smallmatrix}}}}{{{\overline{V}}_{\begin{smallmatrix} dissolved polymer\end{smallmatrix}}}}\frac{D_{solvent}^{self}}{D_{\begin{smallmatrix} dissolved polymer \end{smallmatrix}}^{trace}}$ and ${{\overline{V}}_{X}}$ the molar volume of X. Eq.  is not thought to be general for non-solvent or not plasticizing molecules such as hindered phenols or aromatic amines. The fundamental reason is that such migrants are crystalline (solids) at high concentration. 

It is worth noticing that $D_{P}^{trace}$ can contain additional concentration dependence due to the extra free-volumes brought by the solute itself and the possible shift in the glass transition temperature. For most applications and in the absence of a unified theory, the equation of Fujita can be applied to describe the concentration dependence:

$$
\ln \frac{D_{P}^{mutual}\left( C \right)}{D_{P}^{mutual}\left( {{C}_{ref}} \right)}=-\beta \left( 1-\frac{C}{{{C}_{ref}}} \right)
$$

with $\beta $ a concentration to be determined experimentally.

 

## 5.2   Effect of the geometry of migrants on $D_p$ values

Migrants from polymer are not gas molecules, such as water, dioxygen, carbon dioxide, but heavy molecules larger than voids in the polymer. The smallest additives are monomer residues and solvents commensurable to one or several monomers. At the concentration of use of most of additives and residues (*i.e.*, except plasticizers which obey to mutual diffusion), the many pair contacts between the segments of the polymer and the migrant control the rate of the translation of the most mobile species (the migrant). This configuration is very different from the situation in the food or in liquid food simulants, where the food constituents (water, oil, ethanol, etc.) are much smaller than polymer chains and are usually packed less densely. The correct picture is to consider that the trace diffusion of the migrant in the polymer is smaller than the self-diffusion of food constituents (typically 10^-9^-10^-10^ m^2^⋅s^-1^), but much larger than the self-diffusion of the polymer itself (<10^-22^ m^2^⋅s^-1^).

The exact mechanism of translation of molecules larger than voids in solid polymers has not been fully elucidated yet. They have been investigated independently by two communities: the community of free-volume-theory was interested in the mutual-diffusion of polymer solvents whereas the food packaging community was focused on the development of practical overestimate models for compliance testing. Due to the different nature of the considered substances and the type of considered diffusion coefficients (mutual diffusion measured in ^1^H spin-echo nuclear magnetic resonance and trace diffusion coefficients measured via desorption kinetics), the interactions between the two communities have been limited. Based on the work [136](#_ENREF_136), [185](#_ENREF_185), a unified approach has been sketched and is summarized in [Figure 17](#figure17). Early works [138](#_ENREF_138), [186-189](#_ENREF_186) suggested that the volume of the entire molecule or its rough estimate the molecular mass (see the discussion in Ref.[190](#_ENREF_190)) offered a proper scaling parameter of diffusivities, at least, at the first tier. The use of the entire mass or volume of the migrant is misleading as it covers very different realities for rigid and flexible solutes. Flexible migrants, possibly with large masses, can benefit from the translation of smaller rigid units; whereas a large rigid migrant rarely requires free-volumes matching the shape and size of a larger rigid block. Without paying attention, both the effects of block sizes and their numbers combine, giving an apparent correlation of the logarithm of ${{D}_{P}}$ with the logarithm of the molecular mass $M$(concerted displacements of rigid units or blocks) and with the total volume of the migrant ${{n}_{\begin{smallmatrix} rigid block\end{smallmatrix}}}\cdot V_{\begin{smallmatrix} rigid block\end{smallmatrix}}^{_{vdW}}$. In this description, short *n*-alkanes are partly rigid since the constraints of torsion prevent blocks from translating independently.

<img src="./src/fig17.png" alt="fig17" style="zoom: 67%;" />

<a name="figure17"></a>**Figure 17**. Scaling of Diffusion coefficients between rigid and connected blocks with molar mass and Van der Waals volume in a thermoplastic polymer (groups A and B refer to substances defined in Figure 18).

 

Figure 18 reports the diffusion coefficients of 49 substances (105 ${{D}_{P}}$ values) in low-density polyethylene (LDPE) at 23°C collected by the National Institute of Standards (NIST)[191](#_ENREF_191)  and the European Commission (EC) [192](#_ENREF_192) from literature. The data were selected to reproduce the main features of Figure 17: the decrease of ${{D}_{P}}$ values with the number of rigid blocks and the almost invariance of ${{D}_{P}}$ with the exact shape and flexibility of the solute. The concept of invariance considers here that there is no visible effect if the induced variation is lower than the uncertainty associated with the measurement of ${{D}_{P}}$. As an illustration of Ref.[138](#_ENREF_138) showed that the group of substances denoted A and B presented differences in diffusivities, which could not be explained by the differences between molecules. Three parameters were considered: the van-der-Walls volume (${{V}_{vdW}}$ in  &#8491;^3^), the gyration radius and the shape factor. Linear alkanes scaled as a power law with the molecular mass $M$, ${{D}_{p}}\propto {{M}^{-\alpha }}$ with $\alpha $ of 1.58 (range: 1.09-1.68) for [1] NIST data and 1.84 (range: 0.75-2.01) for [2] EC data. For bulky solutes, the large uncertainty in the experimental ${{D}_{P}}$ values suggests that diffusivities can be considered significantly different only if the ${{V}_{vdW}}$ differences between molecules are larger than 35    3. The correlation with volume is consistent with the free-volume theory of diffusion [136](#_ENREF_136), and as ${{V}_{vdW}}$is also correlated with the number of heavy atoms (see Figure 18e, the correlation depends on the type of substances), a reasonable upper envelope of ${{D}_{P}}$values has been proposed to overestimate diffusivities[185](#_ENREF_185). The equation is coined Equation of *Piringer*:

$$
\ln D=A{{'}_{P}}-0.1351{{M}^{2/3}}+0.003M-\frac{\tau +10454}{RT}
$$
with the key parameters $A{{'}_{P}}=11.5$ and $\tau =0$ for LDPE. The values for other polymers are reported in the EU report [143](#_ENREF_143). The model of *Piringer* provides only a variable overestimation ranging from a factor 0.63 (underestimation) to 100 (see the inset [Figure 18a](#figure18a)).

<img src="./src/fig18.png" alt="fig18" style="zoom:67%;" />

<a name="figure18"></a>**Figure 18**. Scaling of diffusion coefficients of 49 substances (n-alkanes, two groups of molecules A and B with similar     values) in LDPE at 23°C with molar mass, M, and the van-der-Walls volume. Data from: [1] Ref.[193](#_ENREF_193), [2] Ref.[191](#_ENREF_191), [192](#_ENREF_192).

1:
methane [1]; 2: ethane [1]; 3: propane [1]; 4: n-pentane [1,2]; 5: n-hexane
[1,2]; 6: n-heptane [1,2]; 7: n-octane [1,2]; 8: n-decane [1,2,A]; 9:
n-octadecane [1]; 10: n-dodecane [2,A]; 11:
2-trans-3,7-dimethyl-2,6-octadien-8-ole (geraniol) [A]; 12:
3,7-dimethyl-6-octen-1-ol (citronellol) [A]; 13: n-decylaldehyde or n-decanal
(aldehyd c10) [A]; 14: 3,7-dimethyl-1-octanol [A]; 15: decylalcohol or
1-decanol [A]; 16: cis-undecen-8-al (aldehyd c11 inter) [A]; 17: n-undecen-2-al
(aldehyd c11) (2-undecenal) [A]; 18: n-undecylaldehyde (aldehyde c11) [A]; 19:
ethyloctanoate [A]; 20: 2-methoxy-4-propenylanisol (methylisoeugenol) [A]; 21:
citronellyl formate or 6-octen-1-ol, 3,7-dimethyl-, formate [A]; 22:
2-methyl-3-(4-isopropyl)phenylpropanal (cyclamen aldehyde) [A]; 23:
2,6-octadien-1-ol, 3,7-dimethyl-, acetate, (2e)- (geranyl acetate) [A]; 24:
3,7-dimethyl-1,6-octadien-3-ylacetate (linalylacetate) [A]; 25:
allyl-3-cyclohexylpropionate [A]; 26: amylcinnamicaldehyde or
2-phenylmethylene-heptanal [A]; 27: 3-methyl-3-phenylglycidate (aldehyde c16)
[A]; 28: iso-amylsalicilate [A]; 29: benzylbenzoate [A]; 30: diethylphthalate
(dep) [A]; 31: 2-hydroxy-4-methoxybenzophenone (chimassorb 90) [A]; 32:
2-methyl-undecanal (aldehyde c12 mna) [B]; 33: 3,7-dimethyl-6-octen-1-ylacetate
[B]; 34: 3-[4-tert,-buthylphenyl]-2-methylpropanale (lilial) [B]; 35:
2,4-di-t-butylphenol [B]; 36: 2,6-di-t-butylphenol [B]; 37:
4-(2,6,6-trimethyl-2-cyclohexen-1-yle)-3-methyl-3-buten-2-one
(methylionone-gamma) [B]; 38:
5-(2,6,6-trimethyl-2-cyclohexen-1-yle)-3-methyl-3-buten-2-one
(methylionone-alpha) [B]; 39:
4-[4-methyl-4-hydroxyamyl]-3-cyclohexen-carboxaldehyde (lyral) [B]; 40:
2-hexyl-3-phenyl-2-propenal [B]; 41: 2,5-tertbutyl-4-hydroxy-toluene [B]; 42:
2,6-di-tert-butyl-4-methylphenol [B]; 43: 2,6-di-tert-butyl-4-methylphenol
(ionol or BHT) [B]; 44: phenylethylphenylacetate [B]; 45:
nonane-1,3-dioldiacetate (jasmelia) [B]; 46: 2-hydroxy-4-ethanediolbenzophenone
[B]; 47: 2,6-dinitro-1-methyl-3-methoxy-4-tert,-butylbenzole (moschus ambrette)
[B]; 48: 2-hydroxy-4-butoxybenzophenone [B]; 49:
2,4,6-trinitro-1,3-dimethyl-5-tert,-butylbenzene  (moschus xylol) [B]

## 5.3    Effect of the polymer;

The effects of the polymer on the diffusion of migrants are twofold: i) the relaxation of polymer segments affect the renewal of free-volumes around the solute and ii) specific interactions between rigid blocks and the polymer may increase the release time of rigid blocks. In the original free-volume theory of Vrentras and Duda [133](#_ENREF_133), [134](#_ENREF_134), the translation of rigid blocks is assumed to be associated with the local reorganization of special free-volumes, so-called hole free-volume (hFV), whose redistribution is not activated by temperature. As an approximation of the thermal expansion of polymers, the amount of hFV available for the diffusion of migrants is proportional to $T-{{T}_{g}}+{{K}_{\beta }}$ , where ${{T}_{g}}$ is the glass transition temperature and ${{K}_{\beta }}$ a constant possibly dependent on the polymer and the shape of the translating block. This interpretation is illustrated on the diffusion of various migrants in rubber ($T\ge {{T}_{g}}$ ) and glassy ($T<{{T}_{g}}$ ) polymers in [Figure 19](#figure19). The apparent scaling exponent $\alpha \left( T,{{T}_{g}} \right)$ decreases rapidly with $T-{{T}_{g}}$;  for linear solutes, Fang et *al.*[185](#_ENREF_185)  demonstrated that its general expression was:

$$
\alpha \left( T-{{T}_{g}} \right)=1+\frac{{{K}_{\alpha }}}{T-{{T}_{g}}+{{K}_{\beta }}}\;\;  with \;\;T>{{T}_{g}}+{{K}_{\beta }}
$$

with ${{K}_{\alpha }}$ essential a polymer-dependent constant and ${{K}_{\beta }}$ a positive constant for rigid blocks smaller than the cross-section of polymer segments and negative otherwise.

Eq. XXX and Figure 19a demonstrate that mass dependence increase rapidly when $T$ approaches ${{T}_{g}}$ or goes below ${{T}_{g}}$. In elastomers and rubbers (with very low ${{T}_{g}}$) as well as in polymer melts, the mass dependence for linear migrants approaches unity. This theoretical value corresponds to the independent displacements of jumping units. Due to the dependence of ${{K}_{\beta }}$ with the type of rigid blocks, Eq.   describes different mass dependence (see the three series of linear molecules in Figure 17) for different molecules. Although approximative, the technique can be used to extrapolate the results of $\ln {{D}_{P}}$ from one polymer to a new polymer (denoted R) via the correction $\frac{\alpha \left( T,T_{g}^{R} \right)}{\alpha \left( T,{{T}_{g}} \right)}\ln M$. The principles are illustrated in Figure 19e by showing that almost any of the presented polymer (except rubbers) can be used to predict the diffusion coefficients of the same substances in polypropylene only by setting the ${{T}_{g}}$ of the source and destination polymer. The same correction can be used to predict diffusivities in plasticized polymers from their values in non-plasticized ones.



<img src="./src/fig19.png" alt="fig19" style="zoom:67%;" />

<a name="figure19"></a>**Figure 19**. Diffusivities of various substances at 25°C in glassy and rubber polymers: (a) raw values, (b) normalized data to remove polymer effects (standardized to     value of 0°C corresponding to atactic polypropylene PP). Filled symbols correspond to n-alkanes (scaling laws     as dashed lines) and empty symbols to various solutes included gases and plastic additives (scaling laws     as continuous lines). Data from Ref. [81](#_ENREF_81). 

  

## 5.4    Activation of diffusion by temperature

 Increasing temperature affects firstly the structure of the polymer, which, in return, facilitates the translation of the migrants. Temperature activation is consequently higher in polymers with high thermal expansion coefficients (higher at rubber state than at glassy state, higher in plasticized polymer than in non-plasticized one, higher in thermoplastics than in thermosets). Below or near ${{T}_{g}}$, apparent activation energies $Ea=R{{T}^{2}}\frac{\partial \ln {{D}_{p}}}{\partial T}$ are dominated by free-volume effects and are higher than at higher temperatures. The effects of temperature for linear solutes near ${{T}_{g}}$ are shown in [Figure 20](#figure20). If the behavior was strictly Arrhenian, the dependence of $\ln {{D}_{p}}$ would be linear with $1/T$, which is not obviously the case for large diffusants near ${{T}_{g}}$ . Apparent activation energies increase with the size of diffusants as ${{E}_{a}}\propto \ln M$. The interested reader will find a consistent review of free-volume and specific-migrant activation in Ref. [136](#_ENREF_136). In any case, it is not recommended to extrapolate ${{D}_{p}}$ values from low temperatures to high temperatures for polymer recycling, as it would overestimate diffusivities and consequently the rate of decontamination of the polymer.

<img src="./src/fig20.png" alt="fig20" style="zoom:60%;" />

<a name="figure20"></a>**Figure 20**. Arrhenius plot of diffusivities of *n*-alkanes in polyethylene terephthalate. Data from Refs. [188](#_ENREF_188), [189](#_ENREF_189), [194](#_ENREF_194).




# 6	Sorption properties and partition coefficients
-----
**Highlights on the sorption and partition properties**

- Along diffusivities, sorption and partitions coefficients are essential properties to assess the kinetics of migration in food.
- Along with fugacities, partial pressures offer the best representation of the force enabling a substance to move from one phase to the other. Practical approximations are given for all practical cases of mass transfer (from a material to a liquid, from a material to a gas phase, etc.) and possible substances (volatile or not).
- Partition coefficients are defined respectively to activity coefficients.
- Binary and ternary Flory isotherms offer estimations of activity coefficients for all applications of migration modeling: mixture of food simulants, wet or plasticized polymers, homo, and copolymers
- Temperature effects are less important between two condensed phases, but dominant between a condensed and a gas phase.
- Activity coefficients and they activation by temperature can be calculated at atomistic scale for arbitrary solutes (non-charged)

-----

As shown in §4.2.1.1 and §5.1.2, sorption properties and partitioning between materials, food and polymer affect both the kinetics of migration and the distribution of migrants at thermodynamical equilibrium. Previous descriptions were essentially counting the number of substances in a representative elementary volume or a compartment (material, food). The present description will explicit the role of pair interactions between the migrant indexed $i$ and all other constituents ($P$ for the polymer, $F$ for the contact phase, $e$ for ethanol, $w$ for water). Generic phases in equilibrium (in contact or not) are denoted $\alpha $, $\beta $ and $\delta $  regardless they are a solid (polymer or solid food), a liquid contact phase or a gas (headspace, storage environment).

 

## 6.1    Some definitions

### 6.1.1     Chemical potentials, fugacities, and activities

Two phases, $\alpha $ and $\beta $, are said at thermodynamical equilibrium if they are at the same temperature $T$, pressure $P$ and the same chemical potential for any migrant $i$: ${{\mu }_{i,\alpha }}\left( T,P \right)={{\mu }_{i,\beta }}\left( T,P \right)$. The concept of thermodynamical potential was proposed by G. N. Lewis, with initially in mind pure ideal gases verifying ${{\left. \frac{\partial {{\mu }_{i,\alpha }}}{\partial P} \right|}_{T}}={{\bar{V}}_{i}}$ , with ${{\bar{V}}_{i}}=\frac{RT}{P}$ is the molar volume of the gas $i$. After integration at a constant temperature, the change in chemical potential from pressure ${{P}^{ref}}$ to $P$ is $RT\ln \frac{P}{{{P}^{ref}}}$for pure ideal gas. Lewis generalized to phase mixtures (ideal or not) by replacing pressure by a function ${{f}_{i}}$, called fugacity and assessing the capacity of a substance literally to flee:

$$
{{\mu }_{i,\alpha }}-\mu _{i,\alpha }^{ref}=RT\ln \frac{{{f}_{i,\alpha }}}{f_{i,\alpha }^{ref}}
$$

Reference values at the temperature $T$ $\mu _{i,\alpha }^{ref}$ and $f_{i,\alpha }^{ref}$ are physically related, but the choice of the reference chemical potential or the reference fugacity is arbitrary. The ratio ${{a}_{i,\alpha }}=\frac{{{f}_{i,\alpha }}}{f_{i,\alpha }^{ref}}$ defines the activity of the solute $i$ in the phase $\alpha $; it provides namely a measure of the difference of the substance(s chemical potential in $\alpha $ with its reference state.  To describe mass transfer between the phases $\alpha $ and $\beta$ (*i.e.*, with $f_{i,\alpha }^{eq}=f_{\beta ,\alpha }^{eq}$) at equilibrium is convenient to choose the same reference state $\mu _{i,\alpha }^{ref}=\mu _{i,\beta }^{ref}$ in both phases. Following the intuition of Lewis, equivalent partial pressures in an ideal gas offer an almost universal potential to estimate the potential of transfer of any substance (*e.g.*, a volatile organic compound, a liquid plasticizer, a crystalline pigment or a poorly soluble additive or mineral oil residues). The reference fugacity needs to be adapted accordingly as shown in [Table 7](#table7).

<a name="table7"></a>**Table 7**. Expressions of practical partial pressures and saturation concentrations in relationship with the reference state of the substance in the conditions where its migration its studied.

| Application                                                  | Choice for the fugacity: ${{f}_{i}}$                         | Choice for the reference fugacity: $f_{i}^{ref}$             | Relationship with the concentration of substance $i$ in the phase $\alpha $: ${{C}_{i,\alpha   }}$ (SI unit in mol⋅m-3) |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Volatile substance in a gas phase ($\alpha $=gas phase)      | partial pressure: ${{p}_{i}}$                                | total pressure: P                                            | ${{p}_{i}}=RT{{C}_{i,\alpha }}$                              |
| Dissolved substance in   a condensed phase (polymer, liquid) with a liquid reference state ($\alpha $=P   or F) | partial pressure (equivalent partial pressure in the theoretical gas phase $\delta $ in equilibrium   with $\alpha $): ${{p}_{i}}$ | saturation pressure of the pure substance at the same temperature: ${{p}_{i,sat}}\left(   T \right)$ | $\frac{{{p}_{i}}}{{{p}_{i,sat}}}=\gamma _{i,\alpha }^{v}{{\phi }_{i,\alpha   }}=\gamma _{i,\alpha }^{v}{{\bar{V}}_{i}}{{C}_{i,\alpha }}$   with $\gamma _{i,\alpha }^{v}$ the activity coefficient of the substance   $i$ in the phase $\alpha $relatively to   its volume fraction ${{\phi }_{i,\alpha }}$.       At saturation, one gets: $C_{i,\alpha }^{sat}=\frac{1}{\gamma _{i,\alpha   }^{v}{{{\bar{V}}}_{i}}}$ |
| As above but with a solid reference state                    | partial pressure (equivalent partial pressure in the theoretical gas phase $\delta $ in equilibrium   with $\alpha $): ${{p}_{i}}$ | Partial pressure at the surface of the crystal expressed as: ${{\left(   \frac{f_{i,pure}^{S}}{f_{i,pure}^{L}} \right)}_{\left( T   \right)}}{{p}_{i,sat}}\left( T \right)$ with  ${{\left(   \frac{f_{i,pure}^{S}}{f_{i,pure}^{L}} \right)}_{\left( T \right)}}$ being   the ratios of fugacities of the pure solute between solid (pure crystalline)   and molten (pure amorphous) state. | $\frac{{{p}_{i}}}{{{\left( \frac{f_{i,pure}^{S}}{f_{i,pure}^{L}}   \right)}_{\left( T \right)}}{{p}_{i,sat}}}=\gamma _{i,\alpha }^{v}{{\phi   }_{i,\alpha }}=\gamma _{i,\alpha }^{v}{{\bar{V}}_{i}}{{C}_{i,\alpha }}$ |

 

Migrants with a reference solid state in the conditions of migration are very common in plastics and thermosets. They encompass almost all antioxidants, colorants and pigments. They are liquid in conditions of processing of the polymer (*e.g.*, above 160°C), but solid in the conditions of service of the finished material (at $T$and $P$). The partial pressures are lower in this case but non-zero. The variation of fugacities between a pure crystal and its pure amorphous state depends on the temperature difference between $T$ and its melting temperature, ${{T}_{i,m}}$ as shown in [Figure 21](#figure21). 

<img src="./src/fig21.png" alt="fig21" style="zoom: 50%;" />

<a name="figure21"></a>**Figure 21**. Ratio of fugacities between pure solid and amorphous states for 11 model migrants (data from Ref. [31](#_ENREF_31)) and its continuous approximation proposed in Figure S1 of the Supporting Information of Ref. [160](#_ENREF_160).

It is worth noticing that thermodynamical models and simulations calculations may use different reference states (*e.g.*, amorphous reference states even if the reference one is solid). The choices presented here are consistent with the definition of the molar solvation energy:

$$
\begin{align}

  & \Delta G_{i,solvated\ in\ \alpha }^{solvatation}={{\mu }_{i,\alpha }}-\mu _{i,\delta }^{\delta =ideal\ gas}=RT\ln \left( \frac{{{p}_{i,sat}}\left( T \right){{{\bar{V}}}_{i}}}{RT}\frac{\gamma {{_{i,\alpha }^{v}}_{i}}{{\phi }_{i,\alpha }}}{{{C}_{i,\delta }}} \right) \\ 

 & =RT\ln \left( \frac{{{p}_{i,sat}}\left( T \right){{{\bar{V}}}_{i}}}{RT}{{K}_{i,\alpha /\delta }} \right)  

\end{align}
$$

with ${{K}_{i,\alpha /\delta }}$ being the effective partition coefficient between the condensed phase $\alpha $ and the theoretical ideal gas phase $\delta $. Respectively to mass transfer, the use of partition coefficients should be preferred to solvation energies. Free energies need to be preferred to analyze and to interpret quantitatively the nature of the interactions (van-der-Waals, electrostatic, hydrogen bonding) in the phase $\alpha $. Thermodynamic integration at molecular scale offers direct access to $\Delta G_{i,solvated\ in\ \alpha }^{solvatation}$ without requiring a theoretical separation of enthalpic and entropic effects.

### 6.1.2     Effective partition coefficients between P and F

From relationships presented in Table 7, effective partition coefficients between a material P and the contacting phase is given by the ratio of activity coefficients. For semi-crystalline polymers, it is well accepted that crystallites and crystalline phases are impenetrable to migrants. By assuming that no migrating substance has been trapped during processing in the crystalline phase, the effective partition coefficient reads:

$$
{{K}_{i,F/P}}=\frac{1}{\left( 1-c \right)\left( 1-\epsilon  \right)}\frac{{{{\bar{V}}}_{i}}\phi _{i,F}^{eq}}{{{{\bar{V}}}_{i}}\phi _{i,P}^{eq}}=\frac{1}{\left( 1-c \right)\left( 1-\epsilon  \right)}\frac{\gamma _{i,P}^{v}}{\gamma _{i,F}^{v}}
$$

where $\phi _{i,P}^{eq}$ refers to the volume fraction in the amorphous phase of the material (possibly porous and semi-crystalline); $c$and $\epsilon $ are the volume crystallinity and porosity of the material. Thermoplastics are not porous materials and $\epsilon =1$.

For very porous materials, such as papers, it is preferable not to homogenize concentration between all phases (solid crystalline, solid amorphous and gas) and to choose the fibers as P and to adapt the exchange surface area to the shape of fibers. The characteristic dimension of the material should be chosen accordingly commensurable to the half-diameter of fibers.


## 6.2    Sorption isotherms

Sorption isotherms are experimental or theoretical curves, which relate the amounts absorbed by a condensed phase in relationship with an applied activity at constant temperature and pressure. Their use is more general than partition coefficients, when the values of ${{\left\{ \gamma _{i,k}^{v} \right\}}_{k=P,F}}$ depend on concentrations. In this section, we call isotherm a curve relating the activity coefficient or the excess chemical potential (*i.e.*, in excess respectively to the ideal contribution $RT\ln {{\phi }_{i}}$ ) with the composition in the mixture.

### 6.2.1     Linear isotherms

For most applications at infinite dilution, activity coefficients can be assumed independent on composition. The relationship between mass uptake and partial pressure is linear and governed by a Henry isotherm (see Eq. XX). The relationship between the Henry constant ${{\left\{ {{k}_{i,k}} \right\}}_{k=P,F}}$ and ${{\left\{ \gamma _{i,k}^{v} \right\}}_{k=P,F}}$ is given by:

$$
{{k}_{i,k}}={{p}_{i,sat}}\left( T \right){{\bar{V}}_{i}}\gamma _{i,k}^{v}\;\;\;for\;\;\; k=P,F
$$

It is worth noticing that if the condition of infinite dilution is well verified in foods regardless the considered substance (migrations are expected to be low), it may not be verified in the material for substances used close to the saturation (*e.g.*, pigments) and for plasticizing substances (*e.g.* used at weight fractions up to 50%). The proposed description assumes that the substances are well mixed and exclude, by definition, surfactants, and substances causing blooming.

### 6.2.2     Binary Flory isotherms


**General formulation**


Flory-Huggins theory offers a robust framework to account for concentration effects. The theory extends the regular solution theory for liquid mixtures to mixtures with molecules with dissimilar sizes such as solutes mixed with polymers segments or large additives dispersed in a food simulant. Enthalpic and entropic interactions are calculated on a lattice assuming that the mixture is incompressible and that molecules fill space commensurably to their molar volumes. For binary mixtures, $i+P$ or $i+F$ , the activity coefficient is given by:

$$
\ln {{\gamma }_{i,k}}\left( {{\phi }_{i,k}},T \right)=\left( 1-\left( \frac{1}{{{r}_{i,k}}}-n_{k\ around\ i}^{compressible} \right) \right)\left( 1-{{\phi }_{i,k}} \right)+\chi _{i,k}^{\left( T \right)}{{\left( 1-{{\phi }_{i,k}} \right)}^{2}}
$$

where ${{r}_{i,k}}$ is the size of the host molecule respectively to the size of the solute $i$ . In a polymer, the chain is considered infinitely long ${{r}_{i,P}}\to \infty $. In foods and a liquid simulating the food, $r_{i,F}^{-1}$ represents the number of molecules of F displaced by the insertion of the substance $i$ in F. For most of the migrants, $r_{i,F}^{-1}$ is expected to be larger than unity in water and ethanol and lower than unity in oil. ${{r}_{i,F}}$can be approximated as $\frac{{{{\bar{V}}}_{F}}}{{{{\bar{V}}}_{i}}}$ with \[n_{k\ around\ i}^{compressible}\] accounting for the partial compressibility of the molecules of $i$ and$F$. Rigorously, partial molar volumes should be used instead of molar volumes.


**Approximation at infinite dilution**

At infinite dilution and pending estimations of the Flory-Huggins coefficients ${{\left\{ {{\chi }_{i,k}} \right\}}_{k=P,F}}$, Eq.  provides an estimator of activity and partition coefficients in amorphous regions:

$$
\begin{align}

  & \left( 1-c \right)\left( 1-\epsilon  \right)K_{i,F/P}^{\left( T \right)}=\frac{\gamma _{i,P}^{v}}{\gamma _{i,F}^{v}}\approx \frac{\exp \left( 1+\chi _{i,P}^{\left( T \right)} \right)}{\exp \left( 1-{{r}_{i,F}}+\chi _{i,F}^{\left( T \right)}+\text{compressible}\ \text{correction} \right)} \\ 

 & =\exp \left( \chi _{i,P}^{\left( T \right)}-\chi _{i,F}^{\left( T \right)}+{{r}_{i,F}}-\text{compressible}\ \text{correction} \right)  

\end{align}
$$

**Effect of temperature on partitioning**

According to Eq. , the variation of ${{K}_{i,F/P}}$ with temperature depends on the variation of the difference $\chi _{i,P}^{\left( T \right)}-\chi _{i,F}^{\left( T \right)}$  with $T$. When $\chi _{i,P}^{\left( T \right)}$ and $\chi _{i,F}^{\left( T \right)}$ are of the same sign, their variations with $T$ are similar and compensate each other. When  $\chi _{i,P}^{\left( T \right)}\cdot \chi _{i,F}^{\left( T \right)}<0$, there is no compensation and the effect of temperature is maximum (case of polar solute distributed between a polar and an apolar phase). In both cases, the absolute values of ${{\left\{ {{\chi }_{i,k}} \right\}}_{k=P,F}}$ tend to decrease with temperature, so each mixture becomes progressively ideal (*i.e.*, with no enthalpy of mixing). The temperature dependence is predicted via Eqs. -.

### 6.2.3     Ternary Flory isotherms

Eq.  can be generalized to ternary mixtures with two practical applications: i) the estimation of activity coefficients in polar polymers, which contain some amount of water, and ii) the estimation of activity coefficients in water-ethanol mixtures.

The activity coefficient in a wet polymer, associated with a volume fraction of water ${{\phi }_{w}}$ , depends on the three pair of Flory-Huggins coefficients, ${{\chi }_{i,P}}$ , ${{\chi }_{w,P}}$ (water in dry P), ${{\chi }_{i,w}}$ (solute in water), as (see Eq. 9a in Ref. [31](#_ENREF_31)):

$$
\begin{align}

  & \ln {{\gamma }_{i,{{P}_{wet}}}}=\left( 1-{{\phi }_{i}} \right)-{{\phi }_{w}}\frac{{{{\bar{V}}}_{i}}}{{{{\bar{V}}}_{w}}}-\left( 1-{{\phi }_{i}}-{{\phi }_{w}} \right)\frac{1}{{{r}_{i,P}}}+\left( {{\chi }_{w,i}}\frac{{{{\bar{V}}}_{i}}}{{{{\bar{V}}}_{w}}}{{\phi }_{j}}+{{\chi }_{i,P}}\left( 1-{{\phi }_{i}}-{{\phi }_{w}} \right) \right)\left( 1-{{\phi }_{i}} \right) \\ 

 & -{{\chi }_{w,P}}\frac{{{{\bar{V}}}_{i}}}{{{{\bar{V}}}_{w}}}{{\phi }_{w}}\left( 1-{{\phi }_{i}}-{{\phi }_{w}} \right)  

\end{align}
$$

The activity coefficient in water-ethanol mixtures or in any mixture of two miscible liquids F1 and F2 can be estimated similarly [158](#_ENREF_158):

$$
\begin{align}

  & \ln \gamma _{i,{{F}_{1}}+{{F}_{2}}}^{v}=\left( 1-{{\phi }_{i}} \right)-r_{i,{{F}_{1}}}^{-1}{{\phi }_{{{F}_{1}}}}-r_{i,{{F}_{2}}}^{-1}{{\phi }_{{{F}_{2}}}}+\left[ \left( {{\chi }_{i,{{F}_{1}}}}{{\phi }_{{{F}_{1}}}}+{{\chi }_{i,{{F}_{2}}}}{{\phi }_{{{F}_{2}}}} \right)\left( {{\phi }_{{{F}_{1}}}}+{{\phi }_{{{F}_{2}}}} \right) \right] \\ 

 & -{{\chi }_{{{F}_{1}},{{F}_{2}}}}\left( \frac{V_{i}^{vdw}}{V_{{{F}_{2}}}^{P}} \right){{\phi }_{{{F}_{1}}}}{{\phi }_{{{F}_{2}}}}-{{\phi }_{i}}{{\phi }_{{{F}_{2}}}}{{\phi }_{{{F}_{2}}}}\left( \frac{d{{\chi }_{i,{{F}_{2}}}}}{d{{\phi }_{{{F}_{2}}}}} \right)-{{\phi }_{i}}{{\phi }_{{{F}_{1}}}}{{\phi }_{{{F}_{2}}}}\left( \frac{\partial {{\chi }_{i,{{F}_{1}}}}}{\partial {{\phi }_{{{F}_{2}}}}} \right) \\ 

 & -{{\phi }_{i}}\phi _{_{{{F}_{1}}}}^{2}\left( \frac{\partial {{\chi }_{i,{{F}_{1}}}}}{\partial {{\phi }_{{{F}_{2}}}}} \right)-{{\phi }_{i}}\phi _{_{{{F}_{1}}}}^{2}\left( \frac{\partial {{\chi }_{i,{{F}_{1}}}}}{\partial {{\phi }_{{{F}_{1}}}}} \right)-\left( \frac{V_{i}^{vdw}}{V_{{{F}_{2}}}^{P}} \right){{\phi }_{{{F}_{1}}}}\phi _{_{{{F}_{2}}}}^{2}\left( \frac{\partial {{\chi }_{{{F}_{1}},{{F}_{2}}}}}{\partial {{\phi }_{{{F}_{2}}}}} \right) \\ 

 & -\left( \frac{V_{i}^{vdw}}{V_{{{F}_{2}}}^{P}} \right){{\phi }_{{{F}_{1}}}}\phi _{_{{{F}_{2}}}}^{2}\left( \frac{\partial {{\chi }_{{{F}_{1}},{{F}_{2}}}}}{\partial {{\phi }_{{{F}_{1}}}}} \right)-{{\phi }_{i}}{{\phi }_{{{F}_{1}}}}\phi _{_{{{F}_{2}}}}^{2}\left( \frac{\partial {{\chi }_{i,{{F}_{1}},{{F}_{2}}}}}{\partial {{\phi }_{{{F}_{1}}}}} \right)-{{\phi }_{i}}\phi _{_{{{F}_{1}}}}^{2}{{\phi }_{{{F}_{2}}}}\left( \frac{\partial {{\chi }_{i,{{F}_{1}},{{F}_{2}}}}}{\partial {{\phi }_{{{F}_{2}}}}} \right) \\ 

 & -{{\chi }_{i,{{F}_{1}},{{F}_{2}}}}{{\phi }_{{{F}_{1}}}}{{\phi }_{{{F}_{2}}}}\left( 1-2{{\phi }_{i}} \right) \\ 

\end{align}
$$

with ${{\chi }_{i,{{F}_{1}},{{F}_{2}}}}$ being a ternary Flory-Huggins coefficient, whose contribution can be neglected in the absence of specific ternary complex in the solution. The binary Flory-Huggins coefficient between water (F1) and ethanol (F2) is obtained from the molar heat of mixing of the mixture, $\Delta H_{{{F}_{1}}+{{F}_{2}}}^{molar}$ , as ${{\chi }_{{{F}_{1}},{{F}_{2}}}}={{y}_{{{F}_{2}}}}\left( 1-{{\phi }_{2}} \right)\frac{\Delta H_{{{F}_{1}}+{{F}_{2}}}^{molar}}{RT}$. Using the polynomial approximation proposed for water-ethanol mixtures in Ref. [195](#_ENREF_195). [Figure 22](#figure22) reports the value of ${{\chi }_{water,ethanol}}$ with the volume fraction in ethanol as well as the values for common food simulants expressed in alcohol-by-volume (abv) instead of volume fraction. The volume fractions ethanol and water in the mixture are tabulated from their partial molar volumes between 10°C and 60°C in [Table 8](#table8).

<img src="./src/fig22.png" alt="fig22" style="zoom: 67%;" />

<a name="figure22"></a>**Figure 22**. Variation of the binary Flory-Huggins coefficient in water-ethanol mixtures. Data from Ref. [158](#_ENREF_158). “abv” values represent the equivalent alcohol strength (alcohol-by-volume at 20°C and atmospheric pressure). Simulant C: 10% ethanol for alcoholic foodstuffs; simulant D1: 50% ethanol for high alcoholic and milk.

 

<a name="table8"></a>**Table 8**. Density (kg⋅m-3) of water-ethanol mixture and corresponding volume fraction of ethanol ($\phi_{ethanol}$) between 10 and 60°C. Volume fractions are calculated from partial molar volumes. $ abv$=alcohol-by-volume at 20°C.

| $abv$                                                        | 10°C                                                         | 10°C                                                         | 20°C                                                         | 20°C                                                         | 30°C                                                         | 30°C                                                         | 40°C                                                         | 40°C                                                         | 60°C                                                         | 60°C                                                         |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| $abv$                                                        | $\rho_{mixture}$ $\left(kg\cdot m^3\right)$                  | $\phi_{ethanol}$                                             | $\rho_{mixture}$ $\left(kg\cdot m^3\right)$                  | $\phi_{ethanol}$                                             | $\rho_{mixture}$ $\left(kg\cdot m^3\right)$                  | $\phi_{ethanol}$                                             | $\rho_{mixture}$ $\left(kg\cdot m^3\right)$                  | $\phi_{ethanol}$                                             | $\rho_{mixture}$ $\left(kg\cdot m^3\right)$                  | $\phi_{ethanol}$                                             |
| 5.0 %  7.0 %    10%   12 %   15 %   20 %   30 %   40 %   50 %   60 %   70 %   80 %   90 %   95 %   99 % | 992.6087   990.0815    986.5772   984.3969   981.3355   976.6829   967.1218   954.5007   937.5325   917.0055   893.7478   867.6778   837.7669   819.9502   802.7596 | 0.0463   0.0640    0.0903   0.1080   0.1339   0.1773   0.2695   0.3722   0.4773   0.5812   0.6847   0.7889   0.8940   0.9478   0.9900 | 991.0594 988.4460    984.7554   982.3832   978.9701   973.5916   962.2368   948.0405   930.1519   909.1314   885.5598   859.2628   829.2000   811.3845   794.2436 | 0.0465   0.0643    0.0912   0.1092   0.1356   0.1800   0.2739   0.3755   0.4795   0.5828   0.6858   0.7896   0.8946   0.9479   0.9900 | 988.4199   985.7365    981.8636   979.3060   975.5866   969.5825   956.6964   941.2003   922.4957   901.0173   877.1417   850.6228   820.4500   802.6589   785.6217 | 0.0466   0.0646    0.0920   0.1103   0.1372   0.1826   0.2775   0.3789   0.4816   0.5844   0.6869   0.7906   0.8950   0.9481   0.9901 | 984.8805   982.1110    978.0521   975.3416   971.3450   964.7720   950.5944   933.9794   914.5695   892.6443   868.4696   841.7478   811.4700   793.7375   776.8538 | 0.0467   0.0650    0.0928   0.1113   0.1386   0.1847   0.2807   0.3816   0.4837   0.5859   0.6883   0.7915   0.8954   0.9482   0.9901 | 975.3125   972.3164    967.9100   965.0044   960.6287   953.2020   937.1834   918.9400   898.2414   875.3261   850.5234   823.5828   792.9257   774.9002   758.1182 | 0.0473   0.0660    0.0938   0.1125   0.1407   0.1883   0.2857   0.3863   0.4877   0.5894   0.6908   0.7929   0.8971   0.9491   0.9900 |

 

### 6.2.4     Binary Flory-Huggins coefficients in a copolymer AB

Based on calculations at molecular scale, Nguyen *et al.* [*160*](#_ENREF_160) demonstrated that Eq. (12) in Ref. [31](#_ENREF_31) was acceptable for copolymers consisting of repeated blocks shorter than the persistence length of polymer segments. The corresponding Flory-Huggins interaction coefficient in a copolymer AB, ${{\chi }_{i,\left( AB \right)}}$,  reads:

$$
{{\chi }_{i,\left( AB \right)}}={{\chi }_{i,A}}{{\phi }_{A}}+{{\chi }_{i,B}}{{\phi }_{B}}-{{\chi }_{AB}}{{\phi }_{A}}{{\phi }_{B}}
$$

Eq. XXX is based on an averaging of all possible contacts between $i-A$ and $i-B$ , where the term ${{\chi }_{AB}}{{\phi }_{A}}{{\phi }_{B}}$ represents the additional cohesion energy brought by the interactions between A and B. It can be generalized to more complex copolymers while the contacts can be assumed perfectly random (**i.e.**, no macro or microphase separation).

 

## 6.3    High throughput calculations of Flory-Huggins coefficients at atomistic scale

### 6.3.1     Justification and limitations

Molecular modeling offers a good alternative to time-consuming and complex experiments to estimate Flory-Huggins coefficients for various substances: monomers, oligomers, solvent, additives, residues, breakdown products, non-intentionally added substances. The results can be tabulated in advance and used directly with Eqs. - for a broad range of applications. With this respect, they are more intrinsic than partition coefficients. In details, molecular modeling can be seen as an alternative to earlier group contribution methods relying on estimating Flory-Huggins coefficients from solubility parameters (see Ref. [110](#_ENREF_110), [149](#_ENREF_149)). The limits of the approach have been discussed in Ref. [112](#_ENREF_112) and compared with calculations at atomistic scale. In shorts, group contribution methods provide only an average picture of the interactions.  The real conformation of molecules and the distance between interacting chemical groups are, in particular, not preserved. The principles of calculation of $\chi _{i,P}^{\left( T \right)}$ and $\chi _{i,F}^{\left( T \right)}$ have been reviewed in Refs.[112](#_ENREF_112), [158-161](#_ENREF_158). They apply to any homo and copolymer as well as almost to any solute while no net charge is present on the polymer, on the food simulant or the solute. The main limitations are intrinsic to the Flory approximation itself: no energetic barrier should exist in the host system (P and F) so that all states are accessible. It is not true in glassy polymers, where hysteresis is frequent. The contribution of the subcooling and accumulated elastic energy can be *introduced a posteriori* by combining Flory and free-volume theories together[148](#_ENREF_148), [196](#_ENREF_196). These extensions are, nevertheless, beyond the scope of the chapter.

### 6.3.2     Principles

The Flory-Higgins coefficients, ${{\left\{ {{\chi }_{i,k}} \right\}}_{k=P,F}}$ is defined as the dimensionless mixing energy (enthalpy) in excess relatively to pure compounds:

$$
\chi _{i,k}^{\left( {{n}_{k}},T \right)}=\frac{{{\left\langle h_{i+k}^{{{n}_{k}}} \right\rangle }_{T}}+{{\left\langle h_{k+i}^{{{n}_{k}}} \right\rangle }_{T}}-{{\left\langle h_{k+k}^{{{n}_{k}}} \right\rangle }_{T}}-{{\left\langle {{h}_{i+i}} \right\rangle }_{T}}}{2RT}\;\;\;for\;\;\;k=P,F
$$

where ${{\left\langle X \right\rangle }_{T}}$ represents an ensemble-average of $X$, ${{n}_{P}}$ is the length used in the approximation of the polymer and ${{n}_{F}}=1$. In agreement with the original Flory approximation, enthalpies ${{\left\langle {{h}_{A+B}} \right\rangle }_{T}}$ are estimated by summing energies of contact${{\epsilon }_{AB}}$when B is contacting the seed molecule A. ${{\left\langle {{h}_{B+A}} \right\rangle }_{T}}$represents the same energy when B is used as seed molecule. In practice, ${{\epsilon }_{AB}}$ is calculated by choosing an orientation randomly for the contact molecule and by translating it along a random line until at least one point of contact is established between the van-der-Walls surfaces of the contact and seed molecules. The process is repeated for all conformers and stereoisomers considered. Finally, ${{h}_{A+B}}$is estimated as the product of contact energies and the number of neighbors ${{z}_{AB}}$ (number of B molecules surrounding A):

$$
{{\left\langle {{h}_{A+B}} \right\rangle }_{T}}={{\left\langle {{n}_{cooperative}}{{z}_{AB}}{{\varepsilon }_{AB}} \right\rangle }_{T}}\approx {{n}_{cooperative}}\left\langle {{z}_{AB}} \right\rangle {{\left\langle {{\varepsilon }_{AB}} \right\rangle }_{T}}
$$

Eq. XXX assumes that ${{\epsilon }_{AB}}$and ${{z}_{AB}}$ statically independent (zero covariance). For polymers, the property of independence is achieved from a sufficiently large ${{n}_{P}}$ value so that the surface of contact of the polymer is independent of the length of the considered polymer. The main advantage of the whole approach is that there is no need to represent entanglements in the polymer and free-volume. The shape of the backbone of an infinitely long chain with shorter oligomers by preventing head and tail atoms from coming in contact with any van-der-Wall surface. Cooperative hydrogen bonding is accounted by using a value of ${{n}_{cooperative}}$ greater than one. The latent heat of vaporization of water can be correctly approximated by ${{\left\langle {{h}_{water+water}} \right\rangle }_{T}}$ using a value of ${{n}_{cooperative}}$ different to unity. This value depends on the type of forcefield used to simulate water. As an example, the rigid water model governed by the TIP4P forcefield gives an acceptable value with ${{n}_{cooperative}}=1$ whereas ${{n}_{cooperative}}=4$ is required with same forcefield but using to three-point charges (forcefield TIP3P). The number 4 enforces in this case that any water molecule is on average involved in 4 hydrogen bonds of similar strength.

Contact energies are calculated irrespectively to any temperature consideration. The effect of temperature is recovered by weighting the distribution of contact energies with the Boltzmann factor $B\left( \epsilon  \right)$ :

$$
{{\left\langle {{\varepsilon }_{AB}} \right\rangle }_{T}}=\tfrac{\int\limits_{-\infty }^{+\infty }{pr\left( \varepsilon  \right)B\left( \varepsilon  \right)\varepsilon d\varepsilon }}{\int\limits_{-\infty }^{+\infty }{pr\left( \varepsilon  \right)B\left( \varepsilon  \right)d\varepsilon }}\;\;\;with\;\;\;B\left( \varepsilon  \right)=\exp \left( -{\varepsilon }/{\left( RT \right)}\; \right)
$$

At the price of calculating two integrals, Eq.  can be used to estimate $\chi _{i,k}^{\left( {{n}_{k}},T \right)}$ at several temperatures. Conformers need to be generated in a way they are representative of their conformations (radial distribution, constraints of torsion) in the corresponding condensed phase. In practice, they are sampled from molecular dynamics simulations of the equivalent condensed phase.

# 7     Probabilistic modeling of the migration

------
**Highlights on probabilistic modeling of the migration**

- Conventional modeling calculates point estimates, associated to the most likely combination of inputs. Probabilistic modeling calculates the statistic distribution of the amount transferred for all combination of the parameters (likely or not).
- Theoretical results make the calculations of sensitivities and statistic distributions almost as fast as point estimates.
- The combined effects of uncertainty and variability (*e.g.* residence times at specific temperatures) can be analyzed together and used to provide conservative estimates with controlled a risk of underestimation.
- The same approach can identify the most influential parameters acting on the value of the migration or the final decision value.

------


## 7.1      Beyond intuition

Any risk assessment procedure needs to account the possible variabilities in the considered scenario (*e.g.* variable temperature, contact time) and the numerous sources of uncertainties inherent to the limitations of our knowledge and oversimplifications. Variability and uncertainty can be easily recognized and separated by noticing that only uncertainty can be reduced by bringing additional knowledge or refinement. By contrast, variability represents multiple values of several instances (lots, compositions, final use), storage conditions, etc. For compliance testing, conservative assumptions are mandatory, but the relationship between the maximization of parameters (or their minimization depending on the case) and the maximization of the amount transferred is straightforward only in simple configurations: one material or one single layer, one step, no variable conditions. 

The intuitive approach for is illustrated in Figure 8 for a single component and monolayer packaging in contact with food. When the whole food-packaging system is perfectly impervious (no loss to the outside), the cumulated amount leaving the packaging-food interface is a monotonic function of the time, the initial concentration, the diffusion coefficient in the polymer, the chemical affinity for food, the temperature, etc. As a result, choosing a conservative or upper bound for all inputs guarantees an overestimation of the food contamination. In the presence of multiple materials or steps, the property of monotonicity between parameters and inputs is not mathematically verified anymore. In particular, food contamination can be maximum before reaching equilibrium. For example, overestimating all diffusion coefficients or partition coefficients in laminates will spread migrants everywhere instead of bringing them faster to the contacting phase. For laminates, methods described in §4.2.3.5  and theorized in [98](#_ENREF_98)  can be used to recover conservative estimates. The calculation procedure consists in splitting the contribution of $n$ components/materials into $m\le n$independent simulations and in cumulating the concentration in food: ${{\left\{ C_{i,F}^{k} \right\}}_{k=1..m}}$ . The algorithm is sketched for a substance $i$ present simultaneously in the walls (bilayer AB, B is contact) and the cap of the bottle (C). The problem comprises $n=3$ materials and requires $n=2$ simulations:

- One-dimension mass transfer simulation from B to the food (without A and C) to calculate $C_{i,F}^{B}$ 
- One-dimension mass transfer simulation from C to food (without A and B, **i.e.** no walls).

It is worthwhile noticing that the scenario described is also conservative if we assume that the substance has been distributed between A and B before being put in contact with the food. It also covers the case when $i$ is initially located in A and not in B (see Figure 10), but with a higher safety margin.

## 7.2     Epistemic uncertainty

In system engineering, reliability and safety are quantified with respect to some safety margins, defined as the differences between reference values accepted by the regulating body and calculated values. A system is considered safe when the differences calculated for a set of postulated scenarios verifies a minimum distance or when the probability of the distance to be zero or negative is lower than some prescribed value. Introducing conservatism randomly by mixing worst-case bounds may propagate uncertainty fuzzily and leads to uncontrolled overestimation of the amount transferred to the food. At the beginning of the supply chain, the chemical industry and compounders face mainly variability on the different applications of their chemicals and raw materials. On the opposite side of the supply chain, the packaging filler and the retailers face a more different situation with strong uncertainties on the nature of the materials, their thicknesses, their composition. In 2009 and despite the possibilities offered by EU directive 2002/72/EC, we evaluated that migration modeling could be helpful to demonstrate the compliance in finished products only in less than 5% of cases [112](#_ENREF_112). The chief reason was the loss of compositional information along the supply chain. Calculations could be done on part supplies and compounds to produce certificates of food contact compliance, but not on the full system assembled in the intended conditions of use of the packaging. Compositional information is currently better shared in EU and new deformulation techniques provide grounds for spreading calculation practices from the chemical industry to the food industry [114](#_ENREF_114). 

<img src="./src/fig23.png" alt="fig23" style="zoom: 67%;" />

<a name="figure23"></a>**Figure 23**. Illustration of safety margins (SM), overestimations factors (Q) and uncertainty according to the method of calculation: real, likely and very conservative.

 

The best practice for industry relies on sharing minimum information so that a significant safety margin remains for the end-user and keep the utility of migration modeling. As an illustration, [Figure 23](#figure23) compares the safety margins when over-conservative estimates and realistically overestimates are used. We attract the attention to the reader that the location of the “real” value is not known usually so that the overestimation factor cannot be guessed *a priori*. Only the safety margin is directly accessible to calculations. The definition of the safety margin and its usage in various technical guides and supporting risk assessment documents can be inconsistent and confusing. In particular, the concept of safety margin is frequently confused with the concept of overestimation. We used the definition of safety margin applied in medicine to evaluate drugs, in structural engineering, in nuclear engineering, The concept of overestimation applied to some factors including diffusion (${{D}_{P}}$ ) and partition (${{K}_{F/P}}$) coefficients can be misleading. Indeed, their overestimation by a factor ${{Q}_{X}}$ causes an overestimation of the concentration in the medium in contact which is not proportional. For the layer in contact, it varies from$\sqrt{{{Q}_{{{D}_{P}}}}}$ (short contact times) down to 0 (equilibrium) for diffusion ($X={{D}_{P}})$. For $X={{K}_{F/P}}$, it increases from 0 to a value, which depends on the volume of the food.

By evaluating uncertainty associated to realistic estimates, probabilistic modeling offers a robust methodology to assess the effects of the combined sources of uncertainties and finally to have no safety margin at all. The example depicted in shows that the upper limit of likely overestimates including uncertainty (95th percentile) offers a higher safety margin than the very conservative overestimate. The distinction between overestimations and realistic conservative estimates can be exemplified by considering a long contact at a variable temperature (*e.g.* due to transportation). A conservative estimate will calculate migration at the highest temperature, whereas a realistic conservative value will be provided by replacing time with its integral dimensionless version, Fourier number, $Fo=\int\limits_{{{t}_{1}}}^{{{t}_{2}}}{\tfrac{D\left( T\left( t \right) \right)dt}{{{l}^{2}}}}$.  The proposed approximation does not introduce any significant approximation and can be carried out at the same cost as standard simulations (see §4.2.4 and Eq. XXX).



## 7.3      Sensitivity analysis of migration models

### 7.3.1     Local sensitivity analysis

Deterministic modeling and simulation yield the same output (concentration in food, concentration field in the packaging material) for the same set of inputs. The analysis of the sensitivity to input parameters entails evaluating the effect of a modification of each parameter ${{\left\{ {{p}_{k}} \right\}}_{\ddot{\ }k=1..N}}$  (initial concentration, diffusion coefficient, partition coefficient…) on the safety margin, $SM$. Since the parameters have different units, it is convenient to calculate the derivatives of $SM$with respect to the logarithm of each parameter:

$$
\begin{align}

  & \mathbf{J}\left( \mathbf{p} \right)=\mathbf{J}\left( {{\left[ {{p}_{1}}\ldots {{p}_{N}} \right]}^{T}} \right)=\left[ {{\left. \frac{\partial SM}{\partial {{p}_{1}}} \right|}_{{{p}_{2}}..{{p}_{N}}}}{{p}_{1}}\ldots {{\left. \frac{\partial SM}{\partial {{p}_{N}}} \right|}_{{{p}_{1}}..{{p}_{N-11}}}}{{p}_{N}} \right] \\ 

 & =\left[ {{\left. \frac{\partial SM}{\partial \ln {{p}_{1}}} \right|}_{{{p}_{2}}..{{p}_{N}}}}\ldots {{\left. \frac{\partial SM}{\partial \ln {{p}_{N}}} \right|}_{{{p}_{1}}..{{p}_{N-11}}}} \right]  

\end{align}
$$

$\mathbf{J}\left( \mathbf{p} \right)$ is the Jacobian matrix of the migration model and it can be used to evaluate a linear approximation of the safety margin when inputs are changed from ${{\mathbf{p}}_{\mathbf{0}}}$ to $\mathbf{p}$ :

$$
SM\left( \mathbf{p} \right)\approx SM\left( {{\mathbf{p}}_{0}} \right)+\mathbf{J}\left( {{\mathbf{p}}_{0}} \right)\left( \ln \mathbf{p}-\ln {{\mathbf{p}}_{\mathbf{0}}} \right)
$$

Eq. XXX  provides analytical solutions only in simple configurations, but its application is very general. Additionally, implementing the difference $\ln {{p}_{k}}-\ln {{p}_{k,0}}$ as $\ln \frac{{{p}_{k}}}{{{p}_{k,0}}}$prevents the problem of step size known as subtractive cancellation. Its usage is shown for a variant of the realistic but conservative model for monolayer materials presented in 4.2.2 (see Eq. XXX):

$$
{{\left\lceil {{C}_{F}} \right\rceil }_{\left( Fo,{{K}_{F/P}},C_{p}^{0} \right)}}\approx \frac{C_{p}^{0}}{\frac{1}{{{K}_{F/P}}}+\frac{1}{{{L}_{P/F}}}}{{\left\lceil {{{\bar{v}}}^{*}} \right\rceil }_{\left( Fo \right)}}\le \frac{2}{\sqrt{\pi }}\frac{C_{p}^{0}}{\frac{1}{{{K}_{F/P}}}+\frac{1}{{{L}_{P/F}}}}\min \left( \frac{\sqrt{\pi }}{2},\sqrt{Fo} \right)
$$

The evolution of the safety margin $S{{M}_{\left( Fo,{{K}_{F/P}},C_{p}^{0} \right)}}\approx SML-{{\left\lceil {{C}_{F}} \right\rceil }_{\left( Fo,{{K}_{F/P}},C_{p}^{0} \right)}}$ reads when the input values are changed to $\left[ Fo=\frac{{{D}_{P}}t}{l_{P}^{2}},{{K}_{F/P}},C_{p}^{0} \right]$ from an initial set $\left[ \overline{Fo}=\frac{\overline{{{D}_{P}}}t}{l_{P}^{2}},\overline{{{K}_{F/P}}},\overline{C_{p}^{0}} \right]$:

$$
\begin{align}

  & \frac{S{{M}_{\left( Fo,{{K}_{F/P}},C_{p}^{0} \right)}}-S{{M}_{\left( \overline{Fo},\overline{{{K}_{F/P}}},\overline{C_{p}^{0}} \right)}}}{\overline{\left\lceil {{C}_{F}} \right\rceil }}\approx 1-\frac{\left\lceil {{C}_{F}} \right\rceil }{\overline{\left\lceil {{C}_{F}} \right\rceil }} \\ 

 & =-\left[ \frac{1}{\overline{{{K}_{F/P}}}}\frac{1}{\frac{1}{\overline{{{K}_{F/P}}}}+\frac{1}{{{L}_{P/F}}}}\ln \left( \frac{{{K}_{F/P}}}{\overline{{{K}_{F/P}}}} \right)+\ln \left( \frac{C_{P}^{0}}{\overline{C_{P}^{0}}} \right)+\frac{1}{2}\ln \left( \frac{\min \left( \tfrac{\sqrt{\pi }}{2},\sqrt{\text{Fo}} \right)}{\min \left( \tfrac{\sqrt{\pi }}{2},\sqrt{\overline{\text{Fo}}} \right)} \right) \right] \\ 

\end{align}
$$

with $Fo$ $\overline{\left\lceil {{C}_{F}} \right\rceil }=\frac{2\overline{C_{P}^{0}}}{\sqrt{\pi }}\frac{\min \left( \tfrac{\sqrt{\pi }}{2},\sqrt{\overline{\text{Fo}}} \right)}{\left( \frac{1}{\overline{{{K}_{F/P}}}}+\frac{1}{{{L}_{P/F}}} \right)}$ .

Besides showing the interactions and additivity of the different sources of uncertainty, it offers a rapid methodology to identify the main influencing parameters without requiring any simulation or software. Eq. XXX  also applies to numerical simulations, but it requires $N+1$ simulations (the reference one and $N$ variations). When the number of inputs increases, it is preferable to reduce the computational effort by using similitude principles and dimensionless numbers. In the detailed example, the use of a dimensionless time $Fo=\tfrac{{{D}_{p}}t}{l_{p}^{2}}$ enables to test simultaneously the effects of $t$, ${{D}_{P}}$ and ${{l}_{p}}$ . As Eq.  is a first order approximation and, despite the use of a logarithm scale, it is losing accuracy as soon as parameters are tested beyond several factors.

### 7.3.2     Global sensitivity analysis via stochastic simulation

When the number of variables becomes large as well as the intervals to be explored, the statistical sampling of inputs is preferable. Statistical analysis of the outputs can be used to extract the influence of each variable and the probability to have the prescribed threshold exceeded. Each component of the vector $\mathbf{p}$needs to be sampled randomly and uniformly over its interval of interest. The technique is so-called Monte-Carlo trials and its numerical implementation stochastic simulation. By denoting ${{\left\{ {{{\mathbf{\tilde{p}}}}_{k,i}} \right\}}_{i=1..M}}$ $M$ samples chosen around the likely vector $\mathbf{\bar{p}}$ so that only the parameter ${{\left\{ {{p}_{k}} \right\}}_{k=1..N}}$  is modified at a time, and by denoting ${{\left\{ SM\left( {{{\mathbf{\tilde{p}}}}_{k,i}} \right) \right\}}_{\ddot{i}=1..M}}$ the corresponding safety margins. The sample covariance is given by:

$$
{{\mathbf{V}}_{\mathbf{SM}}}=\frac{1}{M-1}\sum\limits_{i=1}^{M}{\left( SM\left( {{{\mathbf{\tilde{p}}}}_{k,i}} \right)-\left\langle SM\left( {{{\mathbf{\tilde{p}}}}_{k,i}} \right) \right\rangle  \right){{\left( SM\left( {{{\mathbf{\tilde{p}}}}_{k,i}} \right)-\left\langle SM\left( {{{\mathbf{\tilde{p}}}}_{k,i}} \right) \right\rangle  \right)}^{T}}}
$$

where $\left\langle SM\left( {{{\mathbf{\tilde{p}}}}_{k,i}} \right) \right\rangle =\frac{1}{M}\sum\limits_{i=1}^{M}{SM\left( {{{\mathbf{\tilde{p}}}}_{k,i}} \right)}$ is the average safety margin, which does not coincide with $SM\left( {\mathbf{\bar{p}}} \right)$ in the general case. $SM\left( {\mathbf{\bar{p}}} \right)$ represents the prediction associated with the 50th percentile.

Eq. XXX generalizes the local sensitivity analysis performed in Eq. XXX, based on small variations and partial derivatives. The concept of covariance enables to screen the whole input spectrum, to identify the interaction and dependency structure on all parameter including the analysis. If the geometry, the temperature and contact time are introduced, the design and the conditions of use can also be explored.

A probabilistic interpretation is achievable; but, as means and covariances provide only the first and second moments, a likely distribution of the safety-margin or the concentration in food is required. If the concentration in food is normally distributed, the problem is fully determined with the first and second moments. This assumption is valid only in the presence of a low range of variabilities and uncertainties. Indeed, the Gaussian distribution is unbounded, and it implies, even with very low probabilities, that the concentration in food could also be negative and the amount transferred could be higher than the amount in the material. The next section removes these limitations for risk assessment and the evaluation of consumer exposure.


## 7.4      Principles of the probabilistic interpretation of mass transfer

Global sensitivity analysis presented in 7.3.2 introduces the first interpretation of mass transfer with a marginal distribution on each input variable, which is assumed to be uniform. The combination of these variables and its interpretation is known as a copula in probability theory and statistics. Copulas describe well the dependence between inputs on the output(s) of a model, but they fail to describe the joint distribution of contamination in realistic situations. The diffusion and partition coefficients, as well as the initial concentration in food are not distributed uniformly. The industry is not applying randomly any concentration value or the molecules do not have random properties. It is our limited knowledge and the variability of practices that spread the inputs around a likely value. In shorts, the sensitivity analysis is a perfect tool to optimize the geometry, the formulation, etc. but it is not appropriated to get a reliable estimate of the probability to have a concentration threshold exceeded.

The principles of probabilistic modeling of migration have been described in [197](#_ENREF_197)  and applied to various cases [96](#_ENREF_96), [97](#_ENREF_97). The central idea is to combine a dimensionless formulation (with a reduced number of input variables) along with random numbers. For the same reason as invoked for the local decomposition in §7.1.3.1, each quantity $X$ is written as the product of a scaling value (with units), denoted $\overline{X}$ (usually the likeliest value), and a random dimensionless number ${{X}^{*}}$distributed around unity. As an illustration, the dimensionless time reads:

$$
Fo=\frac{\left( \bar{D}{{D}^{*}} \right)\left( \bar{t}{{t}^{*}} \right)}{{{\left( \bar{l}{{l}^{*}} \right)}^{2}}}=\frac{\bar{D}\bar{t}}{{{{\bar{l}}}^{2}}}\frac{{{D}^{*}}{{t}^{*}}}{{{l}^{*}}^{2}}=\overline{Fo}F{{o}^{*}}
$$

### 7.4.1     Input distributions

The distributions of ${{X}^{*}}$can be chosen either from experimental measurements or from prior guesses or beliefs. In this second alternative, distributions which have a shape factor and non-negative, denoted ${{f}_{X}}\left( 1,{{s}_{X}} \right)$ , should be preferred (beta, Erlang, exponentially modified Gaussian, exponential, gamma, inverse-gamma, inverse-Gaussian, lognormal, Weibull…). For some quantities, such as concentrations or thicknesses, symmetric distributions are more realistic; truncated normal distributions can be used for this purpose. A non-exhaustive list of practical distributions is given in [Table 9](#table9).

<a name="table9"></a>**Table 9**. Recommended distributions for probabilistic modeling of the migration from monolayer materials. The distributions of ${{\bar{v}}^{*}}$ and $\sqrt{F{{o}^{*}}}$ are *posterior distributions*.

| **random   contribution** | **distribution**            | **recommendations**                                          |                                                              |                   |
| ------------------------- | --------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ----------------- |
| rubber polymers           | glassy polymers             |                                                              |                                                              |                   |
| diffusion coefficient     | ${{\log }_{10}}D_{P}^{*}$   | $Norm\left( O,{{s}_{D}} \right)$                             | ${{s}_{D}}=0.1$                                              | ${{s}_{D}}=0.5$   |
| contact time              | ${{t}^{*}}$                 | $Weib\left( O,{{s}_{t}} \right)$†                            | ${{s}_{t}}~$to be determined                                 |                   |
| initial concentration     | $C{{_{P}^{0}}^{*}}$         | $Norm\left( 1,{{s}_{{{C}_{0}}}} \right)$   truncated†        | $1\le {{s}_{{{C}_{0}}}}\le 5$                                |                   |
| thickness                 | ${{l}^{*}}$                 | $Norm\left( 1,{{s}_{l}} \right)$                             | ${{s}_{l}}~$to be determined                                 |                   |
| mass Biot number          | ${{\log }_{10}}B{{i}^{*}}$  | $Norm\left( O,{{s}_{Bi}} \right)$                            | ${{s}_{Bi}}\to 0$                                            | ${{s}_{Bi}}\to 0$ |
| partition coefficient     | ${{\log }_{10}}K_{F/P}^{*}$ | $Norm\left( O,{{s}_{K}} \right)$                             | $<0.2$                                                       | $<0.2$            |
| Fourier number            | $F{{o}^{*}}^{\frac{1}{2}}$ | $Gamma\left( {{a}_{\text{   }\!\!\Gamma\!\!\text{ }}},{{b}_{\text{ }\!\!\Gamma\!\!\text{ }}} \right)$ | ${{a}_{\text{ }\!\!\Gamma\!\!\text{ }}},{{b}_{\text{   }\!\!\Gamma\!\!\text{ }}}~$: to be calculated |                   |
| concentration in food     | ${{\bar{v}}^{*}}$           | $Beta\left( {{a}_{\text{ }\!\!\beta\!\!\text{   }}},{{b}_{\text{ }\!\!\beta\!\!\text{ }}} \right)$ | ${{a}_{\text{ }\!\!\beta\!\!\text{ }}},{{b}_{\text{ }\!\!\beta\!\!\text{   }}}$ : to be calculated |                   |

† *to be normalized to get a unitary expectation.*

 

### 7.4.2     Estimation of probabilities via Monte-Carlo sampling

Probabilistic modeling aims at determining the cumulated density function (cdf), which can written for monolayer materials as the probability to get a value of ${{C}_{F}}$ lower than arbitrary number $x$ :

$$
pr\left( {{C}_{F}}\le x \right)=F\left( \overline{Fo},\overline{Bi},\overline{{{K}_{F/P}}},\overline{C_{P}^{0}},{{a}_{\Gamma }},{{b}_{\Gamma }},{{s}_{K,}}{{s}_{Bi}},{{s}_{C}} \right)
$$

In the sake of efficiency, $Fo$ (see Eq. XXX) was used instead of $\left( \overline{t},\overline{{{D}_{P}}},\overline{{{l}_{P}}},{{s}_{t}},{{s}_{D}},{{s}_{l}} \right)$. Indeed using Eq. , it can be shown that the distribution *a posteriori* of $\sqrt{F{{o}^{*}}}$ converges in law to a Gamma distribution with parameters $\left( {{a}_{\Gamma }},{{b}_{\Gamma }} \right)$. Ref [114](#_ENREF_114) tabulates their values with ${{s}_{t}}$ , ${{s}_{D}}$ and ${{s}_{l}}$.

The cumulated probability $pr\left( {{C}_{F}}\le x \right)$ (or its complement $pr\left( {{C}_{F}}>x \right)$ ) can be estimated by repeating the simulations for different values of input parameters and by counting the number of occurrences, for which the inequality ${{C}_{F}}\le x$ (or ${{C}_{F}}>x$) is verified. If the intent is to demonstrate that $pr\left( {{C}_{F}}>x \right)$ is low for a sufficiently large $x$, we could think it suffices to apply some worst-case scenarios (bounds of intervals) and to demonstrate that the value $x$ is never exceeded. This approach is correct only if $pr\left( {{C}_{F}}>x \right)=0$, that is for a value of $x$ larger than the one corresponding to a total extraction. In the general case, the intervals of all parameters need to be sampled with each value chosen according to its theoretical prescribed distribution. This technique of picking randomly input values and to launch the corresponding simulation is known as Monte-Carlo simulation.

In practice, for each input quantity $X$, a random number $g$is chosen uniformly between 0 and 1 (function rand() in many programming languages). The specific value to be included in the considered simulation vector will be $F_{X}^{-1}\left( g \right)$, with  ${{F}_{X}}\left( x \right)=pr\left( X\le x \right)$ being the cdf of the variable X. Depending of the size of the intervals, 103 to 105 simulations are required. For multilayer materials, the sampling effort can even higher. The total cost of simulations can be reduced dramatically by tabulating the results in advance for a significant range of dimensionless numbers and by subsequently interpolating the values of interest. These concepts are justified mathematically hereafter. They make the cost of probabilistic modeling close to the cost of deterministic modeling.

  

### 7.4.3     Estimation of joint probabilities via the composition theorem

*The calculation of probability density functions (pdf) associated with a combination of variables ${{X}_{1}},\cdots ,{{X}_{n}}$ (*e.g.*, Eq.* XX)X or with the resolution of partial differential equations is particular computationally expensive and requires a specific treatment of the joint density ${{f}_{{{X}_{1}},\cdots ,{{X}_{n}}}}\left( {{x}_{1}},\cdots ,{{x}_{n}} \right)=pr\left( {{X}_{1}}={{x}_{1}},\cdots ,{{X}_{n}}={{x}_{n}} \right)$. The composition theorem offers a very efficient computational approach for invertible and differentiable transformations:
$$
\begin{align}

  & {{Y}_{i}}={{h}_{i}}\left( {{X}_{1}},\ldots ,{{X}_{n}} \right);\,\,i=1,\ldots ,n \\ 

 & {{X}_{i}}=h_{i}^{-1}\left( {{Y}_{1}},\ldots ,{{Y}_{n}} \right);\,\,i=1,\ldots ,n \\ 

\end{align}
$$

Mathematical functions ${{\left\{ {{h}_{i}} \right\}}_{i=1..n}}$ represent either a variable transformation (*e.g.*, Eq. ) or the full mass transfer model. The joint density of ${{Y}_{1}},\cdots ,{{Y}_{n}}$ is given by the determinant of the Jacobian matrix ${{\left. \partial {{h}_{i}}/\partial {{x}_{i}} \right|}_{{{x}_{j\ne i}}}}$, denoted ${{J}_{g}}$:

$$
\begin{align}

  & {{}_{{{Y}_{1}},\cdots ,{{Y}_{n}}}}\left( {{y}_{1}},\cdots ,{{y}_{n}} \right)= \\ 

 & {{}_{{{X}_{1}},\cdots ,{{X}_{n}}}}\left( h_{1}^{-1}\left( {{y}_{1}},\cdots ,{{y}_{n}} \right),\cdots ,h_{n}^{-1}\left( {{y}_{1}},\cdots ,{{y}_{n}} \right) \right)\cdot {{\left| {{J}_{g}}\left( h_{1}^{-1}\left( {{y}_{1}},\cdots ,{{y}_{n}} \right),\cdots ,h_{n}^{-1}\left( {{y}_{1}},\cdots ,{{y}_{n}} \right) \right) \right|}^{-1}} \\ 

\end{align}
$$

Eq. XXX can be generalized to non-monotonic functions by splitting the transformation into intervals which are locally monotone. For example, in the special case where the dimensionless concentration in food ${{\bar{v}}^{*}}$is not a continuous function of $Fo$ (some multilayer or multicomponent configurations), the pdf of ${{\bar{v}}^{*}}=h\left( Fo \right)$is obtained by the accumulation of the $p$ transformations over $p$ contiguous intervals: 

$$
{{f}_{{{{\bar{v}}}^{*}}}}\left( v \right)=\sum\limits_{k=1}^{p}{{{f}_{Fo}}\left( {{\left. {{{\bar{v}}}^{*}}^{-1} \right|}_{Fo\in {{Y}_{k}}}}\left( v \right) \right){{\left| \frac{d}{dv}{{\left. {{{\bar{v}}}^{*}}^{-1} \right|}_{Fo\in {{Y}_{k}}}}\left( v \right) \right|}^{-1}}}
$$

with ${{\left\{ {{Y}_{k}} \right\}}_{k=1..p}}$ the partitions of $Fo$ where ${{\bar{v}}^{*}}=h\left( Fo \right)$is piecewise monotonic.

### 7.4.2     Some illustrations

Probabilistic modeling must be envisioned as the generalization of deterministic modeling, but its clear definition depends how the normalization of ${{X}^{*}}$ is performed, how deterministic inputs are set (likeliest values or averaged ones) and how the transformation stretches or contracts the probability space. Deterministic modeling provides indeed the unique solution of the initial-value problem. Probabilistic modeling generates the distributions and the confidence tube including the solutions when the initial inputs are changed. The use of dimensionless number, such as a cumulated Fourier number including temperature variations $\int\limits_{0}^{t}{\frac{{{D}_{p}}\left( T\left( t \right) \right)}{l_{p}^{2}}dt}$ or the ones associated to homologous solutes $\frac{{{\left\{ {{D}_{i,p}} \right\}}_{\begin{smallmatrix} i\in \text{homologuous} \text{solutes} \end{smallmatrix}}}t}{l_{p}^{2}}$, offers even a broader interpretation by including random variations of temperatures or uncertainty on diffusion coefficients (see [138](#_ENREF_138) for discussion). In the context of risk assessment, the scaling quantities $\overline{Fo}$ , $\overline{Bi}$ and $\overline{{{K}_{F/P}}}$  as acceptable conservative estimates and the effects of uncontrolled variations and uncertainties can be captured with the distributions of $F{{o}^{*}}$ , $B{{i}^{*}}$ and $K_{F/P}^{*}$ or by choosing them equal to unity.

#### 7.4.2.1      Typical probabilistic migration kinetics

Without loss of generality for estimating multivariate distributions of concentrations, the principles of composition are illustrated in [Figure 24](#figure24) for the dimensionless migration kinetics from monolayer materials (see. §). In this example, $F{{o}^{*}}$ is the only random variable and all other parameters remain fixed (not distributed). The distribution of the dimensionless concentration ${{\bar{v}}^{*}}$is inferred from Eq XXX with $p=1$ (strictly monotonic curve). It can be viewed as the projection of the support of $Fo=\overline{Fo}F{{o}^{*}}$ onto the curve ${{\bar{v}}^{*}}=f\left( F{{o}^{*}} \right)$. When the transformation is repeated for different values for $\overline{Fo}$, the upper and lower percentiles of the migration kinetics can be interpolated continuously with $Fo$. Such curves are not accessible to the direct simulation and are not parallel to the deterministic or 50th percentile. The vertical distance between the extreme percentiles represent the resulting uncertainty in the amount transferred according to the original dispersion in $Fo$values. It is worthwhile noticing that dispersion increases with time but that its effect on ${{\bar{v}}^{*}}$ is decreasing after $Fo>\frac{\sqrt{\pi }}{2}$ . Close to equilibrium, $Fo$ceases to have a significant effect.



<img src="./src/fig24.png" alt="fig24" style="zoom: 80%;" />                  

<a name="figure24"></a>**Figure 24**. Probabilistic modeling of the contamination from a monolayer material via Eqs.  and  : (a) point distribution for $\overline{Fo}=0.5$; (b) corresponding 10th and 90th percentile curves. The likeliest migration curve corresponding to the maximum probability (mode) of the Fo distribution appears in bold.


#### 7.4.2.2     Effect of Bi and sD

The overall mass transfer resistance, $\Re $, of any compartment (either the food or the packaging layer) is evaluated as $\frac{kl}{D}$ or as $\frac{l}{D}$, whether the chemical affinity of the substance is considered or not. For monolayer materials, the packaging-to-food mass transfer resistance ratio is given by ${{K}_{F/P}}B{{i}^{-1}}$ and $B{{i}^{-1}}$ , respectively. $Bi$ values are expected to be large above 103 with well mixed and low viscous liquids in contact with thick or barrier polymers. Low values ranged between 5 and 103 were observed only in polyolefins in contact with liquids [198](#_ENREF_198). In semi-solids, solids and dry foods in contact with polyolefins and plasticized polymers, migration kinetics linearizes with time [199](#_ENREF_199) and $Bi$ approaches unity.

The effects of $Bi$ and ${{s}_{D}}$ are illustrated in [Figure 25](#figure25). Increasing the mass transfer resistance on food side (*i.e.*, decreasing $Bi$) affects non-linearly the dispersion. As justified by Eq. XX and because the contamination is strictly increasing, the dispersion of the contamination is weighted by the term ${{\left| \frac{d}{dv}{{{\bar{v}}}^{*}}^{-1}\left( v \right) \right|}^{-1}}=\frac{d{{{\bar{v}}}^{*}}}{dFo}\left( Fo \right)$. As low $Bi$ values and large $Fo$ ones lead to low slopes of ${{\bar{v}}^{*}}$vs $Fo$ , the dispersion of concentration values are, as expected, maximum for intermediate concentration values far from the initial and equilibrium states. Increasing ${{s}_{D}}$ modifies strongly the shape of the median kinetics, the upper limits and of the overall distribution in food. The depicted example demonstrates that the migration can be overestimated reliably by considering a likely value for $\overline{{{D}_{p}}}$ and by taking for example the 95^th^ percentile of ${{\bar{v}}^{*}}$along with a proper value of ${{s}_{D}}$. Doing the reverse, calculating ${{\bar{v}}^{*}}$from the 95th percentile of ${{D}_{p}}$ values does not guarantee that the value of the contamination is overestimated in 95% of cases. It is not verified in the case when the slope of ${{\bar{v}}^{*}}$changes rapidly with $Fo$. This is true only when $Bi\to 1$, as shown in [Figure 25a](#figure25).



<img src="./src/fig25.png" alt="fig25" style="zoom: 67%;" />  

<a name="figure25"></a>**Figure 25**. (a) Effect of $Bi$on the dimensionless migration kinetics. (b) Effect of for $Bi\to \infty $ (the percentiles are represented as equivalent kinetics; distributions of ${{\bar{v}}^{*}}$ for $\overline{Fo}$ =0.5, 1,2, 3 and 4.

 

# 8     What will be the future?

------
**Highlights on new trends**

- Because migration modeling is cost effective, its application spreads between countries and beyond its original application domain: food packaging in plastics.
- Migration modeling turn the fate of contamination into forecastable phenomena with causes and responsibilities.
- Migration modeling paves the way for a food safety managed in the cloud.

------

 

## 8.1    Extending the legal scope of migration modeling

A robust validation of the macroscopic equations of mass transfer (transport equations and boundary conditions) has been central to the development of the US legal system authorizing migration modeling in the years 90’. The European system focused during the two past decades on diffusion coefficients with much smaller attention on partition coefficients. In both cases, migration modeling is recognized and well accepted but only for thermoplastics. Its application to thermosets, elastomers, paper, and board is comparatively supported by a much small number of scientific and technical publications. Recently, Chinese regulation adopted migration modeling and equivalent calculations without material restriction. The move is a logical step after the substitution of negative lists by positive ones and the adoption of specific migration limits for a broad range of applications.

At the expense of proper significant safety margins, migration modeling has no limited scope and can cover multiple materials and complex conditions of use of food contact materials. The chief difficulty remains in practice the lack of information on the substances initially present in the materials and their amounts. Paradoxically, migration modeling offers an efficient and secure solution to the delicate problem of the declaration of conformity according to the position of the operator in the supply chain. At each step, the essential information to demonstrate the compliance and safety assessment respectively to intentionally added substances can be encoded into almost anonymous numbers (worst-case values of concentrations, properties, and thresholds), without revealing the exact details of formulations, chemical structures, etc. The amount of data could grow along the supply chain by adding new records (association of substances, materials, logistics data, date of production) and refined scenarios, which correspond to the final application of the packaging or the targeted market (*e.g.*, country, food/pharmaceutical). In a very near future, the whole process could be included in a blockchain resisting to falsification and not requiring any third-party authority. Migration modeling software could process any node and verify rapidly hundreds of constraints corresponding to the integrity of all decisions. Date stamps and default expiration dates could alert in advance operators of the needs of revising either risk assessment or risk management decisions. All the effort related to typing inputs, looking into databases and estimating properties could be performed by scanning a QR code or equivalent. The result could be a new QR code verifiable by local authorities.

After several decades of pioneer developments and validation, the chemical industry, producers, converters, recyclers, food industry, retailers, authorities and consumer associations could turn into interoperable systems promoting both safety and economic efficiency. The definition of common standards (CEN, ISO) and training could be the first step. Respectively to modeling, the legal systems initiated for food contact materials and articles could also be used as templates for future evolutions of the regulation of materials used in cosmetic, pharmaceutical, medical, biotechnological, and clothing products.

The possibility of managing cross-contamination via modeling and simulation is expected to impact rapidly good manufacturing and handling practices. Migration modeling can be used, indeed, to establish the causality between practices and migration regardless their position in the supply chain and the distance between the incrimanted material and the food. Demonstrating causality will bring some form of responsibility irrespectively the material or the practice are intended to be in direct contact with food or not. As an example, the contribution of secondary packaging materials acting as large reservoirs of contaminants can be evaluated, as the foreseen effects of corrective actions, including air renewal in storage places, separation of printed and non-printed materials, separation of laminates from monolayer materials…

.

## 8.2    Extending the capacity of migration modeling

The validations of functional barriers and recycled materials have been the main successful applications of migration modeling. The last advance in migration modeling including multiscale modeling offers the only viable solution to evaluate complex problems met by the food industry and the food packaging supply chain:

1.     NIAS: non-intentionally added substances (no need for standards, analytical methods, hypothetical molecules can be accepted);
2.     cross-contamination between materials at any stage of the supply chain (all configurations can be included);
3.      post-consumer contaminations, including misuse
4.     optimization of decontamination step in mechanical recycling processes
5.     materials and articles with repeated use (no need for long experiments)
6.      materials and devices used with flows (no need of any setup)
7.     materials subjected to aging and long-term storage.

NIAS include hypothetical and unknown substances (*e.g.*, breakdown products), but also known impurities and substances intentionally present or added to third-party materials (printed inks, adhesives, lacquers, overpackaging, secondary packaging). All the cases can be evaluated by combining molecular and migration modeling. The same approaches can be used to optimize decontamination conditions (solvent choice, temperature, duration) in recycling processes.

Test conditions for articles associated with repeated use and flows have not been detailed, and total migration scenarios prevail as a general recommendation in most of the regulations. Migration modeling can be used to evaluate articles and devices (gaskets, hoses, tubings, reservoirs, tanks, conveyors…) in more realistic conditions and all along their lifetime. As an illustration, plasticized materials are prone to release substances only for short-time contacts. Once the surface is depleted in plasticizer(s), the glass transition temperature increases sufficiently to transform the contact surface into a temporarily functional barrier. Any period without liquid in contact (*e.g.*, stopping or cleaning period) redistributes the plasticizer(s) uniformly and causes a new risk of leaching.

Long-term storage of materials and material aging redistribute contaminants and may bring new breakdown products. Aging before use is not considered in current regulations. Modeling can offer a very cost-efficient method to evaluate the risk associated with the redistribution of contaminants from tie layers, printing inks, adhesives, lacquers before they are processed in food contact materials.

 

## 8.3    Bridging migration modeling, safe-by-design,  and ecodesign approaches

This chapter encourages to address safety issues at early stages of the design of food packaging and food contact materials before even they become integrated into a finished product. Instead of checking the compliance on the finished product, additives, designs (shapes, surface-to-volume) and conditions of use (shelf-life, storage) can be used to minimize the risk of contamination and cross-contamination. Safe-by-design approaches[181](#_ENREF_181) are particularly relevant for food products devoted to babies and infants as well as for all applications maximizing the amount of recycled materials. The new generation of simulation tools[200](#_ENREF_200) will bridge safe-by-design approaches with eco-design, integrate 3D simulation of mass transfer and explicit food representation. 

 

## 8.4    Online resources ease risk assessment

Safety decisions supported by calculations, modeling and simulation are only acceptable if they demonstrate that the migration amount, estimated by ${{\hat{C}}_{F}}$,  is below the threshold, ${{C}_{thresh}}$ , prescribed by regulation, good manufacturing practices or other considerations. Otherwise, migration testing (experiments) must be applied. Without any threshold to compare, (**i.e.**, non-evaluated substances), ${{\hat{C}}_{F}}$or its conversion into chronic exposure migration modeling loses its interest. The utility of migration increases when ${{\hat{C}}_{F}}$is not too much overestimated (**i.e.**, too conservative) and when the value of ${{C}_{thresh}}$ is not the default value for a very toxic compound. The development of online resources (computational tools, databases, training, guidance,…) offers migration modeling and risk-assessment readily accessible, transparent, and cost-effective.

### 8.4.1     Lower bounds of toxicological thresholds for non-evaluated substances

For non-evaluated substances or in the absence of specific chemical data, the lowest (safest) value of ${{C}_{thresh}}$ is determined by the concept of toxicological threshold of concern (TTC)[176](#_ENREF_176). It has been recently endorsed jointly by EFSA and WHO[201](#_ENREF_201). It is pragmatic ersatz to identify situations, where migration can be presumed to present no appreciable human health risk (risk of cancer lower than 1 in one million). By excluding high potency carcinogens, *i.e.*, aflatoxin-like, azoxy- or N-nitroso-compounds and benzidines, the TTC value for potentially genotoxic compounds reads:

$$
C_{thresh}^{TTC,genotoxic}=0.0025\text{ }\mu g\cdot kg_{BW}^{-1}\cdot da{{y}^{-1}}\times \frac{\text{body}\ \text{weight}}{\text{daily}\ \text{intake}}
$$

Eq.  relies on a chronic exposure estimate: the same food is in contact with the same packaging, and 1 kg  (EU rule) or 3 kg (US rule) of food is daily consumed by an adult. For an adult of 60 kg (EU rule), the default value is 0.15 µg⋅kg^-1^ of food. The default value for infants (EU rule: infant up to 11 months old, weighing 5 kg and eating 0.75 kg) is 0.0167 µg⋅kg^-1^ of food. For infants under three months of age, additional considerations may be considered when ${{\hat{C}}_{F}}$values are close to $C_{thresh}^{TTC,genotoxic}$, including metabolism, frequency, and duration of exposure.

Similar principles apply to substances, which are not part of the exclusionary categories and without structural alert of chemical-specific genotoxicity data suggesting a potential DNA reactivity. The values of $TTC$ and $C_{thresh}^{TTC}$ are tabulated in Table 10 for adults. The classification of substances relies on the *Cramer decision tree* [202](#_ENREF_202), which has been implemented in the program *Toxtree*, commissioned by the *JRC Computational Toxicology and Modelling,* and in the *OECD QSAR Toolbox*. Both tools require chemical structures to be validated independently. Practical tools and databases are listed in [Table 11](#table11).

<a name="table10"></a>**Table 10**. Thresholds concentration values based on $TTC$ approach for an adult of 60 kg and a daily intake of  1 kg.

| *Type of  $TTC$  value*                    | $TTC$   **in   µg/kg body weight per day** | $C_{thresh}^{TTC}$   in µg/kg food   **(assuming 1 kg daily intake)** |
| ------------------------------------------ | ------------------------------------------ | ------------------------------------------------------------ |
| **With structural alter for genotoxicity** | 0.0025                                     | 0.15                                                         |
| **organophosphates and carbamates**        | 0.3                                        | 18                                                           |
| **Cramer class III**                       | 1.5                                        | 90                                                           |
| **Cramer class II**                        | 9.0                                        | 540                                                          |
| **Cramer class I**                         | 30                                         | 1800                                                         |

 

<a name="table11"></a>**Table 11**. Main tools and databases to derive acceptable toxicological threshold.

| **Type of   tools**         | **Online resources**                                         |
| --------------------------- | ------------------------------------------------------------ |
| **TTC tools**               | [OECD toolbox][OECD], [Toxtree][Toxtree], [TTC and related data and scientific opinion of EFSA][TTC] |
| **Toxicological databases** | [CEFIC LRI Toolbox][CEFIC] including RepDose, FeDTex and CEMAS, [Toxnet/ChemIdPlus ][Toxnet],   [Cosmetic Ingredient Review (CIR) database][CIR] and SCCS opinions, [Council of Europe database][CoE] (not publicly available), [Euopean Chemicals Agency][ECHA] final decisions on compliance checks and testing proposals in REACH registration dossiers, [Joint Research Centre][ESIS] ESIS database, [GESTIS substance database][GESTIS] for Occupational Exposure Limits (OELs), [HPV-Program][HPV], [IARC][IARC] list of carcinogenic substances, [NICNAS][NICNAS] (Australia), [NTP][NTP] US/NIH list of carcinogenic substances,  [OSHA][OSHA] EU list of carcinogenic substances, [Toxline][Toxline] from US NIH, [TSCA][TSCA] from US EPA |
| **General databases**       | [Chemspider][Chemspider], [PubChem][PubChem] |

[OECD]: http://www.oecd.org/chemicalsafety/risk-assessment/theoecdqsartoolbox.htm	"OECD toolbox"
[Toxtree]: https://eurl-ecvam.jrc.ec.europa.eu/laboratories-research/predictive_toxicology/qsar_tools/toxtree	"Toxtree maintained by JRC (Ispra, Italy)"
[TTC]: https://www.efsa.europa.eu/en/consultations/call/181112 "results of consultations of EFSA"
[CEFIC]: http://www.cefic-lri.org/lri-toolbox "CEFIC toolbox and related tools"
[Toxnet]: http://chem.sis.nlm.nih.gov/chemidplus/ "Toxnet (NIH, USA)"
[CIR]: http://www.cir-safety.org/ingredients "public database from CIR"
[CoE]: https://fcm.wiv-isp.be/ "not all publicly available"
[ECHA]: https://echa.europa.eu/addressing-chemicals-of-concern "REACH dossiers"
[ESIS]: http://esis.jrc.ec.europa.eu/ "Joint Research Centre"
[GESTIS]: http://www.dguv.de/dguv/ifa/Gefahrstoffdatenbanken/GESTIS-Stoffdatenbank/index-2.jsp "German GESTIS"
[HPV]:  http://webnet.oecd.org/hpv/ui/SponsoredSubstances.aspx "HPV-Program"
[IARC]: https://monographs.iarc.fr/agents-classified-by-the-iarc/ "list of carcinogenic substances"
[NICNAS]: http://www.nicnas.gov.au/chemical-information "Australia"
[NTP]: https://ntp.niehs.nih.gov/pubhealth/roc/index-1.html "US/NIH list of carcinogenic substances"
[OSHA]: https://osha.europa.eu/en/legislation/directives/directive-2004-37-ec-carcinogens-or-mutagens-at-work "EU list of carcinogenic substances"
[Toxline]: http://toxnet.nlm.nih.gov/cgi-bin/sis/htmlgen?TOXLINE "US NIH"
[TSCA]: http://yosemite.epa.gov/oppts/epatscat8.nsf/reportsearch?openform "US EPA"
[Chemspider]: http://www.chemspider.com "free"
[PubChem]: https://pubchem.ncbi.nlm.nih.gov/






### 8.4.2     Migration modeling tools

Several generations of tools have been developed during the last decades using either analytical or numerical schemes. All current tools are capable to describe transfer across multilayer structures and offer the possibility to simulate variable storage conditions (liquid contact, set-off, etc.). Similar results can be obtained with generic solvers of partial differential equations ([PDE](https://en.wikipedia.org/wiki/Partial_differential_equation)), such as [ANSYS](https://www.ansys.com/), [Cast3M](http://www-cast3m.cea.fr/), [Comsol Multiphysics®](https://www.comsol.com/), [Fluent™](https://www.ansys.com/fr-fr/products/fluids/ansys-fluent), [FreeFem](https://freefem.org/), [Modelica](https://www.modelica.org/), [OpenFoam](https://www.openfoam.com/), etc. The difficulty for multilayer structures is the implementation of the internal conditions defined by Eq. XXX , which require a weak formulation of the [PDE](https://en.wikipedia.org/wiki/Partial_differential_equation) problem. 

<a name="table12"></a>**Table 12**. Notable physics-based software and tools to evaluate migration.

| Type of tools                                          | Description                                                  | License                                                      |
| ------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| **Stand-alone compliance testing programs**            | [Migratest© EXP][migratestexp], [AKTS-SML version 6][aktsml] | *commercial*  (demo available)                               |
| **Compliance testing client/server**                   | Client-server [SFPP3][sfpp3login] application   (*SafeFoodPackaging portal version 3*) to be used by one up to 25 simultaneous users. [SFPP3][sfpp3login] includes all public data of the European task force TF-MATHMOD. [Interactive training][propack] on [SFPP3][sfpp3login] tools with case studies is available in French | *freeware*,   partly open-source, online access or standalone  installation. |
| **Expandable preventive and safe-by-design   tools**   | [FMECAengine][fmecaeng] combined with *key2key()* language enabling to simulate supply chains, complex materials and applications *(includes and expands all features implemented in [SFPP3][sfpp3login])* | *fully open-source*  (written in Matlab®, Octave language)   |
| **Online**   **databases**                             | [Thermophysical properties of polymers][poldb]; [database of diffusion and partition coefficients][dbDK] of the European task force TF-MATHMOD | *free access*                                                |
| **Guidance**                                           | [EU guidance][EUrules], [US guidance][USrules], generic guidance from [INRA][gINRA], [Food Packaging Forum][gpackforum] | *free access*                                                |
| **MOOC:   Massive Open Online Course (USA) – 7 hours** | MOOC lectures from the School of Packaging (Michigan State University) :  [overview][MSU0], [diffusion coefficients][MSU1]: $D_{i,P}$, [Partition coefficients][MSU2]: $K_{i,F/P}$, [Migration modeling and decision-making workshop][MSU3] | *no registration*                                            |
| **MOOC (EU) > 50 hours**                               | The [European project ERASMUS+ Fitness][Fitness0] “*Food packaging open courseware for higher education   and staff of companies*” is preparing a MOOC on all aspects of   packaging design including compliance testing, food safety, and risk assessment: [demo site][Fitness1] | *to be released by late 2019*                                |

[migratestexp]: https://www.fabes-online.de/en/software-en/migratest-exp/ "commercial software (demo available)"
[aktssml]: https://www.akts.com/sml-diffusion-migration-multilayer-packaging/download-diffusion-prediction-software.html "commercial software (demo available)"
[sfpp3login]: http://sfpp3.agroparistech.fr:443/cgi-bin/login.cgi "INRA website (account: demouser, password: inramig)"
[propack]: http://rmt-propackfood.actia-asso.eu/ "registration required (free, in French)"
[fmecaeng]: https://github.com/ovitrac/FMECAengine "opensource"
[poldb]: http://polymerdatabase.com/ "external website"
[dbDK]: http://modmol.agroparistech.fr/Database/ "INRA website (account: demouser, password: inramig)"
[gINRA]: http://modmol.agroparistech.fr/home/ "INRA website (account: demouser, password: inramig)"
[gpackforum]: https://www.foodpackagingforum.org/food-packaging-health/migration-modeling, "Food Packaging Forum website"
[EUrules]: http://publications.jrc.ec.europa.eu/repository/handle/JRC98028 "JRC website"
[USrules]: https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P100BCMB.TXT "EPA website"
[Fitness0]: http://www2.agroparistech.fr/FITNess-Project.html "overview (no registration)"
[Fitness1]: http://mol15.agroparistech.fr/fitness/lectures/html/ "main menu (no registration)"
[MSU0]: https://www.canr.msu.edu/cpis/research/downloadable_presentations "overview (no registration)"
[MSU1]: https://mediaspace.msu.edu/media/Dr.+Olivier+Vitrac+presentsA+Diffusion+coefficients+of+organic+solutes+in+polymersA/1_zz20dgt9 "Di,P (no registration)"
[MSU2]: https://mediaspace.msu.edu/media/Dr.+Olivier+Vitrac+presentsA+An+atomistic+Flory-Huggins+formulation+for+the+tailored+prediction+of+activity+and+partition+coefficients/1_uzi6h91k "Ki,F/P, Ki,P1/P2 (no registration)"
[MSU3]: https://mediaspace.msu.edu/media/WorkshopA+Prediction+of+the+migrationA+beyond+conventional+estimates*/1_won1m7aw "workshop (no registration)"


The notable migration modeling software are listed in Table 12 with their corresponding license (open-source, freeware and commercial). It is worth noticing that they implement Fickian diffusion and that they assume uniform and isotropic distribution within each layer. These tools need to be combined with internal or external databases and mathematical relationships in order to take into account the nature of the migrant and the physical properties of the polymer. Open-source software and devoted programming language, as [FMECAengine](https://github.com/ovitrac/FMECAengine) and *key2key()* language, are part of a global effort to foster interoperability, collaborative open-source projects and free templates for common problems. The European Committee for Standardization already identified that “a standardized terminology will improve future exchanges among experts in the entire area of materials modeling, facilitate the exchange with industrial end-users and experimentalists and reduce the barrier utilizing materials modeling” [203](#_ENREF_203).

 

# 9    References

 

<a name="_ENREF_1">\1</a>.         ESA Outgassing data. http://esmat.esa.int/Services/outgassing_data/outgassing_data.html (access on Nov 22th 2018), 

<a name="_ENREF_2">\2</a>.           NASA Outgassing Data for Selecting Spacecraft Materials. https://outgassing.nasa.gov/ (accessed on March 29, 2019), 

<a name="_ENREF_3">\3</a>.           Kortenkamp, A.; Backhaus, T.; Faust, M. State of the Art Report on Mixture Toxicity - Final Report Executive Summary. **2009**.

<a name="ENREF_4">\4</a>           Lee, D.-H. Evidence of the Possible Harm of Endocrine-Disrupting Chemicals in Humans: Ongoing Debates and Key Issues. *Endocrinology and Metabolism (EnM)* **2018,** 33 (1), 44-52.

<a name="ENREF_5">\5</a>           Kortenkamp, A. Ten Years of Mixing Cocktails: A Review of Combination Effects of Endocrine-Disrupting Chemicals. *Environmental Health Perspectives* **2007,** 115 (Suppl 1), 98-105.

<a name="ENREF_6">\6</a>           Wild, C. P. Complementing the Genome with an “Exposome”: The Outstanding Challenge of Environmental Exposure Measurement in Molecular Epidemiology. *Cancer Epidemiology Biomarkers &amp;amp; Prevention* **2005,** 14 (8), 1847.

<a name="ENREF_7">\7</a>           Eicher, A.; Biedermann, M.; Zurfluh, M.; Grob, K. Migration by ‘direct’ or ‘indirect’ food contact? ‘Dry’ and ‘wetting’ foods? Experimental data for ‘touching’ contact of dry foods with paper and board. *Food Additives & Contaminants: Part A* **2015,** 32 (1), 110-119.

<a name="ENREF_8">\8</a>           Muhamad, M. S.; Salim, M. R.; Lau, W. J.; Yusop, Z. A review on bisphenol A occurrences, health effects and treatment process via membrane technology for drinking water. *Environmental Science and Pollution Research* **2016,** 23 (12), 11549-11567.

<a name="ENREF_9">\9</a>           Biedermann-Brem, S.; Biedermann, M.; Grob, K. Required barrier efficiency of internal bags against the migration from recycled paperboard packaging into food: a benchmark. *Food Additives & Contaminants: Part A* **2016,** 33 (4), 725-740.

<a name="ENREF_10">\10</a>        Lommatzsch, M.; Richter, L.; Biedermann-Brem, S.; Biedermann, M.; Grob, K.; Simat, T. J. Functional barriers or adsorbent to reduce the migration of mineral oil hydrocarbons from recycled cardboard into dry food. *European Food Research and Technology* **2016,** 242 (10), 1727-1733.

<a name="ENREF_11">\11</a>        Paul, U. C.; Fragouli, D.; Bayer, I. S.; Mele, E.; Conchione, C.; Cingolani, R.; Moret, S.; Athanassiou, A. Mineral oil barrier sequential polymer treatment for recycled paper products in food packaging. *Materials Research Express* **2017,** 4 (1), 015501.

<a name="ENREF_12">\12</a>        Nriagu, J. O., *Lead and Lead Poisoning in Antiquity*. John Wiley: New-York, USA, 1983; p 437.

<a name="ENREF_13">\13</a>        Patterson, C. C.; Shirahata, H.; Ericson, J. E. Lead in ancient human bones and its relevance to historical developments of social problems with lead. *Science of The Total Environment* **1987,** 61, 167-200.

<a name="ENREF_14">\14</a>        Delile, H.; Blichert-Toft, J.; Goiran, J.-P.; Keay, S.; Albarède, F. Lead in ancient Rome’s city waters. *Proceedings of the National Academy of Sciences* **2014,** 111 (18), 6594.

<a name="ENREF_15">\15</a>        Doumerc, A.; de Leymarie, L., *Législation française et étrangère concernant les falsifications alimentaires*. Rueff et Cie, Editeurs: Paris, 1895; p 258.

<a name="ENREF_16">\16</a>        Hamel, G. T. Méfions-nous des conserves étrangères. *La conserve alimentaire* **1910,** 19, 290-291.

<a name="ENREF_17">\17</a>        Rall, D. P. The Invisible Pollution. *New England Journal of Medicine* **1972,** 287 (22), 1146-1147.

<a name="ENREF_18">\18</a>        Jaeger, R. J.; Rubin, R. J. Plasticizers from Plastic Devices: Extraction, Metabolism, and Accumulation by Biological Systems. *Science* **1970,** 170 (3956), 460.

<a name="ENREF_19">\19</a>        Jaeger, R. J.; Rubin, R. J. Migration of a Phthalate Ester Plasticizer from Polyvinyl Chloride Blood Bags into Stored Human Blood and Its Localization in Human Tissues. *New England Journal of Medicine* **1972,** 287 (22), 1114-1118.

<a name="ENREF_20">\20</a>        Lawrence, W. H.; Turner, J. E.; Autian, J. Reevaluation of plastic tubings currently used in medical and paramedical applications. *Journal of Biomedical Materials Research* **1969,** 3 (2), 291-303.

<a name="ENREF_21">\21</a>        Kestel'man, V. N. e. a. The corrosion of polymeric material under the conditions of the microbiological synthesis of enzymes. *Int. Biodetn Bull.* **1972,** 8 (1), 15-19.

<a name="ENREF_22">\22</a>        Figge, K. Migration of components from plastics-packaging materials into packed goods — test methods and diffusion models. *Progress in Polymer Science* **1980,** 6 (4), 187-252.

<a name="ENREF_23">\23</a>        Giacin, J. R. EVALUATION OF PLASTICS PACKAGING MATERIALS FOR FOOD PACKAGING APPLICATIONS: FOOD SAFETY CONSIDERATIONS. *Journal of Food Safety* **1980,** 2 (4), 257-276.

<a name="ENREF_24">\24</a>        Till, D. E.; Ehntholt, D. J.; Reid, R. C.; Schwartz, P. S.; Schwope, A. D.; Sidman, K. R.; Whelan, R. H. Migration of styrene monomer from crystal polystyrene to foods and food simulating liquids. *Industrial & Engineering Chemistry Fundamentals* **1982,** 21 (2), 161-168.

<a name="ENREF_25">\25</a>        Till, D. E.; Ehntholt, D. J.; Reid, R. C.; Schwartz, P. S.; Sidman, K. R.; Schwope, A. D.; Whelan, R. H. Migration of BHT antioxidant from high density polyethylene to foods and food simulants. *Industrial & Engineering Chemistry Product Research and Development* **1982,** 21 (1), 106-113.

<a name="ENREF_26">\26</a>        Schwope, A. D.; Till, D. E.; Ehntholt, D. J.; Sidman, K. R.; Whelan, R. H.; Schwartz, P. S.; Reid, R. C. Migration of BHT and Irganox 1010 from low-density polyethylene (LDPE) to foods and food-simulating liquids. *Food and Chemical Toxicology* **1987,** 25 (4), 317-326.

<a name="ENREF_27">\27</a>        Schwope, A. D.; Till, D. E.; Ehntholt, D. J.; Sidman, K. R.; Whelan, R. H.; Schwartz, P. S.; Reid, R. C. Migration of Irganox 1010 from ethylene-vinyl acetate films to foods and food-simulating liquids. *Food and Chemical Toxicology* **1987,** 25 (4), 327-330.

<a name="ENREF_28">\28</a>        Gandek, T. P., *PhD Thesis. Migration of Phenolic Antioxidants from Polyolefins to Aqueous Media with Application to Indirect Food Additive Migration*. Massachusetts Institute of Technology, Department of Chemical Engineering: Cambridge, MA, USA, 1986; p 519.

<a name="ENREF_29">\29</a>        Gandek, T. P.; Hatton, T. A.; Reid, R. C. Batch extraction with reaction: phenolic antioxidant migration from polyolefins to water. 1. Theory. *Industrial & Engineering Chemistry Research* **1989,** 28 (7), 1030-1036.

<a name="ENREF_30">\30</a>        Gandek, T. P.; Hatton, T. A.; Reid, R. C. Batch extraction with reaction: phenolic antioxidant migration from polyolefins to water. 2. Experimental results and discussion. *Industrial & Engineering Chemistry Research* **1989,** 28 (7), 1036-1045.

<a name="ENREF_31">\31</a>        Fornasiero, F.; Olaya, M. M.; Wagner, I.; Brüderle, F.; Prausnitz, J. M. Solubilities of nonvolatile solutes in polymers from molecular thermodynamics. *AIChE Journal* **2002,** 48 (6), 1284-1291.

<a name="ENREF_32">\32</a>        Merrill, R. A. FOOD SAFETY REGULATION:Reforming the Delaney Clause. *Annual Review of Public Health* **1997,** 18 (1), 313-340.

<a name="ENREF_33">\33</a>        National Research Council (US) Committee on Diet, N., and Cancer,, *Diet, Nutrition, and Cancer*. National Academies Press Washington, DC, USA, 1982; Vol. Available from: https://www.ncbi.nlm.nih.gov/books/NBK216644/, p 496.

<a name="ENREF_34">\34</a>        Roberts, H. R. Food Additives - A Study in the Evolution of Safety. *Food Drug Cosmetic Law Journal* **1976,** 31, 404.

<a name="ENREF_35">\35</a>        EC. Commission Directive 90/128/EEC of 23 February 1990 relating to plastic materials and articles intended to come into contact with foodstuffs. *Official Journal of the European Union* **1990,** L75, 19-40.

<a name="ENREF_36">\36</a>        EC. Commission Directive 2002/72/EC of 6 August 2002 relating to plastic materials and articles intended to come into contact with foodstuffs. *Official Journal of the European Union* **2002,** L220, 18-56.

<a name="ENREF_37">\37</a>        EC. Commission Regulation (EU) No 10/2011 of 14 January 2011 on plastic materials and articles intended to come into contact with food Text with EEA relevance. *Official Journal of the European Union* **2011,** L12, 1-89.

<a name="ENREF_38">\38</a>        CFDA. National Standard of The People’s Republic Of China GB 9685-2016: National Food Safety Standard Standard for Uses of Additives in Food Contact Materials and Their Products. *National Health and Family Planning Commission* **2016,** Available from: http://samr.cfda.gov.cn/WS01/CL0087/166631.html.

<a name="ENREF_39">\39</a>        FDA. Inventory of Effective Food Contact Substance (FCS) Notifications. **2018,** available from: https://www.fda.gov/food/ingredientspackaginglabeling/packagingfcs/notifications/default.htm.

<a name="ENREF_40">\40</a>        DGCCRF. Arrêté du 9 novembre 1994 relatif aux matériaux et objets en caoutchouc au contact des denrées, produits et boissons alimentaires. *Journal Officiel de la République Française* **1995,**  (2 décembre 1994), 17029-17035.

<a name="ENREF_41">\41</a>        EuPIA Good Manufacturing Practice (GMP) Printing Inks for Food Contact Materials, 4th completely revised version, EuPIA, March 2016. http://www.eupia.org/index.php?id=29 ((access on March 29, 2019)), 

<a name="ENREF_42">\42</a>        EDI Swiss Ordinance of the FDHA on materials and articles intended to come into contact with  foodstuffs (SR 817.023.21) , Eidgenössisches Departement des Innern (EDI), 16.12.2016. https://www.admin.ch/opc/de/classified-compilation/20143393/index.html ((accessed on March 29, 2019)), 

<a name="ENREF_43">\43</a>        EuPIA EuPIA Guidance for Risk Assessment of Non-Intentionally Added Substances (NIAS) and Non-Listed Substances (NLS) in printing inks for food contact materials. http://www.eupia.org/uploads/tx_edm/2016-02-24_EuPIA_Info_Note_Inks_and_Coatings_for_High_Temperature_Applications.pdf (accessed on March 29, 2019), 

<a name="ENREF_44">\44</a>        COE *Paper and board materials and articles intended to come into contact with foodstuffs - version 4 - 12.02.2009 -* [*https://rm.coe.int/16804e4794*](https://rm.coe.int/16804e4794); 2009; p 148.

<a name="ENREF_45">\45</a>        CEPI *Industry Guideline for the Compliance of Paper & Board Materials and Articles for Food Contact. Issue 2 September 2012 -* [*http://www.cepi.org/system/files/public/documents/publications/foodcontact/2012/Industry%20guideline-updated2012final.pdf*](http://www.cepi.org/system/files/public/documents/publications/foodcontact/2012/Industry guideline-updated2012final.pdf); Confederation of European Paper Industries: Brussels, Belgium, 2012; p 30.

<a name="ENREF_46">\46</a>        BfR *XXXVI. Papiere, Kartons und Pappen für den Lebensmittelkontakt -* [*https://bfr.ble.de/kse/faces/resources/pdf/360.pdf*](https://bfr.ble.de/kse/faces/resources/pdf/360.pdf); BfR: Berlin, Gemany, 2017; p 16.

<a name="ENREF_47">\47</a>        EFSA. Scientific Opinion on Mineral Oil Hydrocarbons in Food. *EFSA Journal* **2012,** 10 (6), 2704 (185 p).

<a name="ENREF_48">\48</a>        Bradley, E. L.; Burden, R. A.; Bentayeb, K.; Driffield, M.; Harmer, N.; Mortimer, D. N.; Speck, D. R.; Ticha, J.; Castle, L. Exposure to phthalic acid, phthalate diesters and phthalate monoesters from foodstuffs: UK total diet study results. *Food Additives & Contaminants: Part A* **2013,** 30 (4), 735-742.

<a name="ENREF_49">\49</a>        Bradley, E. L.; Burden, R. A.; Leon, I.; Mortimer, D. N.; Speck, D. R.; Castle, L. Determination of phthalate diesters in foods. *Food Additives & Contaminants: Part A* **2013,** 30 (4), 722-734.

<a name="ENREF_50">\50</a>        MAFF, U. Phthalates in paper and board packaging. *MAFF Food Surveillance Information Sheet.* **1995,** Number 60.

<a name="ENREF_51">\51</a>        MAFF, U. Phthalates in food. *MAFF Food Surveillance Information Sheet.* **1996,** Number 82, [https://webarchive.nationalarchives.gov.uk/20081230231745/http://archive.food.gov.uk/maff/archive/food/infsheet/1996/no82/82phthal.htm](https://webarchive.nationalarchives.gov.uk/20081230231745/http:/archive.food.gov.uk/maff/archive/food/infsheet/1996/no82/82phthal.htm).

<a name="ENREF_52">\52</a>        MAFF, U. Phthalates in infant formulae. *MAFF Food Surveillance Information Sheet.* **1996,** Number 83.

<a name="ENREF_53">\53</a>        MAFF, U. Phthalates in infant formulae – follow-up survey. *MAFF Food Surveillance Information Sheet.* **1998,** Number 168.

<a name="ENREF_54">\54</a>        EFSA. Opinion of the Scientific Panel on food additives, flavourings, processing aids and materials in contact with food (AFC) related to 2,2-bis(4-hydroxyphenyl)propane bis(2,3-epoxypropyl)ether (Bisphenol A diglycidyl ether, BADGE). REF. No 13510 and 39700. *EFSA Journal* **2004,** 2 (7), 86.

<a name="ENREF_55">\55</a>        Six, T.; Feigenbaum, A.; Riquet, A. M. Mechanism of migration from agglomerated cork stoppers: I. An electron spin resonance investigation. *Journal of Applied Polymer Science* **2002,** 83 (12), 2644-2654.

<a name="ENREF_56">\56</a>        Six, T.; Feigenbaum, A. Mechanism of migration from agglomerated cork stoppers. Part 2: Safety assessment criteria of agglomerated cork stoppers for champagne wine cork producers, for users and for control laboratories. *Food Additives & Contaminants* **2003,** 20 (10), 960-971.

<a name="ENREF_57">\57</a>        Campanella, G.; Ghaani, M.; Quetti, G.; Farris, S. On the origin of primary aromatic amines in food packaging materials. *Trends in Food Science & Technology* **2015,** 46 (1), 137-143.

<a name="ENREF_58">\58</a>        Mutsuga, M.; Yamaguchi, M.; Kawamura, Y. Quantification of isocyanates and amines in polyurethane foams and coated products by liquid chromatography–tandem mass spectrometry. *Food Science & Nutrition* **2014,** 2 (2), 156-163.

<a name="ENREF_59">\59</a>        Athenstädt, B.; Fünfrocken, M.; Schmidt, T. C. Migrating components in a polyurethane laminating adhesive identified using gas chromatography/mass spectrometry. *Rapid Communications in Mass Spectrometry* **2012,** 26 (16), 1810-1816.

<a name="ENREF_60">\60</a>        Castle, L.; Sharman, M.; Gilbert, J. Gas chromatographic-mass spectrometric determination of epoxidized soybean oil contamination of foods by migration from plastic packaging. *Journal - Association of Official Analytical Chemists* **1988,** 71 (6), 1183-6.

<a name="ENREF_61">\61</a>        EC *4-Nonylphenol (branched) and nonylphenol -  Risk Assessment - EUR 20387 EN -* [*https://www.bfr.bund.de/cm/343/4_nonylphenol_und_nonylphenol.pdf*](https://www.bfr.bund.de/cm/343/4_nonylphenol_und_nonylphenol.pdf); EUR 20387 EN; European Communities: Luxembourg, 2002; p 14.

<a name="ENREF_62">\62</a>        EC. DIRECTIVE 2003/53/EC OF THE EUROPEAN PARLIAMENT AND OF THE COUNCIL of 18 June 2003 amending for the 26th time Council Directive 76/769/EEC relating to restrictions on the marketing and use of certain dangerous substances and preparations (nonylphenol, nonylphenol ethoxylate and cement). *Official Journal of the European Union* **2003,** L178, 24-27.

<a name="ENREF_63">\63</a>        ECHA. Entry 46: Nonylphenol - Conditions of restriction. *Annex XVII - Restrictions on the manufacture, placing on the market and use of certain dangerous substances, mixtures and articles* **2016,** https://echa.europa.eu/fr/substances-restricted-under-reach/-/dislist/details/0b0236e1807e2b6e.

<a name="ENREF_64">\64</a>        EFSA. Opinion of the Scientific Panel on food additives, flavourings, processing aids and materials in contact with food (AFC) related to Semicarbazide in food. *EFSA Journal* **2005,** 3 (6), 219.

<a name="ENREF_65">\65</a>        Stadler, R. H.; Mottier, P.; Guy, P.; Gremaud, E.; Varga, N.; Lalljie, S.; Whitaker, R.; Kintscher, J.; Dudler, V.; Read, W. A.; Castle, L. Semicarbazide is a minor thermal decomposition product of azodicarbonamide used in the gaskets of certain food jars. *Analyst* **2004,** 129 (3), 276-281.

<a name="ENREF_66">\66</a>        EC. Commission Directive 2004/1/EC of 6 January 2004 amending Directive 2002/72/EC as regards the suspension of the use of azodicarbonamide as blowing agent (Text with EEA relevance). *Official Journal of the European Union* **2004,** C007, 45-46.

<a name="ENREF_67">\67</a>        EFSA. Opinion of the Scientific Panel on food additives, flavourings, processing aids and materials in contact with food (AFC) related to 2-Isopropyl thioxanthone (ITX) and 2-ethylhexyl-4-dimethylaminobenzoate (EHDAB) in food contact materials. *EFSA Journal* **2005,** 3 (12), 293.

<a name="ENREF_68">\68</a>        Johns, S. M.; Gramshaw, J. W.; Castle, L.; Jickells, S. M. Studies on functional barriers to migration. 1. Transfer of benzophenone from printed paperboard to microwaved food. *Deutsche Lebensmittel-Rundschau* **1995,** 91 (3), 69-73.

<a name="ENREF_69">\69</a>        EFSA. Opinion of the Scientific Panel on food additives, flavourings, processing aids and materials in contact with food (AFC) related to 2,2-BIS(4-HYDROXYPHENYL)PROPANE. *EFSA Journal* **2006,** 5 (1), 428.

<a name="ENREF_70">\70</a>        Gundert-Remy, U.; Bodin, J.; Bosetti, C.; FitzGerald, R.; Hanberg, A.; Hass, U.; Hooijmans, C.; Rooney, A. A.; Rousselle, C.; van Loveren, H.; Wölfle, D.; Barizzone, F.; Croera, C.; Putzu, C.; Castoldi, A. F. Bisphenol A (BPA) hazard assessment protocol. *EFSA Supporting Publications* **2017,** 14 (12), 1354E.

<a name="ENREF_71">\71</a>        Castle, L.; Kelly, M.; Gilbert, J. Migration of mineral hydrocarbons into foods. 2. Polystyrene, ABS, and waxed paperboard containers for dairy products. *Food Additives and Contaminants* **1993,** 10 (2), 167-174.

<a name="ENREF_72">\72</a>        EFSA. Scientific Opinion on Mineral Oil Hydrocarbons in Food. *EFSA Journal* **2012,** 10 (6), 2704.

<a name="ENREF_73">\73</a>        EC. Rapid Alert System for Food and Feed - RASFF. https://webgate.ec.europa.eu/rasff-window/portal/index.cfm?event=notificationsList (accessed on Dec 31, 2018), 

<a name="ENREF_74">\74</a>        Begley, T. H. Methods and approaches used by FDA to evaluate the safety of food packaging materials. *Food Additives & Contaminants* **1997,** 14 (6-7), 545-553.

<a name="ENREF_75">\75</a>        Limm, W.; Hollifield, H. C. Modelling of additive diffusion in polyolefins. *Food Additives & Contaminants* **1996,** 13 (8), 949-967.

<a name="ENREF_76">\76</a>        Baner, A.; Brandsch, J.; Franz, R.; Piringer, O. The application of a predictive migration model for evaluating the compliance of plastic materials with European food regulations. *Food Additives & Contaminants* **1996,** 13 (5), 587-601.

<a name="ENREF_77">\77</a>        Baner, A.; Franz, R.; Piringer, O. Alternative methods for the determination and evaluation of migration potential from polymeric food contact materials. *Deutsche Lebensmittel-Rundschau* **1994,** 90 (6), 181-185.

<a name="ENREF_78">\78</a>        Baner, A.; Franz, R.; Piringer, O. Alternative methods for the determination and evaluation of migration potential from polymeric food contact materials. *Deutsche Lebensmittel-Rundschau* **1994,** 90 (5), 137-143.

<a name="ENREF_79">\79</a>        K., H.; O., P. *Evaluation of migration models to used under Directive 90/128/EEC. Final report contract SMT4-CT98-7513 - EUR 20604 EN -* [*https://publications.europa.eu/en/publication-detail/-/publication/07919045-cba5-4ad2-9460-c17c56e2f329/language-en*](https://publications.europa.eu/en/publication-detail/-/publication/07919045-cba5-4ad2-9460-c17c56e2f329/language-en); European Commission, Directorate General for Research: Luxembourg, 2002; p 249.

<a name="ENREF_80">\80</a>        Hoekstra, E.; Brandsch, R.; Dequatre, C.; Mercea, P.; Milana, M. R.; Stoermer, A.; Trier, X.; Vitrac, O.; Schaefer, A.; Simoneau, C. *Practical guidelines on the application of migration modelling for the estimation of specific migration*; European Commission, Joint Research Centre, Institute for Health and Consumers Protection: 2015.

<a name="ENREF_81">\81</a>        Schwope, A. D.; Goydan, R.; Reid, R. *Methodology for Estimating the Migration of Additives and Impurities from Polymeric Materials - EPA 560/5-85-015*; Office of Pesticides and Toxic Substances, U.S. Environmental Protection Agency (EPA): Washington DC, USA - https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P100BCMB.TXT, 1990; p 148.

<a name="ENREF_82">\82</a>        Vitrac, O. An atomistic Flory-Huggins formulation for the tailored prediction of activity and partition coefficients. https://mediaspace.msu.edu/media/Dr.+Olivier+Vitrac+presentsA+An+atomistic+Flory-Huggins+formulation+for+the+tailored+prediction+of+activity+and+partition+coefficients/1_uzi6h91k (accessed on Dec 31, 2019), 

<a name="ENREF_83">\83</a>        Vitrac, O. Diffusion coefficients of organic solutes in polymers: new perspectives of prediction. https://mediaspace.msu.edu/media/Dr.+Olivier+Vitrac+presentsA+Diffusion+coefficients+of+organic+solutes+in+polymersA/1_zz20dgt9 (accessed on Dec 31, 2019), 

<a name="ENREF_84">\84</a>        Vitrac, O. New computational methodologies as enablers of safer food contact articles. https://www.youtube.com/watch?v=v1SMAUy6XJs&t=312s (accessed on Dec 31, 2018), 

<a name="ENREF_85">\85</a>        Vergnaud, J. M., *Liquid transport processes in polymeric materials: modeling and industrial applications*. Prentice Hall International LTD: London, UK, 1991; p 362.

<a name="ENREF_86">\86</a>        Piringer, O. G.; Baner, A. L., *Plastic packaging materials for food: barrier function, mass transport, quality assurance, and legislation*. Wiley-VCH Verlag: Weinheim, Germany, 2000; p 576.

<a name="ENREF_87">\87</a>        Barnes, K.; Sinclair, R.; Watson, D., *Chemical Migration and Food Contact Materials*. CRC Press LCL: Boca-Raton, FL, USA, 2006; p 464.

<a name="ENREF_88">\88</a>       Vergnaud, J. M.; Rosca, I. D., *Assessing Food Safety of Polymer Packaging*. Rapra Technology: Shawbury, UK, 2006; p 273.

<a name="ENREF_89">\89</a>        Piringer, O. G.; Baner, A. L., *Plastic Packaging Interactions with Food and Pharmaceuticals*. Wiley-VCH Verlag: Weinheim, Germany, 2008; p 614.

<a name="ENREF_90">\90</a>        Singh, P.; Wani, A. A.; Langowski, H. C., *Food Packaging Materials: Testing & Quality Assurance*. CRC Press: Boca-Raton, USA, 2017; p 344.

<a name="ENREF_91">\91</a>        Poças, M. F.; Oliveira, J. C.; Oliveira, F. A. R.; Hogg, T. A Critical Survey of Predictive Mathematical Models for Migration from Packaging. *Critical Reviews in Food Science and Nutrition* **2008,** 48 (10), 913-928.

<a name="ENREF_92">\92</a>        Helmroth, E.; Rijk, R.; Dekker, M.; Jongen, W. Predictive modelling of migration from packaging materials into food products for regulatory purposes. *Trends in Food Science & Technology* **2002,** 13 (3), 102-109.

<a name="ENREF_93">\93</a>        Arvanitoyannis, I. S.; Bosnea, L. Migration of substances from food packaging materials to foods. *Critical Reviews in Food Science and Nutrition* **2004,** 44 (2), 63-76.

<a name="ENREF_94">\94</a>        Lau, O.-W.; Wong, S.-K. Contamination in food from packaging material. *Journal of Chromatography A* **2000,** 882 (1–2), 255-270.

<a name="ENREF_95">\95</a>        Gillet, G.; Vitrac, O.; Tissier, D.; Saillard, P.; Desobry, S. Development of decision tools to assess migration from plastic materials in contact with food. *Food Additives and Contaminants Part a-Chemistry Analysis Control Exposure & Risk Assessment* **2009,** 26 (12), 1556-1573.

<a name="ENREF_96">\96</a>        Vitrac, O.; Challe, B.; Leblanc, J. C.; Feigenbaum, A. Contamination of packaged food by substances migrating from a direct-contact plastic layer: Assessment using a generic quantitative household scale methodology. *Food Additives & Contaminants* **2007,** 24 (1), 75-94.

<a name="ENREF_97">\97</a>        Vitrac, O.; Leblanc, J.-C. Consumer exposure to substances in plastic packaging. I. Assessment of the contribution of styrene from yogurt pots. *Food Additives and Contaminants* **2007,** 24 (2).

<a name="ENREF_98">\98</a>        Vitrac, O.; Hayert, M., Design of safe packaging materials under uncertainty. In *Chemical Engineering Research Trends*, Berton, L. P., Ed. Nova Science Publishers: New-York, USA, 2007; pp 251-292.

<a name="ENREF_99">\99</a>        Robertson, G. L., *Food Packaging: Principles and Practice, Third Edition*. CRC Press: Boca-Raton, USA, 2016; p 686.

<a name="ENREF_100">\100</a>      Cussler, E. L., *Diffusion: Mass Transfer in Fluid Systems. 3rd Edition.* Cambridge University Press: Cambridge, UK, 2009; p 631.

<a name="ENREF_101">\101</a>      Vieth, W. R., *Diffusion in and Through Polymers: Principles and Applications*. Carl Hanser GmbH: Munich, Germany, 1991; p 330.

<a name="ENREF_102">\102</a>      Crank, J., *The Mathematics of Diffusion* Oxford University Press: Bristol, UK, 1975; p 414.

<a name="ENREF_103">\103</a>      Stastna, J.; De Kee, D., *Transport Properties in Polymers*. CRC Press: Boca-Raon, USA, 1995; p 303.

<a name="ENREF_104">\104</a>      Mehrer, H., *Diffusion in Solids: Fundamentals, Methods, Materials, Diffusion-Controlled Processes*. Springer-Verlag: Heidelberg, Germany, 2010.

<a name="ENREF_105">\105</a>      Neogi, P., *Diffusion in Polymers*. Marcel Dekker, Inc: New-Yourk, USA, 1996.

<a name="ENREF_106">\106</a>      Ben-Naim, A. Y., *Solvation Thermodynamics*. Springer Science: New-York, USA, 2013; p 246.

<a name="ENREF_107">\107</a>      Vrentas, J. S.; Vrentas, C. M., *Diffusion and Mass Transfer*. CRC Press: Boca-Raton, FL, USA, 2013; p 619.

<a name="ENREF_108">\108</a>      EC. COMMISSION DECISION of 28 January 1997 establishing the identification system for packaging materials pursuant to European Parliament and Council Directive 94/62/EC on packaging and packaging waste (97/ 129/EC). *Official Journal of the European Communities* **1997,** L50, 28-31.

<a name="ENREF_109">\109</a>      GB. GB/T 16288-2008 : Marking of Plastics Products. *NATIONAL STANDARD OF THE PEOPLE’S REPUBLIC OF CHINA* **2008,**  (ICS 83. 140, Y 28), 1-17.

<a name="ENREF_110">\110</a>      van Krevelen, D. W.; te Nijenhuis, K., *Properties of Polymers: Their Correlation with Chemical Structure; their Numerical Estimation and Prediction from Additive Group Contributions*. Elsevier Science: Amsterdam, NL, 2009.

<a name="ENREF_111">\111</a>      Chemical Retrieval on the Web Polymer Properties Database. http://polymerdatabase.com/ (accessed on March 31, 2019), 

<a name="ENREF_112">\112</a>      Gillet, G.; Vitrac, O.; Desobry, S. Prediction of Solute Partition Coefficients between Polyolefins and Alcohols Using a Generalized Flory−Huggins Approach. *Industrial & Engineering Chemistry Research* **2009,** 48 (11), 5285-5301.

<a name="ENREF_113">\113</a>      Gillet, G.; Vitrac, O.; Desobry, S. A Fast Method to Assess the Composition of a Polyolefin: An Application to Compliance Testing of Food Contact Materials. *Journal of Applied Polymer Science* **2011,** 119 (3), 1492-1515.

<a name="ENREF_114">\114</a>      Nguyen, P.-M.; Lyathaud, C.; Vitrac, O. A two-scale pursuit method for the tailored identification and quantification of unknown polymer additives and contaminants by 1H NMR. *Industrial & Engineering Chemistry Research* **2015,** 54 (10), 2667-2681.

<a name="ENREF_115">\115</a>      Masood, S. H.; KeshavaMurthy, V. Development of collapsible PET water fountain bottles. *Journal of Materials Processing Technology* **2005,** 162-163, 83-89.

<a name="ENREF_116">\116</a>      Huang, H.-H.; Chen, L.-W.; Lu, W.-H.; Lin, W.-C.; Chen, Y.-C. Design and Simulation Analysis of Lightweight HDPE Milk Bottle. *Polymers and Polymer Composites* **2018,** 26 (1), 91-98.

<a name="ENREF_117">\117</a>      Hu, Q. C.; Sha, W. J.; Li, Y. H.; Wang, Y. S. Structural Optimization and Lightweight Design of PET Bottle Based on ABAQUS. *Advanced Materials Research* **2012,** 346, 558-563.

<a name="ENREF_118">\118</a>      Demirel, B.; Daver, F. The effects on the properties of PET bottles of changes to bottle-base geometry. *Journal of Applied Polymer Science* **2009,** 114 (6), 3811-3818.

<a name="ENREF_119">\119</a>      Roe, R.-J.; Bair, H. E.; Gieniewski, C. Solubility and diffusion coefficient of antioxidants in polyethylene. *Journal of Applied Polymer Science* **1974,** 18 (3), 843-856.

<a name="ENREF_120">\120</a>      Moisan, J. Y. Diffusion des additifs du polyethylene—I: Influence de la nature du diffusant. *European Polymer Journal* **1980,** 16 (10), 979-987.

<a name="ENREF_121">\121</a>      Ju, S. T.; Duda, J. L.; Vrentas, J. S. Influence of temperature on the diffusion of solvents in polymers above the glass transition temperature. *Industrial & Engineering Chemistry Product Research and Development* **1981,** 20 (2), 330-335.

<a name="ENREF_122">\122</a>      Ehlich, D.; Sillescu, H. Tracer diffusion at the glass transition. *Macromolecules* **1990,** 23 (6), 1600-1610.

<a name="ENREF_123">\123</a>      Mauritz, K. A.; Storey, R. F.; George, S. E. A general free volume-based theory for the diffusion of large molecules in amorphous polymers above the glass temperature. I. Application to di-n-alkyl phthalates in PVC. *Macromolecules* **1990,** 23 (2), 441-450.

<a name="ENREF_124">\124</a>      Arnould, D.; Laurence, R. L. Size effects on solvent diffusion in polymers. *Industrial & Engineering Chemistry Research* **1992,** 31 (1), 218-228.

<a name="ENREF_125">\125</a>      Griffiths, M. C.; Strauch, J.; Monteiro, M. J.; Gilbert, R. G. Measurement of Diffusion Coefficients of Oligomeric Penetrants in Rubbery Polymer Matrixes. *Macromolecules* **1998,** 31 (22), 7835-7844.

<a name="ENREF_126">\126</a>      Hall, D. B.; Hamilton, K. E.; Miller, R. D.; Torkelson, J. M. Translational and Rotational Diffusion of Probe Molecules in Polymer Films near Tg:  Effect of Hydrogen Bonding. *Macromolecules* **1999,** 32 (24), 8052-8058.

<a name="ENREF_127">\127</a>      Vitrac, O.; Hayert, M. Identification of Diffusion Transport Properties from Desorption/Sorption Kinetics:  An Analysis Based on a New Approximation of Fick Equation during Solid−Liquid Contact. *Industrial & Engineering Chemistry Research* **2006,** 45 (23), 7941-7956.

<a name="ENREF_128">\128</a>      von Meerwall, E. D.; Lin, H.; Mattice, W. L. Trace Diffusion of Alkanes in Polyethylene:  Spin-Echo Experiment and Monte Carlo Simulation. *Macromolecules* **2007,** 40 (6), 2002-2007.

<a name="ENREF_129">\129</a>      Grabowski, C. A.; Mukhopadhyay, A. Size Effect of Nanoparticle Diffusion in a Polymer Melt. *Macromolecules* **2014,** 47 (20), 7238-7242.

<a name="ENREF_130">\130</a>      Vagias, A.; Schultze, J.; Doroshenko, M.; Koynov, K.; Butt, H.-J.; Gauthier, M.; Fytas, G.; Vlassopoulos, D. Molecular Tracer Diffusion in Nondilute Polymer Solutions: Universal Master Curve and Glass Transition Effects. *Macromolecules* **2015,** 48 (24), 8907-8912.

<a name="ENREF_131">\131</a>      Begley, T.; Castle, L.; Feigenbaum, A.; Franz, R.; Hinrichs, K.; Lickly, T.; Mercea, P.; Milana, M.; O’Brien, A.; Rebre, S.; Rijk, R.; Piringer, O. Evaluation of migration models that might be used in support of regulations for food-contact plastics AU - Begley, T. *Food Additives & Contaminants* **2005,** 22 (1), 73-90.

<a name="ENREF_132">\132</a>      Hong, S.-U. Predicting ability of free-volume theory for solvent self-diffusion coefficients in rubbers. *Journal of Applied Polymer Science* **1996,** 61 (5), 833-841.

<a name="ENREF_133">\133</a>      Vrentas, J. S.; Vrentas, C. M. Predictive methods for self-diffusion and mutual diffusion coefficients in polymer–solvent systems. *European Polymer Journal* **1998,** 34 (5–6), 797-803.

<a name="ENREF_134">\134</a>      Vrentas, J. S.; Vrentas, C. M.; Faridi, N. Effect of Solvent Size on Solvent Self-Diffusion in Polymer−Solvent Systems. *Macromolecules* **1996,** 29 (9), 3272-3276.

<a name="ENREF_135">\135</a>      Zielinski, J. M.; Duda, J. L. Predicting polymer/solvent diffusion coefficients using free-volume theory. *AIChE Journal* **1992,** 38 (3), 405-415.

<a name="ENREF_136">\136</a>      Fang, X.; Vitrac, O. Predicting diffusion coefficients of chemicals in and through packaging materials. *Critical Reviews in Food Science and Nutrition* **2017,** 57 (2), 275-312.

<a name="ENREF_137">\137</a>      Karlsson, G. E.; Johansson, T. S.; Gedde, U. W.; Hedenqvist, M. S. Determining Limonene Diffusion in Molten Polyethylene from within 0.1 μs Molecular Dynamics Trajectories. *Macromolecules* **2002,** 35 (19), 7453-7459.

<a name="ENREF_138">\138</a>      Vitrac, O.; Lézervant, J.; Feigenbaum, A. Decision trees as applied to the robust estimation of diffusion coefficients in polyolefins. *Journal of Applied Polymer Science* **2006,** 101 (4), 2167-2186.

<a name="ENREF_139">\139</a>      Harmandaris, V. A.; Adhikari, N. P.; van der Vegt, N. F. A.; Kremer, K.; Mann, B. A.; Voelkel, R.; Weiss, H.; Liew, C. Ethylbenzene Diffusion in Polystyrene:  United Atom Atomistic/Coarse Grained Simulations and Experiments. *Macromolecules* **2007,** 40 (19), 7026-7035.

<a name="ENREF_140">\140</a>      Durand, M.; Meyer, H.; Benzerara, O.; Baschnagel, J.; Vitrac, O. Molecular dynamics simulations of the chain dynamics in monodisperse oligomer melts and of the oligomer tracer diffusion in an entangled polymer matrix. *Journal of Chemical Physics* **2010,** 132 (19), 194902.

<a name="ENREF_141">\141</a>      Gautieri, A.; Vesentini, S.; Redaelli, A. How to predict diffusion of medium-sized molecules in polymer matrices. From atomistic to coarse grain simulations. *Journal of Molecular Modeling* **2010,** 16 (12), 1845-1851.

<a name="ENREF_142">\142</a>      Lin, E.; You, X.; Kriegel, R. M.; Moffitt, R. D.; Batra, R. C. Atomistic to coarse grained simulations of diffusion of small molecules into polymeric matrix. *Computational Materials Science* **2017,** 138, 448-461.

<a name="ENREF_143">\143</a>      Hoekstra, E.; Brandsch, R.; Dequatre, C.; Mercea, P.; Milana, M. R.; Stoermer, A.; Trier, X.; Vitrac, O.; Schaefer, A.; Simoneau, C. *Practical guidelines on the application of migration modelling for the estimation of specific migration - EUR 27529 EN -* [*https://ec.europa.eu/jrc/en/publication/practical-guidelines-application-migration-modelling-estimation-specific-migration*](https://ec.europa.eu/jrc/en/publication/practical-guidelines-application-migration-modelling-estimation-specific-migration); 2015; p 36.

<a name="ENREF_144">\144</a>      Kirkwood, J. G.; Buff, F. P. The Statistical Mechanical Theory of Solutions. I. *The Journal of Chemical Physics* **1951,** 19 (6), 774-777.

<a name="ENREF_145">\145</a>      Flory, P. J. Thermodynamics of High Polymer Solutions. *The Journal of Chemical Physics* **1942,** 10 (1), 51-61.

<a name="ENREF_146">\146</a>      Flory, P. J. The Configuration of Real Polymer Chains. *The Journal of Chemical Physics* **1949,** 17 (3), 303-310.

<a name="ENREF_147">\147</a>      Flory, P. J., *Principles of Polymer Chemistry*. Cornell University Press: Iathaca, NY, USA, 1953.

<a name="ENREF_148">\148</a>      Kadam, A.; Karbowiak, T.; Voilley, A.; Bellat, J.-P.; Vitrac, O.; Debeaufort, F. Sorption of n-hexane in amorphous polystyrene. *Journal of Polymer Science Part B: Polymer Physics* **2014,** 52 (19), 1252-1258.

<a name="ENREF_149">\149</a>      Hansen, C. M., *Hansen Solubility Parameters: A User's Handbook, Second Edition*. CRC Press: Boca-Raton, FL, USA, 2007; p 519.

<a name="ENREF_150">\150</a>      Baner, A. L.; Piringer, O. G. Prediction of solute partition coefficients between polyolefins and alcohols using the regular solution theory and group contribution methods. *Industrial & Engineering Chemistry Research* **1991,** 30 (7), 1506-1515.

<a name="ENREF_151">\151</a>      Kontogeorgis, G. M.; Folas, G. K., *Thermodynamic Models for Industrial Applications: From Classical and Advanced Mixing Rules to Association Theories*. John Wiley & Sons Ltd: Chichester, UK, 2009; p 692.

<a name="ENREF_152">\152</a>      De Angelis, M. G.; Boulougouris, G. C.; Theodorou, D. N. Prediction of Infinite Dilution Benzene Solubility in Linear Polyethylene Melts via the Direct Particle Deletion Method. *The Journal of Physical Chemistry B* **2010,** 114 (19), 6233-6246.

<a name="ENREF_153">\153</a>      Boulougouris, G. C. On the Estimation of the Free Energy, From a Single Equilibrium Statistical Ensemble, via Particle Reinsertion. *The Journal of Physical Chemistry B* **2011,** 116 (3), 997-1006.

<a name="ENREF_154">\154</a>      Boulougouris, G. C. Calculation of the Chemical Potential beyond the First-Order Free-Energy Perturbation: From Deletion to Reinsertion†. *Journal of Chemical & Engineering Data* **2010,** 55 (10), 4140-4146.

<a name="ENREF_155">\155</a>      Özal, T. A.; Peter, C.; Hess, B.; van der Vegt, N. F. A. Modeling Solubilities of Additives in Polymer Microstructures: Single-Step Perturbation Method Based on a Soft-Cavity Reference State. *Macromolecules* **2008,** 41 (13), 5055-5061.

<a name="ENREF_156">\156</a>      Hess, B.; van der Vegt, N. F. A. Predictive Modeling of Phenol Chemical Potentials in Molten Bisphenol A−Polycarbonate over a Broad Temperature Range. *Macromolecules* **2008,** 41 (20), 7281-7283.

<a name="ENREF_157">\157</a>      Hess, B.; Peter, C.; Ozal, T.; van der Vegt, N. F. A. Fast-Growth Thermodynamic Integration: Calculating Excess Chemical Potentials of Additive Molecules in Polymer Microstructures. *Macromolecules* **2008,** 41 (6), 2283-2289.

<a name="ENREF_158">\158</a>      Gillet, G.; Vitrac, O.; Desobry, S. Prediction of Partition Coefficients of Plastic Additives between Packaging Materials and Food Simulants. *Industrial & Engineering Chemistry Research* **2010,** 49 (16), 7263-7280.

<a name="ENREF_159">\159</a>      Vitrac, O.; Gillet, G. An Off-Lattice Flory-Huggins Approach of the Partitioning of Bulky Solutes between Polymers and Interacting Liquids. *International Journal of Chemical Reactor Engineering* **2010,** 8.

<a name="ENREF_160">\160</a>      Nguyen, P.-M.; Guiga, W.; Dkhissi, A.; Vitrac, O. Off-lattice Flory-Huggins approximations for the tailored calculation of activity coefficients of organic solutes in random and block copolymers. *Industrial & Engineering Chemistry Research* **2017,** 56 (3), 774–787.

<a name="ENREF_161">\161</a>      Nguyen, P.-M.; Guiga, W.; Vitrac, O. Molecular thermodynamics for food science and engineering. *Food Research International* **2017,** 88, Part A, 91-104.

<a name="ENREF_162">\162</a>      Tribble, A., *Fundamentals of Contamination Control*. SPIE Press: Washington, DC, USA, 2000; p 174.

<a name="ENREF_163">\163</a>      Lee, G. Materials for ultra-high vacuum. . *Illinois* **1989**.

<a name="ENREF_164">\164</a>      Weissler, G.; Carlson, R., *Vacuum Physics and Technology*. Academic Press: New York, 1979.

<a name="ENREF_165">\165</a>      Alexander Stern, S. Polymers for gas separations: the next decade. *Journal of Membrane Science* **1994,** 94 (1), 1-65.

<a name="ENREF_166">\166</a>      Klopffer, M. H.; Flaconneche, B. Transport properties of gases in polymers: Bibliographic review. *Oil Gas Sci. Technol.* **2001,** 56 (3), 223-244.

<a name="ENREF_167">\167</a>      Masaro, L.; Zhu, X. X. Physical models of diffusion for polymer solutions, gels and solids. *Progress in Polymer Science* **1999,** 24 (5), 731-775.

<a name="ENREF_168">\168</a>      Pigeonneau, F.; Jaffrennou, B.; Letailleur, A.; Limouzin, K. Numerical investigation of generalized Graetz problem in circular tube with a mass transfer coupling between the solid and the liquid. *International Journal of Heat and Mass Transfer* **2016,** 96, 381-395.

<a name="ENREF_169">\169</a>      Gao, S.; Wang, Z.-W.; Hu, C.-Y.; Wu, Y.-M. Investigation of Migration Model of Printing Inks on Paper Packaging. *Journal of Food Process Engineering* **2014,** 37 (2), 146-159.

<a name="ENREF_170">\170</a>      Aparicio, J. L.; Elizalde, M. Migration of Photoinitiators in Food Packaging: A Review. *Packaging Technology and Science* **2015,** 28 (3), 181-203.

<a name="ENREF_171">\171</a>      Salafranca, J.; Clemente, I.; Isella, F.; Nerín, C.; Bosetti, O. Influence of oxygen and long term storage on the profile of volatile compounds released from polymeric multilayer food contact materials sterilized by gamma irradiation. *Analytica Chimica Acta* **2015,** 878, 118-130.

<a name="ENREF_172">\172</a>      Szendi, K.; Gyöngyi, Z.; Kontár, Z.; Gerencsér, G.; Berényi, K.; Hanzel, A.; Fekete, J.; Kovács, A.; Varga, C. Mutagenicity and Phthalate Level of Bottled Water Under Different Storage Conditions. *Exposure and Health* **2018,** 10 (1), 51-60.

<a name="ENREF_173">\173</a>      Careghini, A.; Mastorgio, A. F.; Saponaro, S.; Sezenna, E. Bisphenol A, nonylphenols, benzophenones, and benzotriazoles in soils, groundwater, surface water, sediments, and food: a review. *Environmental Science and Pollution Research* **2015,** 22 (8), 5711-5741.

<a name="ENREF_174">\174</a>      Spack, L. W.; Leszczyk, G.; Varela, J.; Simian, H.; Gude, T.; Stadler, R. H. Understanding the contamination of food with mineral oil: the need for a confirmatory analytical and procedural approach. *Food Additives & Contaminants: Part A* **2017,** 34 (6), 1052-1071.

<a name="ENREF_175">\175</a>      FDA. Food Additives; Threshold of Regulation for Substances Used in Food-Contact Articles. *Federal Register* **1995,** 60 (136), 36582-36596.

<a name="ENREF_176">\176</a>      Kroes, R.; Renwick, A. G.; Cheeseman, M.; Kleiner, J.; Mangelsdorf, I.; Piersma, A.; Schilter, B.; Schlatter, J.; van Schothorst, F.; Vos, J. G.; Würtzen, G. Structure-based thresholds of toxicological concern (TTC): guidance for application to substances present at low levels in the diet. *Food and Chemical Toxicology* **2004,** 42 (1), 65-83.

<a name="ENREF_177">\177</a>      Sagiv, A. Exact solution of mass diffusion into a finite volume. *Journal of Membrane Science* **2001,** 186 (2), 231-237.

<a name="ENREF_178">\178</a>      Sagiv, A. Theoretical formulation of the diffusion through a slab—theory validation. *Journal of Membrane Science* **2002,** 199 (1), 125-134.

<a name="ENREF_179">\179</a>      Goujot, D.; Vitrac, O. Extension to nonlinear adsorption isotherms of exact analytical solutions to mass diffusion problems. *Chemical Engineering Science* **2013,** 99, 2-22.

<a name="ENREF_180">\180</a>      Feigenbaum, A.; Dole, P.; Aucejo, S.; Dainelli, D.; Garcia, C. D. L. C.; Hankemeier, T.; N'Gono, Y.; Papaspyrides, C. D.; Paseiro, P.; Pastorelli, S.; Pavlidou, S.; Pennarun, P. Y.; Saillard, P.; Vidal, L.; Vitrac, O.; Voulzatis, Y. Functional barriers: Properties and evaluation AU - Feigenbaum, A. *Food Additives & Contaminants* **2005,** 22 (10), 956-967.

<a name="ENREF_181">\181</a>      Nguyen, P.-M.; Goujon, A.; Sauvegrain, P.; Vitrac, O. A Computer-Aided Methodology to Design Safe Food Packaging and Related Systems. *Aiche Journal* **2013,** 59 (4), 1183-1212.

<a name="ENREF_182">\182</a>      Vitrac, O.; Hayert, M. Effect of the distribution of sorption sites on transport diffusivities: A contribution to the transport of medium-weight-molecules in polymeric materials. *Chemical Engineering Science* **2007,** 62 (9), 2503-2521.

<a name="ENREF_183">\183</a>      Courgneau, C.; Vitrac, O.; Ducruet, V.; Riquet, A.-M. Local demixion in plasticized polylactide probed by electron spin resonance. *Journal of Magnetic Resonance* **2013,** 233, 37-48.

<a name="ENREF_184">\184</a>      Vrentas, J. S.; Vrentas, C. M. A new equation relating self-diffusion and mutual diffusion coefficients in polymer-solvent systems. *Macromolecules* **1993,** 26 (22), 6129-6131.

<a name="ENREF_185">\185</a>      Fang, X.; Domenek, S.; Ducruet, V.; Refregiers, M.; Vitrac, O. Diffusion of Aromatic Solutes in Aliphatic Polymers above Glass Transition Temperature. *Macromolecules* **2013,** 46 (3), 874-888.

<a name="ENREF_186">\186</a>      Ewender, J.; Welle, F. Determination of the activation energies of diffusion of organic molecules in poly(ethylene terephthalate). *Journal of Applied Polymer Science* **2013,** 128 (6), 3885-3892.

<a name="_ENREF_187">\187</a>.      Welle, F. A new method for the prediction of diffusion coefficients in poly(ethylene terephthalate). *Journal of Applied Polymer Science* **2013,** 129 (4), 1845-1851.

<a name="_ENREF_188">\188</a>.      Ewender, J.; Welle, F. Determination and Prediction of the Lag Times of Hydrocarbons through a Polyethylene Terephthalate Film. *Packaging Technology and Science* **2014,** 27 (12), 963-974.

<a name="_ENREF_189">\189</a>.      Ewender, J.; Welle, F. Determination of the activation energies of diffusion of organic molecules in poly(ethylene terephthalate). *Journal of Applied Polymer Science* **2018**.

<a name="_ENREF_190">\190</a>.      Zhao, Y. H.; Abraham, M. H.; Zissimos, A. M. Fast Calculation of van der Waals Volume as a Sum of Atomic and Bond Contributions and Its Application to Drug Compounds. *The Journal of Organic Chemistry* **2003,** 68 (19), 7368-7373.

<a name="_ENREF_191">\191</a>.      NIST NIST Chemistry WebBook. https://webbook.nist.gov/chemistry/ (accessed on March 31, 2019), 

<a name="_ENREF_192">\192</a>.      Hinrichs, K.; Piringer, O. Evaluation of Migration Models to Be Used Under Directive 90/128/EEC; Final repot Contract SMT4-CT9867513. EUR 20604 EN. *European Commission, Directorate General for Research* **2002**.

<a name="_ENREF_193">\193</a>.      Flynn, J. H. A collection of kinetic data for the diffusion of organic compounds in polyolefins. *Polymer* **1982,** 23 (9), 1325-1344.

<a name="_ENREF_194">\194</a>.      Ewender, J.; Welle, F. Functional Barrier Performance of a Polyamide-6 Membrane Towards n-Alkanes and 1-Alcohols. *Packaging Technology and Science* **2016,** 29 (6), 277-287.

<a name="_ENREF_195">\195</a>.      Boyne, J. A.; Williamson, A. G. Enthalpies of mixture of ethanol and water at 25.degree.C. *Journal of Chemical & Engineering Data* **1967,** 12 (3), 318.

<a name="_ENREF_196">\196</a>.      Krüger, K.-M.; Sadowski, G. Fickian and Non-Fickian Sorption Kinetics of Toluene in Glassy Polystyrene. *Macromolecules* **2005,** 38 (20), 8408-8417.

<a name="_ENREF_197">\197</a>.      Vitrac, O.; Hayert, M. Risk assessment of migration from packaging materials into foodstuffs. *AIChE Journal* **2005,** 51 (4), 1080-1095.

<a name="_ENREF_198">\198</a>.      Vitrac, O.; Mougharbel, A.; Feigenbaum, A. Interfacial mass transport properties which control the migration of packaging constituents into foodstuffs. *Journal of Food Engineering* **2007,** 79 (3), 1048-1064.

<a name="_ENREF_199">\199</a>.      Till, D.; Schwope, A. D.; Ehntholt, D. J.; Sidman, K. R.; Whelan, R. H.; Schwartz, P. S.; Reid, R. C.; Rainey, M. L. Indirect Food Additive Migration from Polymeric Food Packaging Materials. *CRC Critical Reviews in Toxicology* **1987,** 18 (3), 215-243.

<a name="_ENREF_200">\200</a>.      Zhu, Y.; Guillemat, B.; Vitrac, O. Rational design of packaging: towards safer and eco-designed food packaging systems. *accepted in Frontiers in Chemistry. Chemical Engineering.* **2019**.

<a name="_ENREF_201">\201</a>.      EFSA; WHO. Review of the Threshold of Toxicological Concern (TTC) approach and development of new TTC decision tree. *EFSA Supporting Publications* **2016,** 13 (3), 1006E.

<a name="_ENREF_202">\202</a>.      Cramer, G. M.; Ford, R. A.; Hall, R. L. Estimation of toxic hazard—A decision tree approach. *Food and Cosmetics Toxicology* **1976,** 16 (3), 255-276.

<a name="_ENREF_203">\203</a>.      CEN *Materials modelling - Terminology, classification and metadata - CWA 17284 -* [*https://www.cen.eu/news/workshops/Pages/WS-2017-012.aspx*](https://www.cen.eu/news/workshops/Pages/WS-2017-012.aspx); EUROPEAN COMMITTEE FOR STANDARDIZATION: Brussels, Belgium, 2018; pp 1-18.

 