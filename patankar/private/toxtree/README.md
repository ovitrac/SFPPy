# Toxtree for SFPPy 🍏⏩🍎

This directory contains the private copy of **Toxtree** used by SFPPy. Toxtree is not supplied initially and you need to install it by yourself in SFPPy installation folder.



> 1) Download the last version of Toxtree from https://sourceforge.net/projects/toxtree/.
> 2) choose the Zip archive
> 3) Unzip it in `patankar/private/toxtree`
> 4) Download the 18 plugins in `patankar/private/toxtree/ext`
> 5) Edit `toxtree-plugins.properties`based on available plugins.



## 📁 Files and Subdirectories
- `ext/` 📁: Plugins folder
- `toxtree-plugins.properties` 📏: Plugin configuration file
- `Toxtree-3.1.0.1851.jar` 🔬: the main executable of Toxtree.

## 🔹 Example of `toxtree-plugins.properties`
```bash
#Folder to look for plugins
folder=ext
#Plugin jars; set true to load; false to skip
01=toxTree.tree.cramer.CramerRules
02=toxtree.tree.cramer3.RevisedCramerDecisionTree
03=toxtree.plugins.kroes.Kroes1Tree
04=verhaar.VerhaarScheme
05=toxtree.plugins.verhaar2.VerhaarScheme2
06=mutant.BB_CarcMutRules
07=toxtree.plugins.ames.AmesMutagenicityRules
08=sicret.SicretRules
09=eye.EyeIrritationRules
10=toxtree.plugins.skinsensitisation.SkinSensitisationPlugin
#11=michaelacceptors.MichaelAcceptorRules
12=com.molecularnetworks.start.BiodgeradationRules
13=toxtree.plugins.smartcyp.SMARTCYPPlugin
14=mic.MICRules
15=toxtree.plugins.func.FuncRules
16=toxtree.plugins.proteinbinding.ProteinBindingPlugin
17=toxtree.plugins.dnabinding.DNABindingPlugin
#999=toxTree.tree.demo.SubstructureTree
#18=toxtree.tree.cramer3.RevisedCramerDecisionTree
20=cramer2.CramerRulesWithExtensions
```

***

<div style="border: 2px solid #4CAF50; border-radius: 8px; padding: 10px; background: linear-gradient(to right, #4CAF50, #FF4D4D); color: white; text-align: center; font-weight: bold;">
  <span style="font-size: 20px;">🍏⏩🍎 <strong>SFPPy for Food Contact Compliance and Risk Assessment</strong></span><br>
  Contact <a href="mailto:olivier.vitrac@gmail.com" style="color: #fff; text-decoration: underline;">Olivier Vitrac</a> for questions |
  <a href="https://github.com/ovitrac/SFPPy" style="color: #fff; text-decoration: underline;">Website</a> |
  <a href="https://ovitrac.github.io/SFPPy/" style="color: #fff; text-decoration: underline;">Documentation</a>
</div>

