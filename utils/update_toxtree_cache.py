#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This small script update ToxTree cache for getting toxicological data when Toxtree is not installed

Created on Fri Apr  4 11:25:13 2025

@author: olivi
"""

from patankar.private.EUFCMannex1 import EuFCMannex1
from patankar.private.USFDAfcn import USFDAfcn
from patankar.private.GBappendixA import GBappendixA

from patankar.loadpubchem import migrantToxtree
dbeu = EuFCMannex1()
dbus = USFDAfcn()
dbgb = GBappendixA()


# %% General refresh function
def refreshTox(CAS):
    CAS = r["CAS"]
    if isinstance(CAS,list):
        for c in CAS:
            try:
                m = migrantToxtree(c)
            except:
                print(f'no data for {c}')
    else:
        try:
            m = migrantToxtree(CAS)
        except:
            print(f'no data for {CAS}')



# %% 🇪🇺 10/2011/EC
for i,r in enumerate(dbeu):
    print(f"EU{i} cid={r.cid}")
    if r.cid:
        refreshTox(r["CAS"])


# %% 🇺🇸 FCN list
for i,r in enumerate(dbus):
    print(f"US{i} cid={r.cid}")
    if r.cid:
        refreshTox(r["CAS"])


# %%  🇨🇳 GB 9685-201
for i,r in enumerate(dbus):
    print(f"GB{i} cid={r.cid}")
    if r.cid:
        refreshTox(r["CAS"])
