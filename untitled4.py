# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 10:08:05 2015
Huang Parameter File
@author: Samuel
"""

import numpy as np

Avogadro = 6.0221413e23 # Particles per mole
Kb = 0.0019872041 # Kcal/mol/K
T = 300.0 # Kelvin
Beta = 1.0/(T*Kb)

ChainLength = 150
NumChains = 300
NumPCBM = 0
NumTotal = NumChains*ChainLength + NumPCBM

Density = .005 # g/mol

Total_Moles = NumTotal/ Avogadro

Monomer_MW = 166.3

Molar_Volume = Monomer_MW/ Density

Volume = Molar_Volume*Total_Moles

Box_Length_CM = Volume**(1./3.)
Box_Length = Box_Length_CM*100000000