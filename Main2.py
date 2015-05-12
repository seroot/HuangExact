import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import PolyModelFunctionsHuang as poly
import Huangparameters as Huang
import HuangData

def main():
	PDF, CDF = poly.Gen_PDF_CDF_Multi(Huang.P1P1P1P1, Huang.Beta)
	#Position_M, Position_S1, Position_S2, Bases = poly.generatepolymer( 100, Huang.BondLengths, Huang.Angles, CDF, Huang.Box_Length)
	
	Position_M, Position_S1, Position_S2, Bases = poly.Gen_Many_Polymers( 100, 10, Huang.BondLengths, Huang.Angles, CDF, Huang.SigmaM_M, 1000.)
	poly.Scatter_Plot(Position_M, Position_S1, Position_S2, Huang.Box_Length)
	return

if __name__=='__main__': main()

