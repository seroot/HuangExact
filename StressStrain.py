# -*- coding: utf-8 -*-
"""
Created on Mon Apr 06 10:43:23 2015
Functions for analyzing thermodynamic output of lammps simulations
@author: Samuel
"""

import numpy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os


def Poly10(x,a,b,c,d,e,f,h,g,i,k):
	return a*x**10 + b*x**9 + c*x**8 + d*x**7 + e*x**6 + f*x**5 +h*x**4 + g*x**3 + i*x**2 + k*x
 
def Poly10Val(x,Pop):
	return Pop[0]*x**10 + Pop[1]*x**9 + Pop[2]*x**8 + Pop[3]*x**7 + Pop[4]*x**6 + Pop[5]*x**5 +Pop[6]*x**4 + Pop[7]*x**3 + Pop[8]*x**2 + Pop[9]*x
 
 
def Moving_Average(interval, window_size):
	window = numpy.ones(int(window_size))/float(window_size)
	return numpy.convolve(interval, window, 'same')
 
def quadratic(x, m, b, c):
    return m*x**2 + b*x + c

def Linear(x,m):
	return m*x
	
 

def Plot_Stress_Strain( Filename ):
    """ This function is for plotting and analyzing the stress strain curves outputed from a tensile deformation simulation run in lammps. 
    
    Input:
        Filename = name of stress-strain file to plot, ( strainx strainy strainz stressx )
    
    """
    
    
    File = open(Filename)
    Strainx = []
    Strainy = []
    Strainz = []
    Stress = []
    
    for line in File:
        line = line.split()
        try:
            if (float(line[0]>0.0)):
                Strainx.append(float(line[0]))
                Strainy.append(float(line[1]))
                Strainz.append(float(line[2]))
                Stress.append(float(line[3]))
        except:
            continue
        
    
    Strainx = numpy.asarray(Strainx)
    Strainy = numpy.asarray(Strainy)
    Strainz = numpy.asarray(Strainz)
    Stress = numpy.asarray(Stress)
    window = 1000
    StressAVG = Moving_Average(Stress, window)
    Pop, Pco = curve_fit(Poly10, Strainx, Stress)
    Popt, Pcov = curve_fit(Linear, Strainx, Strainy)
    
    plt.figure(1)
    plt.subplot(211)
    plt.xlim((0,3.0))
    plt.ylim((0,.22))
    plt.plot(Strainx, StressAVG,'g:', label = 'Moving Average Window = %d' % window)
    plt.plot(Strainx, Poly10Val(Strainx, Pop),'r', linewidth = 4, label = '10th Order Polynomial fit')
    plt.plot(Strainx, Pop[-1]*Strainx, linewidth = 4, label = 'Linear Regime E = %.2f GPa ' % Pop[-1])
    plt.ylabel('Stress (GPa)', fontsize = 20)
    plt.title('Uniaxial Tensile Deformation', fontsize = 30)
    plt.legend(loc= 'lower right')
    
    plt.subplot(212)
    plt.plot(Strainx, Strainy, label = 'Y Strain')
    plt.plot(Strainx, -.37*Strainx, linewidth = 4,label = 'Linear Regime v = %.2f ' % .37)
    plt.ylabel( 'Transverse Strain', fontsize = 20)
    plt.xlabel('Axial Strain', fontsize = 30)
    plt.xlim((0,3.0))
    plt.ylim((0,-.4))
    plt.legend()
    plt.show()
    
    return
    

os.chdir('P3HTHuang150_100')

Plot_Stress_Strain('stresstrainP3HTHuang.txt')

        

    
    