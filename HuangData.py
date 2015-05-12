# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 16:44:34 2015
File containing a function for preparing a LAMMPS Data file to implement the huang model
@author: Samuel
"""

import sys
import numpy as np
from Huangparameters import *
import os
import matplotlib.pyplot as plt

def P3HTHuangDataLammps(Position_M, Position_S1, Position_S2, Filename):
    
    File = open(Filename, 'w')
    
    NumAtoms = NumChains*ChainLength*3
    NumBonds = NumChains*(ChainLength - 1 + ChainLength*2)
    NumAngles = NumChains*(ChainLength + 3*ChainLength - 4)
    NumDih = NumChains*(ChainLength - 3 + ChainLength - 1 + ChainLength - 1 + ChainLength - 1)
    #NumImproper = NumChains*(ChainLength - 2)
    NumImproper = 0
    
    File.write('Created by Sam Root\n\n')   
    File.write('\t%d\tatoms\n' % NumAtoms)
    File.write('\t%d\tbonds\n' % NumBonds)
    File.write('\t%d\tangles\n' % NumAngles)
    File.write('\t%d\tdihedrals\n'% NumDih)
    File.write('\t%d\timpropers\n\n' % NumImproper)
    
    File.write('\t3\tatom types\n')
    File.write('\t3\tbond types\n')
    File.write('\t4\tangle types\n')
    File.write('\t4\tdihedral types\n')
    File.write('\t0\timproper types\n\n')
    
    File.write('\t%.6f\t%.6f xlo xhi\n' %(0., Box_Length))
    File.write('\t%.6f\t%.6f ylo yhi\n' %(0., Box_Length))
    File.write('\t%.6f\t%.6f zlo zhi\n\n' %(0., Box_Length))

    # Declare Masses
    File.write('Masses\n\n')
    File.write('\t1 %.4f\n' % MP1)
    File.write('\t2 %.4f\n' % MP2)
    File.write('\t3 %.4f\n\n' % MP3)
    
    # Declare Bond Coefficients 
    File.write('Bond Coeffs\n\n')
    File.write('\t1\t%.4f\t%.4f\n' %(P1P1[1], P1P1[0]))
    File.write('\t2\t%.4f\t%.4f\n' %(P1P2[1], P1P2[0]))
    File.write('\t3\t%.4f\t%.4f\n' %(P2P3[1], P2P3[0]))
    
    # Declare Angle Coefficients 
    File.write('Angle Coeffs\n\n')
    File.write('\t1\t%.4f\t%.4f\n' %( 542.08, 171.38) )
    File.write('\t2\t%.4f\t%.4f\n' %( 448.25, 157.74) )
    File.write('\t3\t%.4f\t%.4f\n' %( 505.58, 110.08) )
    File.write('\t4\t%.4f\t%.4f\n\n' %( 10., P2P1P1[0]) )
    
    """
    # Declare BondAngle Coeffs (Style = class2)
    File.write('BondAngle Coeffs\n\n')
    File.write('1 0.0 0.0 0.0 0.0\n')
    File.write('2 0.0 0.0 0.0 0.0\n')
    File.write('3 0.0 0.0 0.0 0.0\n')
    File.write('4 0.0 0.0 0.0 0.0\n\n')
    
    # Declare BondBond Coeff (Style = class2)
    File.write('BondBond Coeffs\n\n')
    File.write('1 0.0 0.0 0.0\n')
    File.write('2 0.0 0.0 0.0\n')
    File.write('3 0.0 0.0 0.0\n')
    File.write('4 0.0 0.0 0.0\n\n')
    """
    
    # Declare Dihedral Coeffs
    File.write('Dihedral Coeffs\n\n')
    File.write('\t1\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n' % (P1P1P1P1[0], P1P1P1P1[1], P1P1P1P1[2], P1P1P1P1[3], P1P1P1P1[4]))
    File.write('\t2\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n' % (P2P1P1P2OPT[0], P2P1P1P2OPT[1], P2P1P1P2OPT[2], P2P1P1P2OPT[3], P2P1P1P2OPT[4]))
    File.write('\t3\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n' % (P1P1P2P3[0], P1P1P2P3[1], P1P1P2P3[2], P1P1P2P3[3], P1P1P2P3[4]))
    File.write('\t4\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n\n' % (P3P2P1P1[0], P3P2P1P1[1], P3P2P1P1[2], P3P2P1P1[3], P3P2P1P1[4]))
    
    
     # Declare atoms initial conditions
    File.write('Atoms\n\n')
    Atom_id = 0
    Mol_id = 0
    k = 0
    for i in range(NumChains):
        Mol_id += 1
        for j in range(ChainLength):
            Atom_id += 1
            File.write('%d %d 1 %.2f %.2f %.2f\n' % ( Atom_id, Mol_id, Position_M[k,0], Position_M[k,1], Position_M[k,2]))
            Atom_id += 1
            File.write('%d %d 2 %.2f %.2f %.2f\n' % ( Atom_id, Mol_id,  Position_S1[k,0], Position_S1[k,1], Position_S1[k,2]))
            Atom_id += 1
            File.write('%d %d 3 %.2f %.2f %.2f\n' % ( Atom_id, Mol_id,  Position_S2[k,0], Position_S2[k,1], Position_S2[k,2]))
            k += 1
            
    # Declare Bonding Topology
    File.write('\n\nBonds\n\n')
    E = 0
    B =1 
    for i in range(NumChains):
        for j in range(ChainLength - 1):
            E+= 1
            File.write('%d 1 %d %d\n' % (E, B, B+3))
            E += 1
            File.write('%d 2 %d %d\n' % (E, B, B+1))
            E += 1
            File.write('%d 3 %d %d\n' % (E, B+1, B+2))
            B += 3
        E += 1
        File.write('%d 2 %d %d\n' % (E, B, B+1) )
        E += 1
        File.write('%d 3 %d %d\n' % (E, B+1, B+2))
        B+= 3
        
    
    # Declare Angle Topology
    File.write('\n\nAngles\n\n')
    E = 0
    B = 1
    for i in range(NumChains):
        for j in range(ChainLength - 2):
            E += 1
            File.write('%d 1 %d %d %d\n' % (E, B, B+3, B+6))
            E +=1 
            File.write('%d 2 %d %d %d\n' % (E, B, B+1, B+2))
            E +=1
            File.write('%d 3 %d %d %d\n' % (E, B, B+3, B+4))
            E += 1
            File.write('%d 4 %d %d %d\n' % (E, B+1, B, B+3))
            B+=3
        E += 1
        File.write('%d 2 %d %d %d\n' % (E, B, B+1, B+2))
        E += 1 
        File.write('%d 3 %d %d %d\n' % (E, B, B+3, B+4))
        E += 1
        File.write('%d 4 %d %d %d\n' % (E, B+1, B, B+3))
        E += 1
        File.write('%d 2 %d %d %d\n' % (E, B+3, B+4, B+5))
        B+= 6
        
        
        
    # Declare Dihedral Topology
    File.write('\n\nDihedrals\n\n')
    E = 0
    B = 1
    for i in range(NumChains):
        for j in range(ChainLength - 3):
            E += 1 
            File.write('%d 1 %d %d %d %d\n' % (E, B, B+3, B+6, B+9) )
            E += 1 
            File.write('%d 2 %d %d %d %d\n' % (E, B+1, B, B+3, B+4))
            E += 1
            File.write('%d 3 %d %d %d %d\n' % (E, B, B+3, B+4, B+5))
            E += 1
            File.write('%d 4 %d %d %d %d\n' % (E, B+2, B+1, B, B+3))
            B +=3
        E += 1
        File.write('%d 2 %d %d %d %d\n' % (E, B+1, B, B+3, B+4))
        E += 1
        File.write('%d 3 %d %d %d %d\n' % (E, B, B+3, B+4, B+5))
        E += 1
        File.write('%d 4 %d %d %d %d\n' % (E, B+2, B+1, B, B+3))
        B+= 3
        E+= 1
        File.write('%d 2 %d %d %d %d\n' % (E, B+1, B, B+3, B+4))
        E+= 1
        File.write('%d 3 %d %d %d %d\n' % (E, B, B+3, B+4, B+5))
        E += 1
        File.write('%d 4 %d %d %d %d\n' % (E, B+2, B+1, B, B+3))
        B += 6
        
    File.close()
    
    return


def Run_Huang_Equil(In_File, Data_File, Pair_Table, Restart_Out, Dump = 'None', Num_Steps= 10000, Time_Step = 1, Temp = 423.0):
    """
    Function for generating an input file for an NVT equilibration simulation of Huang Model
    
    Input Arguments:
        Required:
        In_File = name of the input file to be created
        Restart_Out = name of restart file to be outputed
        Data_File | Restart_In = File for declaring initial position and topology (atleast on of these is required)
        Pair_table = name of pair table file
    Optional Arguments:
        Dump = name of dumpfile to be outputed, Default is none
        Num_Steps = number of timesteps to integrate equations of motion
        Time_Step = size of time increment
        Temp = Temperature to run simulation at
    """
    File = open(In_File, 'w')
    File.write( '# Input file for running an NVT Equilibration in LAMMPS, Filename = %s\n\n' % In_File)
    File.write('units real\n')
    File.write('atom_style molecular\n')
    File.write('boundary p p p\n')
    File.write('bond_style harmonic\n')
    File.write('pair_style table linear 187 \n')
    File.write('angle_style harmonic\n')
    File.write('dihedral_style multi/harmonic\n')
    File.write('special_bonds lj 0 0 0\n')
    File.write('improper_style none\n')
    File.write('kspace_style none\n')
    File.write('read_data %s\n' % Data_File)
    File.write('pair_coeff 1 1 HuangPair2.table P1P1\n')
    File.write('pair_coeff 1 2 HuangPair2.table P1P2\n')
    File.write('pair_coeff 1 3 HuangPair2.table P1P3\n')
    File.write('pair_coeff 2 2 HuangPair2.table P2P2\n')
    File.write('pair_coeff 2 3 HuangPair2.table P2P3\n')
    File.write('pair_coeff 3 3 HuangPair2.table P3P3\n')
    if (Dump != 'none'):
        File.write('dump 1 all custom 200 %s id type mol xs ys zs vx vy vz\n' % Dump)
    File.write('neighbor 10.0 bin\n')
    File.write('neigh_modify every 10 delay 0 one 1000\n')
    File.write('velocity all create %f 1223\n' % Temp)
    File.write('thermo_style custom step temp press etotal density\n')
    File.write('fix 1 all nve/limit 0.5\n')
    File.write('fix 2 all temp/berendsen %f %f 15\n' % (Temp, Temp));  File.write('fix 3 all deform 1 x scale .2 y scale .2 z scale .2 remap x \n')
    File.write('timestep %d\n' % Time_Step)
    File.write('thermo 100\n')
    File.write('run %d\n' % Num_Steps)
    File.write('unfix 1\n')
    File.write('unfix 2\n')
    File.write('write_restart %s\n' % Restart_Out)
    File.close()
    
    os.system('mpiexec -np 8 lmp_mpi -in %s' % In_File )
    
    return
    
    
def EquilibrateNPT( In_File, Restart_In, Restart_Out, Dump = 'none', Num_Steps = 100000, Time_Step = 1, Temp = 423, Press = 1):
    """ Function for equilibrating an MD Simulation in an NPT Ensemble

    Required Inputs:
            In_File = Name of simulation input file
            
            Restart_In = Name of Restart file to read in
            
            Restart_Out = Name of Restart file to write out
            
            
    Optional Inputs:
            Dump = name of dump file to output, Default = none
            
            Num_Steps = number of steps to integrate over
            
            Time_Step = time increment
            
            Temp = Simulation Temperature
            
            Press = Simulation Pressure
    """
    File = open( In_File, 'w')
    
    File.write('# Input file for running an NPT Equilibration with LAMMPS, File name = %s\n\n' % In_File)
    File.write('read_restart %s\n' % Restart_In)
    File.write('pair_style table linear 187 \n');File.write('angle_style table linear 200\n')
    File.write('bond_style table linear 200\n');File.write('special_bonds lj 0 0 0\n')
    File.write('pair_coeff 1 1 HuangPair2.table P1P1\n')
    File.write('pair_coeff 1 2 HuangPair2.table P1P2\n')
    File.write('pair_coeff 1 3 HuangPair2.table P1P3\n')
    File.write('pair_coeff 2 2 HuangPair2.table P2P2\n')
    File.write('pair_coeff 2 3 HuangPair2.table P2P3\n')
    File.write('pair_coeff 3 3 HuangPair2.table P3P3\n');File.write('bond_coeff 1 HuangBond.table Bond1\n' )
    File.write('bond_coeff 2 HuangBond.table Bond2\n' )
    File.write('bond_coeff 3 HuangBond.table Bond3\n' )
    File.write('angle_coeff 1 HuangAngle.table Angle1\n')
    File.write('angle_coeff 2 HuangAngle.table Angle2\n')
    File.write('angle_coeff 3 HuangAngle.table Angle3\n')
    File.write('angle_coeff 4 HuangAngle.table Angle4\n')
    File.write('log log.NPT\n')
    
    File.write('fix 1 all npt temp %d %d 100 iso 0.0 %d 1000 drag 2\n' % (Temp, Temp, Press))
    File.write('fix 2 all momentum 1 linear 1 1 1\n')
    File.write('thermo_style custom step temp press etotal density\n')
    
    if (Dump != 'none'):
        File.write('dump 1 all custom 200 %s id type mol xs ys zs vx vy vz\n' % Dump)
    
    File.write('thermo 200\n')
    File.write('timestep %d\n' % Time_Step)
    File.write('run %d\n' % Num_Steps)
    File.write('unfix 1\n')
    File.write('unfix 2\n')
    File.write('write_restart %s\n' % Restart_Out )
    File.close()
    
    os.system('mpiexec -np 8 lmp_mpi -in %s' % In_File)
    
    return

def Extract_Trajectory(Filename, N, T, K):
    """
    Function to extract the trajectories from a sorted lammpstrj file
    inputs: Filename
            N = number of particles of each type
            T = number of timesteps
            k = # of types of particles
            
    outputs numpy arrays
    """
    
    File = open(Filename)
    Positions = np.zeros((K, N*T, 3))
    i= 0
    j= 0
    k = 0
    for line in File:
        line = line.split()
        if len(line) == 9:
            if line[1] == '1':
                Positions[0, i] = [float(line[3]), float(line[4]), float(line[5])]
                i += 1
            if line[1] == '2':
                Positions[1, j] = [float(line[3]), float(line[4]), float(line[5])]
                j+=1
            if line[1] == '3':
                Positions[2, k] = [float(line[3]), float(line[4]), float(line[5])]
                k+=1
    return Positions
    
def Order_Parameter(Positions, N, T, k, chainlength, avgT):
	Num_Chains = N/chainlength
	Num_angles = (chainlength-2)*Num_Chains
	Orientationx = np.zeros((Num_angles, T))
	for i in range(T):
		a=0
		b=0
		for j in range(Num_Chains):
			for k in range(chainlength - 2):
				rjk= (Positions[0, a*i] - Positions[0, a*i+2])
				rjk = rjk/np.linalg.norm(rjk)
				Orientationx[b, i] = 1.5*(rjk[0]**2) - .5
				a += 1
				b += 1
			a += 2
			
			
			
	if avgT:
		OrientationxAvgT = np.mean(Orientationx,axis=1)
		return OrientationxAvgT
	else:
		return Orientationx
    
    
def Persistence_Length( Positions, N, T, k, chainlength):
    Num_Chains = N/chainlength
    Num_bonds = (chainlength -1)*Num_Chains
    Orientation = np.zeros((T, Num_Chains, Num_Bonds))
    for i in range(T):
        a =0
        for j in range(Num_Chains):
            U0 = Positions[0, a*i] - Positions[0,a*i + 1]
            U0 = np.linalg.norm(U0)
            for k in range(Num_bonds):
                Uk = Position[0, a*i] - Positions[0, a*i  + 1]
                Uk = np.linalg.norm(Uk)
                Orientation[i, j, k] = U0[0]*Uk[0] + U0[1]*Uk[1] + U0[2]*Uk[2]
                a += 1
            a+=2
    
    OrientationAVGT = np.mean(Orientation, axis = 0)
    OrientationAVGChain = np.mean(OrientationAVGT, axis = 0)
    
    return Orientation, OrientationAVGT, OrientationAVGChain
    
                
                
                
    
def Calculate_LP( In_FILE, Restart_In, Restart_out, Num_Steps=10000, Time_Step = 1, Temp= 300, Press = 1):
    """
    Function for calculating the persistence length of Huang's model of P3HT
    """
    os.mkdir('Persistence Length')
    os.chdir('Persistence Length')
    
    Traj_File = 'Persistance_Length.lammpstrj'
    
    File = open(In_FILE, 'w')
    File.write('# Input file for calculating persistence length of P3HT in melt\n\n')
    File.write('read_restart %s\n' % Restart_In)
    File.write('pair_style table linear 187 \n')
    File.write('pair_coeff 1 1 HuangPair2.table P1P1\n')
    File.write('pair_coeff 1 2 HuangPair2.table P1P2\n')
    File.write('pair_coeff 1 3 HuangPair2.table P1P3\n')
    File.write('pair_coeff 2 2 HuangPair2.table P2P2\n')
    File.write('pair_coeff 2 3 HuangPair2.table P2P3\n')
    File.write('pair_coeff 3 3 HuangPair2.table P3P3\n')
    File.write('fix 1 all npt temp %d %d 100 iso %d %d 1000 drag 2\n' % (Temp, Temp, Press, Press))
    File.write('fix 2 all momentum 1 linear 1 1 1\n')
    File.write('thermo_style custom step temp press etotal density\n')
    File.write('thermo 500\n')
    File.write('timestep %d\n' % Time_Step)
    File.write('dump 1 all custom 100 %s id type mol xs ys zs\n' % Traj_File)
    File.write('dump_modify 1 sort id\n')
    File.write('run %d\n' % Num_Steps)
    File.write('write_restart %s\n' % Restart_out)
    
    File.close()
    os.system('mpiexec -np 8 lmp_mpi -in %s' % In_FILE)
    
    N = NumChains*ChainLength
    T = Num_Steps 
    
    
    Positions = Extract_Trajectory(Traj_File, N, 3)
    
    Corr, CorrT, CorrCh = Persistence_Length( Positions, N, T, 3, ChainLength)
    
    plt.figure()
    plt.plot(CorrCh)
    plt.show()    
    
    return
    
def Strain(In_File, Restart_In, Restart_Out, Stress_Strain_File , Dump = 'none', Strain_Rate = .000001, Num_Steps = 10000, Time_Step = 1, Press = 0, Temp = 423):
    """
    Function for running a uniaxial tensile deformation simulation with LAMMPS
    
    Required Inputs:
    In_File = name of simulation input file
    
    Restart_In = Name of restart file to read in
    
    Restart_Out = Name of restart file to write out
    
    Stress_Strain_File = Name of file to output stress/strain data (Standard Format)
    
    Optional Inputs:
    Strain_Rate = Rate of deformation 
    
    Num_Steps = Number of steps to integrate equations of motion
    
    Time_Step = Time incremant for integrations
    
    Press = externally applied transverse pressure
    
    Temp = Simulation Temperature
    
    """
    
    File = open( In_File, 'w')
    
    File.write('# Input file for running a uniaxial tensile deformation simulation, Filename = %s\n\n' % In_File)
    File.write('read_restart %s\n' % Restart_In)
    File.write('pair_style table linear 187 \n'); File.write('angle_style table linear 200\n')
    File.write('bond_style table linear 200\n')
    File.write('pair_coeff 1 1 HuangPair2.table P1P1\n')
    File.write('pair_coeff 1 2 HuangPair2.table P1P2\n')
    File.write('pair_coeff 1 3 HuangPair2.table P1P3\n')
    File.write('pair_coeff 2 2 HuangPair2.table P2P2\n')
    File.write('pair_coeff 2 3 HuangPair2.table P2P3\n')
    File.write('pair_coeff 3 3 HuangPair2.table P3P3\n');     File.write('bond_coeff 1 HuangBond.table Bond1\n' )
    File.write('bond_coeff 2 HuangBond.table Bond2\n' )
    File.write('bond_coeff 3 HuangBond.table Bond3\n' )
    File.write('angle_coeff 1 HuangAngle.table Angle1\n')
    File.write('angle_coeff 2 HuangAngle.table Angle2\n')
    File.write('angle_coeff 3 HuangAngle.table Angle3\n')
    File.write('angle_coeff 4 HuangAngle.table Angle4\n')
    File.write('variable tmp equal "lx"\n')
    File.write('variable L0 equal ${tmp}\n')
    File.write('variable strainx equal "(lx-v_L0)/v_L0"\n')
    File.write('variable strainy equal "(ly -v_L0)/v_L0"\n')
    File.write('variable strainz equal "(lz -v_L0)/v_L0"\n')
    File.write('variable stress equal "-pxx/10000*1.01325"\n')
    File.write('variable strainrate equal "%f"\n' % Strain_Rate)
    File.write('thermo_style custom step temp lx ly lz pxx density etotal\n')
    File.write('thermo 300\n')
    
    if (Dump != 'none'):
        File.write('dump 1 all custom 200 %s id type mol xs ys zs vx vy vz\n' % Dump)
        
    File.write('reset_timestep 0\n')
    File.write('timestep %d\n' % Time_Step)
    File.write('fix 1 all npt temp %d %d 100 y %d %d 1000 z %d %d 1000 drag 2\n' % ( Temp, Temp, Press, Press, Press, Press))
    File.write('fix 2 all deform 1 x erate %f remap x units box\n' % Strain_Rate)
    File.write('fix def1 all print 1 "${strainx}  ${strainy}  ${strainz}  ${stress}" file %s screen no\n' % Stress_Strain_File)
    File.write('run %d\n' % Num_Steps)
    File.write('write_restart %s\n' % Restart_Out)
    File.write('unfix 1\n')
    File.write('unfix 2\n')
    File.close()
    
    os.system('mpiexec -np 8 lmp_mpi -in %s' % In_File)
    return
    
    
    
    
    
    
    
    