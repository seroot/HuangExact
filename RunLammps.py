# -*- coding: utf-8 -*-
"""
Created on Sat Apr 04 11:52:08 2015

Module containing functions for running lammps simulations with given input parameters

@author: Samuel
"""

import os

def EquilibrateNVT(In_File, Restart_Out,  Data_File = 'none', Restart_In = 'none', Dump = 'none', Num_Steps = 10000, Time_Step = 1, units = 'real',  Temp = 423.0):
    """
    Function for generating an input file for an NVT equilibration simulation
    
    Input Arguments:
        Required:
        In_File = name of the input file to be created
        Restart_Out = name of restart file to be outputed
        Data_File | Restart_In = File for declaring initial position and topology (atleast on of these is required)
    Optional Arguments:
        Dump = name of dumpfile to be outputed, Default is none
        Num_Steps = number of timesteps to integrate equations of motion
        Time_Step = size of time increment
        Units = units to be used
        Temp = Temperature to run simulation at
    """
        
    
    File = open(In_File, 'w')
    File.write( '# Input file for running an NVT Equilibration in LAMMPS, Filename = %s\n\n' % In_File)
    
    File.write('units %s\n' % units)
    File.write('atom_style molecular\n')
    File.write('boundary p p p\n')
    File.write('bond_style harmonic\n')
    File.write('pair_style lj/cut 25\n')
    File.write('angle_style harmonic\n')
    File.write('dihedral_style opls\n')
    File.write('special_bonds lj 0 1 1\n')
    File.write('improper_style none\n')
    File.write('kspace_style none\n')
    
    if (Data_File == 'none' and Restart_In == 'none'):
        print '##Error: You must specify either a data file or a restart file. Try again'
    
    if (Data_File == 'none' and Restart_In != 'none'):
        File.write('read_restart %s\n' % Restart_In)
    
    if (Restart_In == 'none' and Data_File != 'none' ):
        File.write('read_data %s\n' % Data_File)
    
    File.write('thermo_style custom step temp press ke pe etotal density\n')
    
    if (Dump != 'none'):
        File.write('dump 1 all custom 200 %s id type mol xs ys zs vx vy vz\n' % Dump)
    
    File.write('neighbor 10.0 bin\n')
    File.write('neigh_modify every 1 delay 0 one 1000\n')
    
    File.write('fix 1 all nve/limit 0.5\n')
    
    if (Restart_In == 'none' and Data_File != 'none' ):
        File.write('velocity all create %f 1223\n' % Temp)
    
    File.write('fix 2 all temp/berendsen %f %f 15\n' % (Temp, Temp))
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
    File.write('fix 1 all npt temp %d %d 100 iso %d %d 1000 drag 2\n' % (Temp, Temp, Press, Press))
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
    
    

def RampTemp( In_File, Restart_In, Restart_Out, Tstart, Tend, Dump = 'none', Num_Steps = 10000, Time_Step = 1, Press = 1):
    """
    Function for ramping the temperature of a simulation in LAMMPS
    
    Required Inputs:
        In_File = Name of simulation input file
        
        Restart_In = Name of restart file to read in
        
        Restart_Out = Name of restart file to output
        
        Tstart = starting temperature
        
        Tend = ending temperature
    
    Optional Inputs: 
        Dump = name of dumpfile to output
        
        Num_Steps = number of steps to integrate equations of motion
        
        Time_Step = Time increment to integrate equations of motion
        
        Press = external pressure 
        
    """
    
    File = open(In_File, 'w')
    
    File.write('# Input file for LAMMPS SImulation to ramp the temperature, Filename = %s\n\n' % In_File)
    
    File.write('read_restart %s\n' % Restart_In)
    File.write('fix 1 all npt temp %d %d 100 iso %d %d 1000 drag 2\n' % (Tstart, Tend, Press, Press))
    File.write('fix 2 all momentum 1 linear 1 1 1\n')
    File.write('thermo_style custom step temp press etotal density vol\n')
    File.write('thermo 200\n')
    
    if (Dump != 'none'):
        File.write('dump 1 all custom 200 %s id type mol xs ys zs vx vy vz\n' % Dump)
        
    File.write('timestep %d\n' % Time_Step)
    File.write('run %d\n' % Num_Steps)
    File.write('unfix 1\n')
    File.write('unfix 2\n')
    File.write('write_restart %s' % Restart_Out)
    File.close()
    
    os.system('mpiexec -np 8 lmp_mpi -in %s' % In_File)
    
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
    
    

    
    
    
    
    
    
        

    
    