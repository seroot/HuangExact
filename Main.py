# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 15:46:22 2015

@author: Samuel
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os




import PolyModelFunctionsHuang as poly
import Huangparameters as Huang
import HuangData
import RunLammps



def main():
    
    
    #poly.Compare_Dihedral(Huang.P2P1P1P2)
    
    PDF, CDF = poly.Gen_PDF_CDF_Multi(Huang.P1P1P1P1, Huang.Beta)
    #Position_M, Position_S1, Position_S2, Bases = poly.generatepolymer(Huang.ChainLength, Huang.BondLengths, Huang.Angles, CDF, Huang.Box_Length)
    #print Position_M
    
    Position_M, Position_S1, Position_S2, Bases  =  poly.Gen_Many_Polymers( Huang.ChainLength, Huang.NumChains, Huang.BondLengths, Huang.Angles, CDF, Huang.SigmaM_M, Huang.Box_Length)
    In_Filename = 'in.P3HTHuangE%d_%d' % ( Huang.ChainLength, Huang.NumChains)
    Data_Filename = 'data.P3HTHuang%d_%d' % ( Huang.ChainLength, Huang.NumChains)
    Restart_Filename = 'restart.P3HTHuang%d_%d' % (Huang.ChainLength, Huang.NumChains)
    Dump_Filename = 'P3HTHuangCheck.lammpstrj'
    Pair_Table = 'HuangPair.table'
     
    os.mkdir('P3HTHuang%d_%d' %(Huang.ChainLength, Huang.NumChains))
    os.system('copy HuangPair2.table C:\Users\Samuel\Desktop\Sam\Simulation_Code\Publication1\ApproximateHuangModel\P3HTHuang%d_%d"' % (Huang.ChainLength, Huang.NumChains)); os.system('copy HuangBond.table C:\Users\Samuel\Desktop\Sam\Simulation_Code\Publication1\ApproximateHuangModel\P3HTHuang%d_%d"' % (Huang.ChainLength, Huang.NumChains)); os.system('copy HuangAngle.table C:\Users\Samuel\Desktop\Sam\Simulation_Code\Publication1\ApproximateHuangModel\P3HTHuang%d_%d"' % (Huang.ChainLength, Huang.NumChains))
	
    os.chdir('P3HTHuang%d_%d' %(Huang.ChainLength, Huang.NumChains))
    
    HuangData.P3HTHuangDataLammps(Position_M, Position_S1, Position_S2, Data_Filename)
    
    HuangData.Run_Huang_Equil(In_Filename, Data_Filename, Pair_Table, Restart_Filename, Dump = Dump_Filename, Num_Steps= 100000, Time_Step = 2)
    
    infileNPT = 'in.P3HTHuang%d_%d_NPT' % (Huang.ChainLength, Huang.NumChains)
    Restart_Out = 'restart.P3HTHuang%d_%d_E' %(Huang.ChainLength, Huang.NumChains)
    Dump_Filename2 = 'P3HTHuang%d_%d_compress.lammpstrj' %(Huang.ChainLength, Huang.NumChains)
    HuangData.EquilibrateNPT(infileNPT, Restart_Filename, Restart_Out, Dump=Dump_Filename2, Time_Step = 2, Num_Steps = 1000000);"""
  
    InFileLP = 'in.PersistenceLength'
    RestartLP = 'restart.P3HTHuangLP'
    
    HuangData.Calculate_LP(InFileLP , Restart_Out, RestartLP , Time_Step = 4)
    

    
   
    (PDF, CDF) = poly.Gen_PDF_CDF_OPLS(Lee.DIH_OPLS_LEE, Lee.Beta)
    Position, Bases = poly.Gen_Many_Lin_Polymers(Lee.ChainLength, Lee.NumChains, Lee.r0, Lee.th0, CDF, Lee.SigmaPP, Lee.Box_Length)
    
    Data_Filename = 'data.P3HTELee100_100'
    In_Filename = 'in.P3HTLeeE100_100'
    
    Restart_Filename = 'restart.P3HTLeeE100_100'
    Dump_Filename = 'P3HTLEECheck.lammpstrj'
    
    os.mkdir('P3HTLee%d_%d' % (Lee.ChainLength, Lee.NumChains))
    os.chdir('P3HTLee%d_%d' % (Lee.ChainLength, Lee.NumChains))
    
    Wlammps.WriteToLammps1(Position, Data_Filename)
    
    RunLammps.EquilibrateNVT(In_Filename, Restart_Filename, Data_File = Data_Filename, Dump = Dump_Filename, Time_Step = 10)
    
    FilenameC = 'in.P3HTLeeC%d_%d' % (Lee.ChainLength, Lee.NumChains)
    Restart_OutC = 'restart.P3HTLeeC%d_%d' % (Lee.ChainLength, Lee.NumChains)
    Dump_FileC = 'P3HTLee%d_%dC.lammpstrj' % (Lee.ChainLength, Lee.NumChains)
    RunLammps.EquilibrateNPT( FilenameC, Restart_Filename, Restart_OutC, Dump= Dump_FileC, Time_Step = 10, Num_Steps = 1000000)
    
    FilenameA = 'in.P3HTLeeA%d_%d' % (Lee.ChainLength, Lee.NumChains)
    Restart_OutA = 'restart.P3HTLee%d_%dA' % (Lee.ChainLength, Lee.NumChains)
    RunLammps.RampTemp(FilenameA, Restart_OutC, Restart_OutA, 423, 300, Press = 0, Time_Step = 10, Num_Steps = 1000000)
    
    FilenameS = 'in.P3HTLeeS%d_%d' % (Lee.ChainLength, Lee.NumChains)
    Restart_OutS = 'restart.P3HTLeeS%d_%d' % (Lee.ChainLength, Lee.NumChains)
    StressStrain = 'P3HTLee%d_%d_SS.txt' % (Lee.ChainLength, Lee.NumChains)
    DumpFile = 'P3HTLee%d_%d_Strain.lammpstrj' % (Lee.ChainLength, Lee.NumChains)
    
    RunLammps.Strain( FilenameS, Restart_OutA, Restart_OutS, StressStrain, Dump = DumpFile, Temp = 300, Time_Step = 10, Num_Steps = 500000)
    """
    
    
        
    
    
    
    """
    Wlammps.Write_Lammps_Equil_Input1(In_Filename, Data_Filename, Restart_Filename, Dump_Filename)
    os.system('mpiexec -np 8 lmp_mpi -in %s'% In_Filename)
    
    FilenameC = 'in.P3HTLeeC100_100'
    Restart_OutC = 'restart.P3HTLee100_100C'
    Dump_FileC = 'P3HTLee100_100C.lammpstrj'
    
    Wlammps.Write_Lammps_Condense(FilenameC, Restart_Filename, Restart_OutC, Dump_FileC)
    
    Wlammps.Write_Lammps_Condense(FilenameC, Restart_Filename, Restart_OutC, Dump_FileC)
    os.system('mpiexec -np 8 lmp_mpi -in %s' % FilenameC)
    """
    
    
    

if __name__=='__main__': main()