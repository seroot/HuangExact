#! C:\Python27
import numpy as np
import math
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


"""
Auxillary Functions for Generating Initial Configuration
"""

def Gen_PDF_CDF(Kph, Beta):
	dx = 0.00001
	x = np.arange(0, 6.28, dx) # discretized Phi Values (0,2*pi)
	U = np.full_like(x, 0.0) # Potential Energy
	P = np.full_like(x, 0.0) # Boltzmann Probability
	L = len(x.tolist()) # Range
	CDF_NN = np.full_like(x, 0.0) # Non Normalized CDF
	CDF = np.full_like(x, 0.0) # Normalized CDF
	norm = 0
	C = np.array([6.3082, -0.69189991, 0.56789949, 3.56899986, -9.75319933], float)
	U = C[0]*np.cos(x)**0 + C[1]*np.cos(x)**1 + C[2]*np.cos(x)**2 + C[3]*np.cos(x)**3 + C[4]*np.cos(x)**4
	P = np.exp(-U*Beta)
	
	for i in range(L-1):
		CDF_NN[i+1] = CDF_NN[i] + P[i]*dx
	for i in range(L):
		P[i] = P[i]/CDF_NN[-1]
		norm += P[i]*dx
	for i in range(L-1):
		CDF[i+1] = CDF[i] + P[i]*dx
    
    return P, CDF
	
	
		
def Gen_Random_Dih(CDF, dx):
	"""
	Function for generating a random Dihedral angle by inverting the CDF
	input: Numpy array CDF
	Output: Random Float [0,2*pi]
	"""
	
	L = len(CDF.tolist())
	h = random.random()
	for i in range(L-1):
		if h > CDF[i] and h < CDF[i+1]:
			rand = i*dx
	return rand
 
 
 
def Demonstrate_InvertCDF(PDF, CDF, dx):
    L=10000
    freq = np.zeros([L])
    for i in range(L):
        freq[i] = Gen_Random_Dih(CDF, dx)
    x = np.arange(0,6.28,dx)
    print x
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title('Inverting the CDF')
    plt.xlabel('Dihedral Angle (Radians)')
    plt.ylabel('Frequency')
    n, bins, patches = ax.hist(freq, 50, normed=1, facecolor='green', alpha=0.75)
    plt.plot(x,PDF, 'r', linewidth=5)
    plt.xlim((0,x[-1]))
    #plt.plot(x,CDF,'b')
    plt.show()
    return
         



def Normalize(A):
	"""
	Function for normalizing a Numpy 3-array
	"""
	norm = 0
	for i in range(len(A.tolist())):
		norm += A[i]*A[i]
	B = A/(math.sqrt(norm))
	return B
	
def dotproduct(A,B):
	# Function to compute dotproduct
	C = A[0]*B[0] + A[1]*B[1] + A[2]*B[2]
	return C
	
def crossproduct(A,B):
	# Function to compute cross product
	C= np.zeros(3,float)
	C[0] = A[1]*B[2] - A[2]*B[1]
	C[1] = A[2]*B[0] - A[0]*B[2]
	C[2] = A[0]*B[1] - A[1]*B[0]
	return C
 
def Apply_PBC(Position, Box_Length):
    """ Apply Periodic Boundary Conditions to a set of particle positions
        
        input: 
               Position: Numpy array [N,3] Containing position coordinates of particles
               
               Box_Length: Length of Simulation box (float)
               
        Output: 
               PositionPBC: Numpy array [N,3] Containing remapped coordinates of particles
    """
    N = np.size(Position,0)
    M = np.size(Position,1)
    PositionPBC = np.zeros((N,3), float)
    for i in range(N):
        for j in range(M):
            PositionPBC[i,j] = Position[i,j] - Box_Length*math.floor(Position[i,j]/Box_Length)
    
    return PositionPBC


def Gen_rand_orthonormal_Basis():
    """
    This function generates a random orthonormal basis in the form of a 3x3 numpy array
    """
    xrsq = random.random()
    yrsq= (1 - xrsq)*random.random()
    zrsq = 1 - xrsq -yrsq
    yr = (random.randint(0,1)*2 - 1)*math.sqrt(yrsq)
    zr = (random.randint(0,1)*2 - 1)*math.sqrt(zrsq)
    xr = (random.randint(0,1)*2 - 1)*math.sqrt(xrsq)
    Basis = np.zeros((3,3))
    Basis[0] = [xr, yr, zr]
    x2r = random.random()
    y2r = -xr*x2r*random.random()
    z2r = -(xr*x2r + yr*y2r)/zr
    Basis[1] = [x2r,y2r, z2r]
    Basis[2] = crossproduct(Basis[0], Basis[1])
    return Basis
    
	

def generatepolymer( chainlength, R0, SigmaB, theta0, Sigmat, CUM2, dx, Box_Length):
    Bases = np.zeros([chainlength,3,3],float)
    Position_M = np.zeros([chainlength,3],float)
    Position_S1 = np.zeros([chainlength,3],float)
    Position_S2 = np.zeros([chainlength,3],float)
    Position_M[0] = [ Box_Length*random.random(), Box_Length*random.random(), Box_Length*random.random()]
    Bases[0] = [[1,0,0],[0,1,0],[0,0,1]]
    #Bases[0] = Gen_rand_orthonormal_Basis()
    Position_S1[0] = Position_M[0] + random.gauss(R0,SigmaB)*Bases[0,2]
    Position_S2[0] = Position_S1[0] + random.gauss(R0,SigmaB)*Bases[0,2]

    for i in range(chainlength-1):
        R = random.gauss(R0,SigmaB)
        Theta = 3.14 - random.gauss(theta0, Sigmat)
        Phi = Gen_Random_Dih(CUM2, dx)
        Position_M[i+1] = Position_M[i]
        Position_M[i+1] += R*(math.cos(Theta))*Bases[i,0] 
        Position_M[i+1] += R*(math.cos(Phi))*(math.sin(Theta))*Bases[i,1]
        Position_M[i+1] += R*(math.sin(Theta))*(math.sin(Phi))*Bases[i,2]
        Bases[i+1,0] = Normalize(Position_M[i+1] - Position_M[i])
        Bases[i+1,1] = Normalize(crossproduct(Bases[i+1,0],Bases[i,0]))
        Bases[i+1,2] = crossproduct(Bases[i+1,0],Bases[i+1,1])
        Position_S1[i+1] = Position_M[i+1] + random.gauss(R0,SigmaB)*Bases[i+1,2]
        Position_S2[i+1] = Position_S1[i+1] + random.gauss(R0,SigmaB)*Bases[i+1,2]
            #print dotproduct(Bases[i+1,0], Bases[i+1,1])
    #Position_M = Apply_PBC(Position_M, Box_Length)
    #Position_S1 = Apply_PBC(Position_S1, Box_Length)
    #Position_S2 = Apply_PBC(Position_S2, Box_Length)
            #print dotproduct(Bases[i+1,0], Bases[i+1,2])
		#print dotproduct(Bases[i+1,1], Bases[i+1,2])

    return Position_M, Position_S1, Position_S2, Bases

	

def Quiver_Plot(Position, Bases, Box_Length):
    """
	Function that takes an array of positions and their associated orthonormal basis sets to 
	generate a 3D vector plot of the polymer chains along with their respective orientations
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection= '3d')
    ax.quiver(Position[:,0], Position[:,1], Position[:,2], Bases[:,0,0], Bases[:,0,1], Bases[:,0,2], length = 3, color='r')
    ax.quiver(Position[:,0], Position[:,1], Position[:,2], Bases[:,1,0], Bases[:,1,1], Bases[:,1,2], length = 3, color='g')
    ax.quiver(Position[:,0], Position[:,1], Position[:,2], Bases[:,2,0],Bases[:,2,1], Bases[:,2,2], length = 3)
    #ax.scatter(Position[:,0], Position[:,1], Position[:,2],'k')
    #ax.grid(on=False)
    plt.title('Non-Overlapping Random Configuration, PBC\'s')
    #plt.xlim((0,Box_Length))
    #plt.ylim((0,Box_Length))
    ax.set_xlim(0,Box_Length)
    ax.set_ylim(0,Box_Length)
    ax.set_zlim(0,Box_Length)
    plt.show()
    
    return
    
def Scatter_Plot(Position_M, Position_P,Position_S, Box_Length):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(Position_M[:,0], Position_M[:,1], Position_M[:,2],c=u'g', s= 70)
    ax.scatter(Position_P[:,0], Position_P[:,1], Position_P[:,2], c=u'r', s= 30 )
    ax.scatter(Position_S[:,0], Position_S[:,1], Position_S[:,2], c=u'r', alpha = .2, s=30 )
    #ax.set_xlim(0,Box_Length)
    #ax.set_ylim(0,Box_Length)
    #ax.set_zlim(0,Box_Length)
    plt.show()
    
    return
    
    
def Center_of_Mass(Position):
    """
    Function to calculate center of mass coordinates of a single polymer chain
    input: Numpy array with positions of all the monomers
    output: COM = center of mass coordinates [3]
            ORIENT = Center of mass orientation [3][3]
            POSCOM = Positions of Monomers wrt COM [chainlength][3]
    """
    chainlength = len(Position.tolist())
    COM = np.mean(Position,axis=0)
    Orient = np.array([[1,0,0],[0,1,0],[0,0,1]], float)
    POSCOM = np.zeros(chainlength, 3)
    for i in range(chainlength):
        POSCOM[i] = Position[i] - COM
    return COM, Orient, POSCOM

            
            
    
def Check_For_Overlap(P1, P2, radius):
    """ 
    Function to check for overlap between two polymer chains
    input: P1 = Existing polymer chain
            P2 = Potential New Polymer Chain
            radius = radius... DUH
    Output: True indicates no overlap
            False indicates overlap
    """
    L1 = np.size(P1,0)
    L2 = np.size(P2,0)
    Rij = np.zeros(3)
    radiusSQ = radius*radius
    for i in range(L1):
        for j in range(L2):
            Rij = P1[i] - P2[j]
            RijSQ = dotproduct(Rij,Rij)
            if (RijSQ < radiusSQ):
                return False
    return True           
    
        

def Gen_Many_Polymers( ChainLength, NumChains, R0, SigmaB, theta0, Sigmat, CDF, dx, SigmaM_M, Box_Length):
    """ 
    Function to generate many non-overlapping polymer chains
    """
    Position_M, Position_S1, Position_S2, Bases = generatepolymer(ChainLength, R0, SigmaB, theta0, Sigmat, CDF, dx, Box_Length)
    NChains = 1
    while (NChains < NumChains):
        PositionNew_M, PositionNew_S1, PositionNew_S2, BasesNew = generatepolymer(ChainLength, R0, SigmaB, theta0, Sigmat, CDF, dx, Box_Length)
        Bool = Check_For_Overlap(Position_M, PositionNew_M, SigmaM_M)
        if (Bool):
            Position_M = Apply_PBC( np.concatenate((Position_M, PositionNew_M), axis=0), Box_Length)
            Position_S1 = Apply_PBC(np.concatenate((Position_S1, PositionNew_S1), axis=0), Box_Length)
            Position_S2 = Apply_PBC(np.concatenate((Position_S2, PositionNew_S2), axis=0), Box_Length)
            Bases = np.concatenate((Bases, BasesNew), axis=0)
            NChains += 1
            print "Polymer Deposited"
    return Position_M, Position_S1, Position_S2, Bases
    
    
def Deposit_PCBM( Position_M, SigmaM_P, SigmaP_P, NumPCBM, Box_Length, Bases):
    """
    Function to deposit PCBM Molecules randomly so that they dont overlap with the polymer molecules
    """
    
    NS = 0
    Pos = np.zeros(3,float)
    Position_C60 = np.zeros((NumPCBM,3), float)
    Position_C2H2 = np.zeros((NumPCBM,3), float)
    Position_C2H4 = np.zeros((NumPCBM,3), float)
    Position_Cnyl = np.zeros((NumPCBM,3), float)
    Position_C6r = np.zeros((NumPCBM,3), float)
    NumM = np.size(Position_M,0)
    SigmaM_PSQ = SigmaM_P*SigmaM_P
    SigmaP_PSQ = SigmaP_P*SigmaP_P
    while (NS<NumPCBM):
         Pos = [Box_Length*random.random() - 20, Box_Length*random.random() - 20,Box_Length*random.random() - 20 ]
         deposit = True
         for i in range(NumM):
             RijMP = Pos - Position_M[i]
             RijSQMP = dotproduct(RijMP,RijMP)
             if (RijSQMP < SigmaM_PSQ):
                 deposit = False
                 #print "Solvent Polymer Overlap"
         if (deposit and NS>0):
            for i in range(NS):
                RijPP = Pos - Position_C60[i]
                RijPPSQ = dotproduct(RijPP,RijPP)
                if (RijPPSQ < SigmaP_PSQ):
                    deposit = False
                    #print "Solvent Solvent Overlap"
        
         if (deposit):
            Basis = Bases[random.randint(0, NumM-1)]
            Position_C60[NS] = Pos
            Position_C2H2[NS] = Pos + 5.573*Basis[0]
            Position_C2H4[NS] = Position_C2H2[NS] + 3.*Basis[0]
            Position_Cnyl[NS] = Position_C2H4[NS] + 3.*Basis[1]
            Position_C6r[NS] = Position_C2H2[NS] - 3.*Basis[1]
            NS += 1
            print "PCBM Deposited"    
    return Apply_PBC(Position_C60, Box_Length), Apply_PBC(Position_C2H2, Box_Length), Apply_PBC(Position_C2H4, Box_Length), Apply_PBC(Position_Cnyl, Box_Length), Apply_PBC(Position_C6r, Box_Length)
            
def Deposit_Solvent(Position_M, Position_P, SigmaS_S, SigmaS_M, SigmaS_P, NumSolvent, Box_Length):
    
    NS = 0
    Pos = np.zeros(3,float)
    Position_S = np.zeros((NumSolvent,3),float)
    NumM = np.size(Position_M,0)
    NumP = np.size(Position_P,0)
    SigmaS_MSQ = SigmaS_M*SigmaS_M
    SigmaS_PSQ = SigmaS_P*SigmaS_P
    SigmaS_SSQ = SigmaS_S*SigmaS_S
    while (NS<NumSolvent):
        Pos = [Box_Length*random.random(), Box_Length*random.random(),Box_Length*random.random() ]
        deposit = True
        for i in range(NumM):
            RijMS = Pos - Position_M[i]
            RijMSSQ = dotproduct(RijMS,RijMS)
            if (RijMSSQ<SigmaS_MSQ):
                deposit = False
                #print "Solvent Polymer Overlap"
        if (deposit):
            for i in range(NumP):
                RijPS = Pos - Position_P[i]
                RijPSSQ = dotproduct(RijPS,RijPS)
                if (RijPSSQ<SigmaS_PSQ):
                    deposit = False
                    #print "Solvent PCBM Overlap"
        if (deposit and NS>0):
            for i in range(NS):
                RijSS = Pos - Position_S[i]
                RijSSSQ = dotproduct(RijSS,RijSS)
                if (RijSSSQ< SigmaS_SSQ):
                    deposit = False
                    #print "Solvent Solvent Overlap"
                    
        if (deposit):
            Position_S[NS]= Pos
            NS += 1
            print "Frac done", float(NS)/float(NumSolvent)
            
    return Position_S
        
        
    
         
    
    
	
	
	
	