import numpy as np
import math 
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D

def Gen_PDF_CDF_OPLS( V, Beta ):
	"""
	This function takes in a numpy array V that containes the energetic coefficients for the OPLS style dihedral potential of the form:
			U = (1/2)V1(1+cos(phi)) + (1/2)V2(1-cos(2phi)) + (1/2)V3(1+cos(3phi)) + ....
	It then uses Boltzmann statistics along with the inverse temperature Beta to generate a PDF and CDF of the dihedral angle
	
	The output is two numpy arrays that represent the PDF and CDF associated with this potential energy function
	"""
	dx = 0.0001
	x = np.arange(0, 6.28, dx) # Generate discretized values for x (phi)
	U = np.full_like(x, 0.0) # Initialize potential energy array
	PDF = np.full_like(x, 0.0) # Initialize PDF array
	CDF_NN = np.full_like(x, 0.0) # Initialize non-normalized CDF array
	CDF = np.full_like(x, 0.0) # Initialize normalized CDF array
	norm = 0
	L = len(x.tolist()) 
	U = 0.5*(V[0]*(1 + np.cos(x)) + V[1]*(1 - np.cos(2*x)) + V[2]*(1 + np.cos(3*x)) + V[3]*(1 - np.cos(4*x)))
	PDF = np.exp(-U*Beta)
	
	for i in range(L-1):
		CDF_NN[i+1] = CDF_NN[i] + PDF[i]*dx
	
	for i in range(L):
		PDF[i] = PDF[i]/CDF_NN[-1]
		norm += PDF[i]*dx
	
	for i in range(L-1):
		CDF[i+1] = CDF[i] + PDF[i]*dx
		
	return PDF, CDF
 
def Gen_PDF_CDF_Multi( V, Beta ):
	"""
	This function takes in a numpy array V that containes the energetic coefficients for the Multi style dihedral potential of the form
	It then uses Boltzmann statistics along with the inverse temperature Beta to generate a PDF and CDF of the dihedral angle
	
	The output is two numpy arrays that represent the PDF and CDF associated with this potential energy function
	"""
	dx = 0.0001
	x = np.arange(0, 6.28, dx) # Generate discretized values for x (phi)
	U = np.full_like(x, 0.0) # Initialize potential energy array
	PDF = np.full_like(x, 0.0) # Initialize PDF array
	CDF_NN = np.full_like(x, 0.0) # Initialize non-normalized CDF array
	CDF = np.full_like(x, 0.0) # Initialize normalized CDF array
	norm = 0
	L = len(x.tolist()) 
	U = V[0] + V[1]*np.cos(x) + V[2]*np.cos(x)**2 + V[3]*np.cos(x)**3 + V[4]*np.cos(x)**4
	PDF = np.exp(-U*Beta)
	
	for i in range(L-1):
		CDF_NN[i+1] = CDF_NN[i] + PDF[i]*dx
	
	for i in range(L):
		PDF[i] = PDF[i]/CDF_NN[-1]
		norm += PDF[i]*dx
	
	for i in range(L-1):
		CDF[i+1] = CDF[i] + PDF[i]*dx
		
	return PDF, CDF
	
	
	
	
def Plot_Dihedral( V, style):
    """ 
    This function plots the dihedral energy function from 0 to 2*pi
    inputs:
        V is a numpy array containing the parameters for the functional form of the dihedral potential
        style can either be 'OPLS' or 'MULTI'
    returns nothing
    """
    dx = 0.0001
    x = np.arange(0, 6.28, dx)
    U = np.full_like(x, 0.0)    
    
    if style == 'OPLS':
        U = 0.5*(V[0]*(1 + np.cos(x)) + V[1]*(1 - np.cos(2*x)) + V[2]*(1 + np.cos(3*x)) + V[3]*(1 - np.cos(4*x)))
        
    elif style == 'MULTI':
        print 'Under Construction! Sorry'
        
    
    plt.figure()
    plt.plot(x, U)
    plt.title('Dihedral Potential')
    plt.xlabel('Dihedral Angle (rad)')
    plt.ylabel('Energy (kcal/mol)')
    plt.xlim((0,6.28))
    plt.show()
    return
    
def Compare_Bond( Bond ):
    dx = 0.0001
    N = len(Bond) - 1
    r = np.arange(Bond[0]-.5, Bond[0]+.5, dx)
    U1 = np.full_like(r, 0.0)
    U2 = np.full_like(r, 0.0)
    
    for i in range(1,N+1):
        U1 += Bond[i]*(r - Bond[0])**(i+1)
        print i+1
        
    U2 = Bond[1]*(r - Bond[0])**2
    
    plt.figure()
    plt.plot(r, U1, label = 'Full Potential')
    plt.plot(r, U2, label = 'Harmonic Approximation')
    plt.xlim((Bond[0]-.5,Bond[0]+.5))
    plt.ylim((U1.min(), U1.max()))
    plt.title('P2P3', fontsize = 30)
    plt.xlabel('Bond Length (Angstrom)', fontsize = 20)
    plt.ylabel('Potential Energy (Kcal/mol)', fontsize = 20)
    plt.legend()
    plt.show()
    return
    
def Compare_Angle( Angle):
    dx = 0.0001
    M = [0, 2 , 3, 4 , 5 , 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
    N = len(Angle) - 1
    Th0 = Angle[0]*(3.1415/180.)
    dTh = .5
    Th = np.arange(Th0-dTh, Th0 + dTh, dx)
    U1 = np.full_like(Th, 0.0)
    U2 = np.full_like(Th, 0.0)
    
    for i in range(1, N+1):
        U1 += Angle[i]*(Th - Th0)**M[i-1]
    
    U2 = Angle[2]*(Th - Th0)**2 + Angle[3]*(Th - Th0)**3 + Angle[4]*(Th-Th0)**4
    plt.figure()
    plt.plot(Th, U1, label = 'Full Potential')
    plt.plot(Th, U2, label = 'Class 2 Approximation')
    plt.ylim((U2.min(), U1.max() ))
    plt.xlim((Th0-dTh, Th0 + dTh))
    plt.title('P2P1P1')
    plt.xlabel('Angle (Radians)', fontsize = 20)
    plt.ylabel('Potential Energy (Kcal/mol)', fontsize = 20)
    plt.legend()
    plt.show()

    return  

def Compare_Dihedral(Dihedral):
    dx = 0.0001
    N = len(Dihedral)
    Dih = np.arange(0, 2*3.1415, dx)
    U1 = np.full_like(Dih, 0.0)
    U2 = np.full_like(Dih, 0.0)
    
    for i in range(N):
        U1 += Dihedral[i]*np.cos(Dih)**i
    
    for i in range(4):
        U2 += Dihedral[i]*np.cos(Dih)**i
    
    plt.figure()
    plt.plot(Dih, U1, label = 'Full Potential')
    plt.plot(Dih, U2, label = 'Truncated Approximation')
    plt.ylim((U2.min(), U2.max()))
    plt.xlim((0.0, 2*3.1415))
    plt.title('P2P1P1P2', fontsize = 30)
    plt.xlabel('Dihedral Angle (radians)', fontsize = 20)
    plt.ylabel('Energy (Kcal/mol)', fontsize = 20)
    plt.legend()
    plt.show()
    return
    
            
        
    

def Gen_Random_Dih(CDF):
	"""
	Function for generating a random Dihedral angle by inverting the CDF
	input: Numpy array CDF
	Output: Random Float [0,2*pi]
	"""
	dx = 0.0001
	L = len(CDF.tolist())
	h = random.random()
	for i in range(L-1):
		if h > CDF[i] and h < CDF[i+1]:
			rand = i*dx
	return rand
 
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
    
def Gen_Linear_Polymer(ChainLength, r0, th0, CDF, Box_Length):
    Bases = np.zeros([ChainLength,3,3], float)
    Position = np.zeros([ChainLength,3], float)
    Position[0] = [ Box_Length*random.random(), Box_Length*random.random(), Box_Length*random.random()]
    Bases[0] = Gen_rand_orthonormal_Basis()
    
    
    for i in range(ChainLength-1):
        Phi = Gen_Random_Dih(CDF)
        Position[i+1] = Position[i]
        Position[i+1] += r0*(math.cos(th0+3.14))*Bases[i,0]
        Position[i+1] += r0*(math.sin(th0+3.14))*(math.cos(Phi))*Bases[i,1]
        Position[i+1] += r0*(math.sin(th0))*(math.sin(Phi))*Bases[i,2]
        Bases[i+1, 0] = Normalize(Position[i+1] - Position[i])
        Bases[i+1, 1] = Normalize(crossproduct(Bases[i+1,0], Bases[i,0]))
        Bases[i+1, 2] = crossproduct(Bases[i+1,0], Bases[i+1,1])
        
    Position = Apply_PBC(Position, Box_Length)
    
    return Position, Bases
    
def generatepolymer( chainlength, BondLengths , Angles, CUM2, Box_Length):
    Bases = np.zeros([chainlength,3,3],float)
    Position_M = np.zeros([chainlength,3],float)
    Position_S1 = np.zeros([chainlength,3],float)
    Position_S2 = np.zeros([chainlength,3],float)
    Position_M[0] = [ Box_Length*random.random(), Box_Length*random.random(), Box_Length*random.random()]
    #Bases[0] = [[1,0,0],[0,1,0],[0,0,1]]
    Bases[0] = Gen_rand_orthonormal_Basis()
    Position_S1[0] = Position_M[0] + BondLengths['P1P2']*Bases[0,2]
    Position_S2[0] = Position_S1[0] + BondLengths['P2P3']*Bases[0,2]

    for i in range(chainlength-1):
        R = BondLengths['P1P1']
        Theta = 3.14 - Angles['P1P1P1']*(3.14/180.0)
        Phi = Gen_Random_Dih(CUM2)
        Position_M[i+1] = Position_M[i]
        Position_M[i+1] += R*(math.cos(Theta))*Bases[i,0] 
        Position_M[i+1] += R*(math.cos(Phi))*(math.sin(Theta))*Bases[i,1]
        Position_M[i+1] += R*(math.sin(Theta))*(math.sin(Phi))*Bases[i,2]
        Bases[i+1,0] = Normalize(Position_M[i+1] - Position_M[i])
        Bases[i+1,1] = Normalize(crossproduct(Bases[i+1,0],Bases[i,0]))
        Bases[i+1,2] = crossproduct(Bases[i+1,0],Bases[i+1,1])
        Position_S1[i+1] = Position_M[i+1] + BondLengths['P1P2']*Bases[i+1,2]
        Position_S2[i+1] = Position_S1[i+1] + BondLengths['P2P3']*Bases[i+1,2]
            #print dotproduct(Bases[i+1,0], Bases[i+1,1])
    Position_M = Apply_PBC(Position_M, Box_Length)
    Position_S1 = Apply_PBC(Position_S1, Box_Length)
    Position_S2 = Apply_PBC(Position_S2, Box_Length)
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
    ax.scatter(Position[:,0], Position[:,1], Position[:,2],'k')
    #ax.grid(on=False)
    plt.title('Non-Overlapping Random Configuration, PBC\'s')
    plt.xlim((0,Box_Length))
    plt.ylim((0,Box_Length))
    plt.show()
    return



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
                print 'Overlap'
                return False, Rij
    return True     
    
def Gen_Many_Lin_Polymers( ChainLength, NumChains, r0, th0, CDF, Sigma, Box_Length):
    """ 
    Function to generate many non-overlapping polymer chains
    """
    Position, Bases = Gen_Linear_Polymer(ChainLength, r0, th0, CDF,  Box_Length)
    NChains = 1
    while (NChains < NumChains):
        PositionNew, BasesNew = Gen_Linear_Polymer(ChainLength, r0, th0, CDF, Box_Length)
        Bool = Check_For_Overlap(Position, PositionNew, Sigma)
        if (Bool):
            Position = Apply_PBC( np.concatenate((Position, PositionNew), axis=0), Box_Length)
            Bases = np.concatenate((Bases, BasesNew), axis=0)
            NChains += 1
            print "Polymer Deposited %d " % NChains
    return Position, Bases
    


def Gen_Many_Polymers( ChainLength, NumChains, BondLengths, Angles, CUM2, SigmaM_M, Box_Length):
    """ 
    Function to generate many non-overlapping polymer chains
    """
    Position_M, Position_S1, Position_S2, Bases = generatepolymer( ChainLength, BondLengths , Angles, CUM2, Box_Length)
    NChains = 1
    while (NChains < NumChains):
        PositionNew_M, PositionNew_S1, PositionNew_S2, BasesNew = generatepolymer( ChainLength, BondLengths , Angles, CUM2, Box_Length)
        Bool = Check_For_Overlap(Position_M, PositionNew_M, SigmaM_M)
        if (Bool):
            Position_M = Apply_PBC( np.concatenate((Position_M, PositionNew_M), axis=0), Box_Length)
            Position_S1 = Apply_PBC(np.concatenate((Position_S1, PositionNew_S1), axis=0), Box_Length)
            Position_S2 = Apply_PBC(np.concatenate((Position_S2, PositionNew_S2), axis=0), Box_Length)
            Bases = np.concatenate((Bases, BasesNew), axis=0)
            NChains += 1
            print "Polymer Deposited %d" % NChains
        
    return Position_M, Position_S1, Position_S2, Bases
    
	