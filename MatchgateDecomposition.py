import numpy as np
import numpy.linalg as npl
import math
from scipy.stats import ortho_group,unitary_group
from scipy.linalg import expm
import HoffmanDecomposition as hof


def Tr(M):
    return sum(np.diag(M))

def RotationToMatchgates(T): 
    """
    Takes an element T of SO(2N) and outputs a list of nearest-neighbor matchgates that implement T.
    Output is given by a list of ((i,j),U), where (i,j) denote two nearest-neighbor sites and U denotes the matchgate to apply to them.
    """
    N = int(T.shape[0]/2)
    I,X,Y,Z=np.matrix(np.identity(2)),np.matrix([[0j,1],[1,0]]),np.matrix([[0,-1j],[1j,0]]),np.matrix([[1,0],[0,-1]]) # Pauli matrices
    SE=np.matrix(expm(1j*math.pi/4*np.kron(Y,X)))
    SO=np.matrix(expm(-1j*math.pi/4*np.kron(X,Y)))

    anglesList = hof.FindAngles(T)
    matchgateList = []
    sampleMajoranas = [np.kron(X,np.kron(I,I))/2**1.5,np.kron(Y,np.kron(I,I))/2**1.5,np.kron(Z,np.kron(X,I))/2**1.5,np.kron(Z,np.kron(Y,I))/2**1.5,np.kron(np.kron(Z,Z),X)/2**1.5,np.kron(np.kron(Z,Z),Y)/2**1.5]
    for q in range(2*N-1,-1,-1):
        for p in range(q-1,-1,-1):
            theta=anglesList[q][p]
            A,B = int(p/2),int(q/2) # The qubit indices corresponding to matchgate indices p and q
            if A==B: # If the indices are on the same site, no swaps are needed
                if A==0:
                    U = np.kron(np.matrix(np.diag([math.e**(1j*theta/2),math.e**(-1j*theta/2)])),np.identity(2))
                    pair = (0,1)
                    rot = np.block([[np.matrix([[math.cos(theta),math.sin(theta)],[-math.sin(theta),math.cos(theta)]]),np.zeros((2,2))],[np.zeros((2,2)),np.identity(2)]])
                else:
                    U = np.kron(np.identity(2),np.matrix(np.diag([math.e**(1j*theta/2),math.e**(-1j*theta/2)])))
                    pair = (A-1,A)
                    rot = np.block([[np.identity(2),np.zeros((2,2))],[np.zeros((2,2)),np.matrix([[math.cos(theta),math.sin(theta)],[-math.sin(theta),math.cos(theta)]])]])
                matchgateList.append((pair,U))
            else:
                for i in range(A,B-1): # Add swap operators to bring the rotation to neighboring sites
                    if p%2==0:
                        matchgateList.append(((i,i+1),SE))
                    else:
                        matchgateList.append(((i,i+1),SO))
                # Now, perform the nearest-neighbor swap. Note each swap operator results in sending theta to (-theta), so we include a sign in our exponential
                if p%2==0 and q%2==0:
                    U = expm((-1)**(B-A)*1j*theta/2*np.kron(Y,X))
                elif p%2==0 and q%2==1:
                    U = expm((-1)**(B-A)*1j*theta/2*np.kron(Y,Y))
                elif p%2==1 and q%2==0:
                    U = expm((-1)**(B-A-1)*1j*theta/2*np.kron(X,X))
                elif p%2==1 and q%2==1:
                    U = expm((-1)**(B-A-1)*1j*theta/2*np.kron(X,Y))
                matchgateList.append(((B-1,B),U))
                for i in range(B-2,A-1,-1): # Apply the inverse swap operators
                    if p%2==0:
                        matchgateList.append(((i,i+1),SE.H))
                    else:
                        matchgateList.append(((i,i+1),SO.H))
    return matchgateList 

def FindMajoranaOperators(N):
    """
    Given N sites, finds the 2*N Majorana operators associated to the sites
    """
    I,X,Y,Z=np.matrix(np.identity(2)),np.matrix([[0j,1],[1,0]]),np.matrix([[0,-1j],[1j,0]]),np.matrix([[1,0],[0,-1]])
    majoranaList = [] # Create the list of normalized 2N majorana operators
    for i in range(N):
        c_2i = np.matrix([[1]])/2**(N/2)
        c_2iplus1= np.matrix([[1]])/2**(N/2)
        for j in range(i):
            c_2i=np.kron(Z,c_2i)
            c_2iplus1=np.kron(Z,c_2iplus1)
        c_2i=np.kron(c_2i,X)
        c_2iplus1=np.kron(c_2iplus1,Y)
        for j in range(i+1,N):
            c_2i=np.kron(c_2i,I)
            c_2iplus1=np.kron(c_2iplus1,I)
        majoranaList.append(c_2i)
        majoranaList.append(c_2iplus1)
    return majoranaList
 
def MatchgatesToRotation(matchgateList,majoranaList=[]):
    """
    Given a list of matchgates U and sites (i,j), produce the rotation matrix
    This can either take a precomputed list of Majorana operators (efficient if run multiple times)
    or will compute the Majorana operators at runtime (efficient if only run once)
    """
    N = max([max(m[0]) for m in matchgateList])+1
    I,X,Y,Z=np.matrix(np.identity(2)),np.matrix([[0j,1],[1,0]]),np.matrix([[0,-1j],[1j,0]]),np.matrix([[1,0],[0,-1]])
    if len(majoranaList)!=2*N: # If there is no majoranaList given, find the majoranaList
        majoranaList = FindMajoranaOperators(N)
   
    totalUnitary = np.matrix(np.identity(2**N))
    for (a,b),M in matchgateList:
        U = np.matrix([[1]])
        for i in range(a):
            U = np.kron(I,U)
        U = np.kron(U,M)
        for i in range(b+1,N):
            U = np.kron(U,I)
        totalUnitary=totalUnitary*U
    rotationMatrix = np.matrix([[sum(np.diag(c2.H*totalUnitary*c1*totalUnitary.H)) for c1 in majoranaList] for c2 in majoranaList])
    return rotationMatrix

def TestMatchgateFunctions(numSamples,N):
    """
    Generates numSamples elements of SO(2*N), and checks that MatchgatesToRotation(RotationToMatchgates(T))=T
    """
    majoranaList = FindMajoranaOperators(N)
    for i in range(numSamples): 
        T = ortho_group.rvs(2*N)
        T = np.matrix(np.diag([npl.det(T)]+[1]*(2*N-1)))*T #Converts an element of O(2*N) to SO(2*N)
        assert (np.max(abs(MatchgatesToRotation(RotationToMatchgates(T),majoranaList=majoranaList)-T))<1e-10)

TestMatchgateFunctions(100,3)
