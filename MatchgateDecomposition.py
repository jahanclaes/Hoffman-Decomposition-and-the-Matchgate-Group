import numpy as np
import numpy.linalg as npl
import math
from scipy.stats import ortho_group,unitary_group
import SpecialOrthogonalDecomposition as sod


S = np.matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,-1]])
for i in range(N-1):
    SI = S
    for j in range(i):
        SI = np.kron(np.identity(2),SI)
    for j in range(i+2,N):
        SI = np.kron(SI,np.identity(2))
    SList.append(SI)
 
def RotationToUnitary(T): #Takes an element T of SO(2N) and outputs a list of nearest-neighbor matchgates that implement T. Output is given by a list of ((i,j),U), where (i,j) denote two nearest-neighbor sites and U denotes the matchgate to apply to them.
    anglesList = FindAngles(T)
    matchgateList = []
    for q in range(2*self.N-1,-1,-1):
        for p in range(q-1,-1,-1):
            theta=anglesList[q][p]
            A,B = int(p/2),int(q/2)
            if A==B:
                if A==0:
                    U = np.kron(np.matrix(np.diag([math.e**(1j*theta/2),math.e**(-1j*theta/2)])),np.identity(2))
                    pair = (0,1)
                else:
                    U = np.kron(np.identity(2),np.matrix(np.diag([math.e**(1j*theta/2),math.e**(-1j*theta/2)])))
                    pair = (A-1,A)
                matchgateList.append((pair,U))
            else:
                for i in range(A,B-1):
                    totalUnitary=totalUnitary*self.SList[i]
                if p%2==0 and q%2==0:
                    U = expm(-1j*theta/2*np.kron(self.Y,self.X))
                elif p%2==0 and q%2==1:
                    U = expm(-1j*theta/2*np.kron(self.Y,self.Y))
                elif p%2==1 and q%2==0:
                    U = expm(1j*theta/2*np.kron(self.X,self.X))
                elif p%2==1 and q%2==1:
                    U = expm(1j*theta/2*np.kron(self.X,self.Y))
                for i in range(B-1):
                    U = np.kron(np.identity(2),U)
                for i in range(B+1,self.N):
                    U = np.kron(U,np.identity(2))
                totalUnitary=totalUnitary*U
                for i in range(B-2,A-1,-1):
                    totalUnitary=totalUnitary*self.SList[i]
    return totalUnitary 

