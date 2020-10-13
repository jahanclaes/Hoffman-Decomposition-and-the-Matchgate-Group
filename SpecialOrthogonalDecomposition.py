import numpy as np
import numpy.linalg as npl
import math
from scipy.stats import ortho_group

"""
This code implements the Hoffman decomposition of a matrix T in SO(N) into the product of N(N-1)/2 rotations that only involve 2 basis elements each. This decomposition requires describing a matrix T by its generalized Euler angles, {theta(k,nu) : nu=2...N,k=1...nu-1}.

The parameterization of a matrix using generalized Euler angles, and the formula for constructing T from the list of Euler angles, was first described in Raffenetti and Rudenberg, Int. J. Quantum Chem.4, 625 (1969).

The inverse formula for constructing Euler angles from a matrix T, and the decomposition of T into rotations that only involve 2 basis elements, was first described in Hoffman, Raffenetti, and Rudenberg, J. Math. Phys.13, 528 (1972).

"""


def FindAngles(T): #Find the Euler angles. Outputs [[],[theta(1,2)],...,[theta(1,N),...,theta(N-1,N)]]
    N = T.shape[0]
    anglesList = []    
    TNu=T
    for nu in range(N-1,0,-1):
        thetaK = np.arctan2(TNu[nu-1,nu],TNu[nu,nu])
        anglesListNu = [thetaK]
        for k in range(nu-2,-1,-1):
            thetaK = np.arctan2(TNu[k,nu],TNu[k+1,nu]/math.sin(anglesListNu[0]))
            if np.isnan(thetaK):
                thetaK=0.
            anglesListNu = [thetaK]+anglesListNu
        FNu = np.matrix(np.zeros((N,N)))
        for n in range(N):
            FNu[nu,n]=TNu[nu,n]
            for k in range(nu-1,-1,-1):
                FNu[k,n]=math.sin(anglesListNu[k])*TNu[k,n]+math.cos(anglesListNu[k])*FNu[k+1,n]
        TNuMinusOne = np.matrix(np.zeros((N,N)))
        for n in range(N):
            for k in range(nu):
                TNuMinusOne[k,n]=math.cos(anglesListNu[k])*TNu[k,n]-math.sin(anglesListNu[k])*FNu[k+1,n]
            if n>=nu:
                TNuMinusOne[n,n]=1
        TNu=TNuMinusOne
        anglesList = [anglesListNu]+anglesList
    return [[]]+anglesList

def FindMatrix(anglesList): #Given input [[],[theta(1,2)],...,[theta(1,N),...,theta(N-1,N)]], output the corresponding matrix T
    N = len(anglesList)
    TNu = np.matrix(np.identity(N))
    for nu in range(1,N):
        FNu = np.matrix(np.zeros((N,N)))
        FNu[0,nu]=1
        for k in range(nu):
            for n in range(N):
                FNu[k+1,n]=-math.sin(anglesList[nu][k])*TNu[k,n]+math.cos(anglesList[nu][k])*FNu[k,n]
        TNuPlusOne = np.matrix(np.zeros((N,N)))
        for n in range(N):
            for k in range(nu):
                TNuPlusOne[k,n] = math.cos(anglesList[nu][k])*TNu[k,n]+math.sin(anglesList[nu][k])*FNu[k,n]
            TNuPlusOne[nu,n]=FNu[nu,n]
            if n>nu:
                TNuPlusOne[n,n]=1
        TNu=TNuPlusOne 
    return TNu

def FindHoffmanDecomposition(T): #Given a matrix T, output a list of rotations [A_(N(N-1)/2),...,A_1] such that T = A_(N(N-1)/2)x...xA_1 and each A_i only involves two basis elements
    N = T.shape[0]
    anglesList = FindAngles(T)
    rotationsList=[]
    for q in range(N-1,-1,-1):
        for p in range(q-1,-1,-1):
            A= np.matrix(np.identity(N))
            theta=anglesList[q][p]
            A[p,p]=math.cos(theta)
            A[q,q]=math.cos(theta)
            A[p,q]=math.sin(theta)
            A[q,p]=-math.sin(theta)
            rotationsList.append(A)
    return rotationsList

def TestAnglesFunctions(numSamples,N): #Generates numSamples elements of SO(N), and tests FindMatrix(FindAngles(T))=T
    for i in range(numSamples):
        T = ortho_group.rvs(N)
        T = np.matrix(np.diag([npl.det(T)]+[1]*(N-1)))*T #Converts an element of O(N) to SO(N)
        assert (np.max(abs(FindMatrix(FindAngles(T))-T))<1e-10)

def TestHoffmanDecomposition(numSamples,N): #Generates numSamples elements of SO(N), and tests that the Hoffman decomposition of T multiplies to T
    for i in range(numSamples):
        T = ortho_group.rvs(N)
        T = np.matrix(np.diag([npl.det(T)]+[1]*(N-1)))*T #Converts an element of O(N) to SO(N)
        TNew = np.identity(N)
        for A in FindHoffmanDecomposition(T):
            TNew = TNew*A
        assert (np.max(abs(TNew-T))<1e-10)

TestAnglesFunctions(100,5)
TestHoffmanDecomposition(100,5)
