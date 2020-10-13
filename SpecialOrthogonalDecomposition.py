import numpy as np
import math

def FindAngles(T): #Find a list of angles [[],[theta_]
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

def FindMatrix(anglesList):
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

def FindMatrixSimpleRotations(anglesList):
    N = len(anglesList)
    T= np.matrix(np.identity(N))
    for q in range(N-1,-1,-1):
        for p in range(q-1,-1,-1):
            A= np.matrix(np.identity(N))
            theta=anglesList[q][p]
            A[p,p]=math.cos(theta)
            A[q,q]=math.cos(theta)
            A[p,q]=math.sin(theta)
            A[q,p]=-math.sin(theta)
            T=T*A
    return T

