# right dispersion

import math
import numpy as np
import numba as nb
import cmath

pi=math.pi

@nb.njit
def AS_core(g,nu,X,Y,nmedia,c,z):

    H = np.zeros(shape=(g.shape[0],g.shape[1],g.shape[2]), dtype=nb.types.complex64)
    U = np.zeros(shape=(g.shape[0],g.shape[1],g.shape[2]), dtype=nb.types.float64)

    fx=np.linspace(-g.shape[1]/(2*X),g.shape[1]/(2*X)-1/X,g.shape[1])
    fy=np.linspace(-g.shape[2]/(2*Y),g.shape[2]/(2*Y)-1/Y,g.shape[2])

    # transfer function H
    for k in range(g.shape[0]):
        for i in range(g.shape[1]):
            for j in range(g.shape[2]):
                U[k,i,j] = 1 - ((c * fx[i]/nu[k]/nmedia[k]) ** 2) - ((c * fy[j]/nu[k]/nmedia[k]) ** 2 )
                H[k,i,j] = np.exp(-1j*z*(2*pi*nu[k]*nmedia[k]/c)*np.sqrt(U[k,i,j]))

    return U,H


@nb.njit
def ASz_core(g,nu,X,Y,nmedia,c,z):

    Hz = np.zeros(shape=(g.shape[0],g.shape[1],g.shape[2]), dtype=nb.types.complex64)
    U = np.zeros(shape=(g.shape[0],g.shape[1],g.shape[2]), dtype=nb.types.float64)

    fx=np.linspace(-g.shape[1]/(2*X),g.shape[1]/(2*X)-1/X,g.shape[1])
    fy=np.linspace(-g.shape[2]/(2*Y),g.shape[2]/(2*Y)-1/Y,g.shape[2])

    # transfer function H
    for k in range(g.shape[0]):
        for i in range(g.shape[1]):
            for j in range(g.shape[2]):
                U[k,i,j] = 1 - ((c * fx[i]/nu[k]/nmedia[k]) ** 2) - ((c * fy[j]/nu[k]/nmedia[k]) ** 2 )
                Hz[k,i,j] = np.exp(-1j*z*(2*pi*nu[k]*nmedia[k]/c)*np.sqrt(U[k,i,j]))*fx[i]*c/(nu[k]*nmedia[k]*np.sqrt(U[k,i,j]))

    return U,Hz


@nb.njit
def ASz_polar_core(gx,gy,nu,X,Y,nmedia,c,z):

    Hxz = np.zeros(shape=(gx.shape[0],gx.shape[1],gx.shape[2]), dtype=nb.types.complex64)
    Hyz = np.zeros(shape=(gx.shape[0],gx.shape[1],gx.shape[2]), dtype=nb.types.complex64)
    U = np.zeros(shape=(gx.shape[0],gx.shape[1],gx.shape[2]), dtype=nb.types.float64)

    fx=np.linspace(-gx.shape[1]/(2*X),gx.shape[1]/(2*X)-1/X,gx.shape[1])
    fy=np.linspace(-gx.shape[2]/(2*Y),gx.shape[2]/(2*Y)-1/Y,gx.shape[2])

    # transfer function H
    for k in range(gx.shape[0]):
        for i in range(gx.shape[1]):
            for j in range(gx.shape[2]):
                U[k,i,j] = 1 - ((c * fx[i]/nu[k]/nmedia[k]) ** 2) - ((c * fy[j]/nu[k]/nmedia[k]) ** 2 )
                Hxz[k,i,j] = np.exp(-1j*z*(2*pi*nu[k]*nmedia[k]/c)*np.sqrt(U[k,i,j]))*fx[i]*c/(nu[k]*nmedia[k]*np.sqrt(U[k,i,j]))
                Hyz[k,i,j] = np.exp(-1j*z*(2*pi*nu[k]*nmedia[k]/c)*np.sqrt(U[k,i,j]))*fy[j]*c/(nu[k]*nmedia[k]*np.sqrt(U[k,i,j]))

    return U,Hxz,Hyz


@nb.njit
def ASeva_core(g,nu,X,Y,nmedia,c,z):

    Heva = np.zeros(shape=(g.shape[0],g.shape[1],g.shape[2]), dtype=nb.types.complex64)
    U = np.zeros(shape=(g.shape[0],g.shape[1],g.shape[2]), dtype=nb.types.float64)

    fx=np.linspace(-g.shape[1]/(2*X),g.shape[1]/(2*X)-1/X,g.shape[1])
    fy=np.linspace(-g.shape[2]/(2*Y),g.shape[2]/(2*Y)-1/Y,g.shape[2])

    # transfer function H
    for k in range(g.shape[0]):
        for i in range(g.shape[1]):
            for j in range(g.shape[2]):
                U[k,i,j] = 1 - ((c * fx[i]/nu[k]/nmedia[k]) ** 2) - ((c * fy[j]/nu[k]/nmedia[k]) ** 2 )
                if U[k,i,j]>=0:
                    Heva[k,i,j]=0
                else:
                    Heva[k,i,j] = cmath.exp(-1j*z*(2*pi*nu[k]*nmedia[k]/c)*cmath.sqrt(U[k,i,j]))*(-1)

    return U,Heva


@nb.njit
def C_core(g,nu,X,Y,nmedia,c,z):

    h = np.zeros(shape=(g.shape[0],g.shape[1],g.shape[2]), dtype=nb.types.complex64)

    Nx=g.shape[1]
    Ny=g.shape[2]
    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    
    # transfer function H
    for k in range(g.shape[0]):
        for i in range(g.shape[1]):
            for j in range(g.shape[2]):
                r=np.sqrt((x[i]*x[i])+(y[j]*y[j])+(z*z))
                h[k,i,j]=((X/Nx)**2)*np.exp(-1j*r*2*pi*nu[k]*nmedia[k]/c)*z*nu[k]*nmedia[k]/(-1j*c*r*r)

    return h


@nb.njit
def Cz_core(g,nu,X,Y,nmedia,c,z):

    U = np.zeros(shape=(g.shape[0],g.shape[1],g.shape[2]), dtype=nb.types.float64)
    h = np.zeros(shape=(g.shape[0],g.shape[1],g.shape[2]), dtype=nb.types.complex64)
    kxkz = np.zeros(shape=(g.shape[0],g.shape[1],g.shape[2]), dtype=nb.types.float64)

    Nx=g.shape[1]
    Ny=g.shape[2]
    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    fx=np.linspace(-Nx/(2*X),Nx/(2*X)-1/X,Nx)
    fy=np.linspace(-Ny/(2*Y),Ny/(2*Y)-1/Y,Ny)

    # transfer function H
    for k in range(g.shape[0]):
        for i in range(g.shape[1]):
            for j in range(g.shape[2]):
                r=np.sqrt((x[i]*x[i])+(y[j]*y[j])+(z*z))
                U[k,i,j] = 1 - ((c * fx[i]/nu[k]/nmedia[k]) ** 2) - ((c * fy[j]/nu[k]/nmedia[k]) ** 2 )
                h[k,i,j]=((X/Nx)**2)*np.exp(-1j*r*2*pi*nu[k]*nmedia[k]/c)*z*nu[k]*nmedia[k]/(-1j*c*r*r)
                kxkz[k,i,j]= fx[i]*c/(nu[k]*nmedia[k]*np.sqrt(U[k,i,j]))

    return U,h,kxkz


@nb.njit
def Cz_polar_core(gx,gy,nu,X,Y,nmedia,c,z):

    U = np.zeros(shape=(gx.shape[0],gx.shape[1],gx.shape[2]), dtype=nb.types.float64)
    h = np.zeros(shape=(gx.shape[0],gx.shape[1],gx.shape[2]), dtype=nb.types.complex64)
    kxkz = np.zeros(shape=(gx.shape[0],gx.shape[1],gx.shape[2]), dtype=nb.types.float64)
    kykz = np.zeros(shape=(gx.shape[0],gx.shape[1],gx.shape[2]), dtype=nb.types.float64)

    Nx=gx.shape[1]
    Ny=gx.shape[2]
    x=np.linspace(-X/2,X/2-X/Nx,Nx)
    y=np.linspace(-Y/2,Y/2-Y/Ny,Ny)
    fx=np.linspace(-Nx/(2*X),Nx/(2*X)-1/X,Nx)
    fy=np.linspace(-Ny/(2*Y),Ny/(2*Y)-1/Y,Ny)

    # transfer function H
    for k in range(gx.shape[0]):
        for i in range(gx.shape[1]):
            for j in range(gx.shape[2]):
                r=np.sqrt((x[i]*x[i])+(y[j]*y[j])+(z*z*nmedia[k]*nmedia[k]))
                U[k,i,j] = 1 - ((c * fx[i]/nu[k]/nmedia[k]) ** 2) - ((c * fy[j]/nu[k]/nmedia[k]) ** 2 )
                h[k,i,j]=((X/Nx)**2)*np.exp(-1j*r*2*pi*nu[k]*nmedia[k]/c)*z*nu[k]*nmedia[k]/(-1j*c*r*r)
                kxkz[k,i,j]= fx[i]*c/(nu[k]*nmedia[k]*np.sqrt(U[k,i,j]))
                kykz[k,i,j]= fy[j]*c/(nu[k]*nmedia[k]*np.sqrt(U[k,i,j]))

    return U,h,kxkz,kykz





