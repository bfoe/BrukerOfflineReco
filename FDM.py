# The following implements the Finite Difference Method (FDM) denoising 
# originally described in:
# Song SM, Napel S, Glover GH, Pelc NJ.
# "Noise reduction in three-dimensional phase-contrast MR velocity measurements"
# J Magn Reson Imaging. 1993 Jul-Aug;3(4):587-96. doi: 10.1002/jmri.1880030407
#
# This Python port is based on Frank Ong's 2013 Mathlab code, found at 
# http://people.eecs.berkeley.edu/~mlustig/software/DivFreeWavelet_code.zip
# discussed (among various methods) in:
# Ong F, Uecker M, Tariq U, Hsiao A, Alley MT, Vasanawala SS, Lustig M.
# "Robust 4D Flow Denoising using Divergence-free Wavelet Transform"
# Magn Reson Med. 2015 Feb;73(2):828-42. doi: 10.1002/mrm.25176
#
# You are encouraged to modify/distribute this code. 
# Please acknowledge this code and cite the papers appropriately.
#

import numpy as np

# Calculates divergence
def div (vx,vy,vz,res):
    dx = (vx-np.roll(vx,1,axis=0))/res[0]
    dy = (vy-np.roll(vy,1,axis=1))/res[1]
    dz = (vz-np.roll(vz,1,axis=2))/res[2]
    return  dx+dy+dz

# calculates gradient (adjoint to div)    
def grad(inp,res):
    fx = (np.roll(inp,-1,axis=0)-inp)/res[0];
    fy = (np.roll(inp,-1,axis=1)-inp)/res[1];
    fz = (np.roll(inp,-1,axis=2)-inp)/res[2];
    return fx,fy,fz

# First order Poisson solver using the FFT with periodic extension
# Solves \laplacian u = f
def poisson(f,res):
    [nx,ny,nz] = f.shape
    [X,Y,Z] = np.mgrid[0:nx,0:ny,0:nz]
    # lambdas are Fourier Transform of Discrete Laplacian: [-2,1,0...0,1];
    lambdax= -4*np.square(np.sin(X*np.pi/nx))/res[0]**2 
    lambday= -4*np.square(np.sin(Y*np.pi/ny))/res[1]**2
    lambdaz= -4*np.square(np.sin(Z*np.pi/nz))/res[2]**2     
    mu = ( lambdax + lambday + lambdaz)
    mu[0,0,0] = 1   
    u = np.divide(np.fft.fftn(f),mu)    
    u[0,0,0]=0
    u = np.real(np.fft.ifftn(u))
    return u
    
# Finite difference method denoising 
# with Hodge decomposition using the fast Poisson solver
def fdmDenoise(vx,vy,vz,res):
    d = div(vx,vy,vz,res)           # d = div(v) = div(grad(u) = laplacian(u)
    u = poisson(d,res)              # u = laplacian^-1(u)
    [vxCF,vyCF,vzCF] = grad(u,res); # vCF = grad(u)
    vxDF = vx-vxCF;
    vyDF = vy-vyCF;
    vzDF = vz-vzCF;
    return vxDF, vyDF, vzDF
    

   
        
''' Python testcode below
      
ux = np.asarray([[[0.1,0.2,0.3],[0.4,0.5,0.6],[0.7,0.8,0.9]],\
                 [[1.1,1.2,1.3],[1.4,1.5,1.6],[1.7,1.8,1.9]],\
                 [[2.1,2.2,2.3],[2.4,2.5,2.6],[2.7,2.8,2.9]]])
                 
uy = np.asarray([[[3.1,3.2,3.3],[3.4,3.5,3.6],[3.7,3.8,3.9]],\
                 [[4.1,4.2,4.3],[4.4,4.5,4.6],[4.7,4.8,4.9]],\
                 [[5.1,5.2,5.3],[5.4,5.5,5.6],[5.7,5.8,5.9]]])
                 
uz = np.asarray([[[6.1,6.2,6.3],[6.4,6.5,6.6],[6.7,6.8,6.9]],\
                 [[7.1,7.2,7.3],[7.4,7.5,7.6],[7.7,7.8,7.9]],\
                 [[8.1,8.2,8.3],[8.4,8.5,8.6],[8.7,8.8,8.9]]])
                 
res = np.asarray([0.2,0.5,0.8])

#print ("ux[:,:,0] =");print(ux[:,:,0]);print
#print ("ux[:,:,1] =");print(ux[:,:,1]);print
#print ("ux[:,:,2] =");print(ux[:,:,2]);print
#print ("uy[:,:,0] =");print(uy[:,:,0]);print
#print ("uy[:,:,1] =");print(uy[:,:,1]);print
#print ("uy[:,:,2] =");print(uy[:,:,2]);print
#print ("uz[:,:,0] =");print(uz[:,:,0]);print
#print ("uz[:,:,1] =");print(uz[:,:,1]);print
#print ("uz[:,:,2] =");print(uz[:,:,2]);print

d = div (ux,uy,uz,res)
#print ("d[:,:,0] =");print(d[:,:,0]);print
#print ("d[:,:,1] =");print(d[:,:,1]);print
#print ("d[:,:,2] =");print(d[:,:,2]);print

fx,fy,fz = grad (d,res)
#print ("fx[:,:,0] =");print(fx[:,:,0]);print
#print ("fx[:,:,1] =");print(fx[:,:,1]);print
#print ("fx[:,:,2] =");print(fx[:,:,2]);print
#print ("uy[:,:,0] =");print(uy[:,:,0]);print
#print ("uy[:,:,1] =");print(uy[:,:,1]);print
#print ("uy[:,:,2] =");print(uy[:,:,2]);print
#print ("uz[:,:,0] =");print(uz[:,:,0]);print
#print ("uz[:,:,1] =");print(uz[:,:,1]);print
#print ("uz[:,:,2] =");print(uz[:,:,2]);print

p=poisson(ux,res)
#print ("p[:,:,0] =");print(p[:,:,0]);print
#print ("p[:,:,1] =");print(p[:,:,1]);print
#print ("p[:,:,2] =");print(p[:,:,2]);print

vxDF,vyDF,vzDF=fdmDenoise (ux,uy,uz,res)
#print ("vxDF[:,:,0] =");print(vxDF[:,:,0]);print
#print ("vxDF[:,:,1] =");print(vxDF[:,:,1]);print
#print ("vxDF[:,:,2] =");print(vxDF[:,:,2]);print
#print ("vyDF[:,:,0] =");print(vyDF[:,:,0]);print
#print ("vyDF[:,:,1] =");print(vyDF[:,:,1]);print
#print ("vyDF[:,:,2] =");print(vyDF[:,:,2]);print
#print ("vzDF[:,:,0] =");print(vzDF[:,:,0]);print
#print ("vzDF[:,:,1] =");print(vzDF[:,:,1]);print
#print ("vzDF[:,:,2] =");print(vzDF[:,:,2]);print
'''

''' Mathlab testcode below (extract to .m file)

ux = zeros(3, 3, 3);
ux(1,:,:) = [0.1,0.2,0.3; 0.4,0.5,0.6; 0.7,0.8,0.9];
ux(2,:,:) = [1.1,1.2,1.3; 1.4,1.5,1.6; 1.7,1.8,1.9];
ux(3,:,:) = [2.1,2.2,2.3; 2.4,2.5,2.6; 2.7,2.8,2.9];

uy = zeros(3, 3, 3);
uy(1,:,:) = [3.1,3.2,3.3; 3.4,3.5,3.6; 3.7,3.8,3.9];
uy(2,:,:) = [4.1,4.2,4.3; 4.4,4.5,4.6; 4.7,4.8,4.9];
uy(3,:,:) = [5.1,5.2,5.3; 5.4,5.5,5.6; 5.7,5.8,5.9];

uz = zeros(3, 3, 3);
uz(1,:,:) = [6.1,6.2,6.3; 6.4,6.5,6.6; 6.7,6.8,6.9];
uz(2,:,:) = [7.1,7.2,7.3; 7.4,7.5,7.6; 7.7,7.8,7.9];
uz(3,:,:) = [8.1,8.2,8.3; 8.4,8.5,8.6; 8.7,8.8,8.9];

res = [0.2,0.5,0.8];
d = div (ux,uy,uz,res);
[fx,fy,fz] = grad (d,res);
p=poisson(ux,res);
[vxDF,vyDF,vzDF]=fdmDenoise (ux,uy,uz,res)
'''