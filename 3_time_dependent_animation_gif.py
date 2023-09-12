# -*- coding: utf-8 -*-
"""
Created on Wed May 25 16:09:15 2022
@author: NanSong52
Q3: Analysing the time-dependent Schrödinger equation for the perpetual particle in the box.
"""

import imageio
import os
import numpy as np
import matplotlib.pyplot as plt

path=os.getcwd()
if not os.path.exists('mov'):
    os.makedirs('mov')
os.chdir(str(path)+'/mov')

def integral(f,a,b,steps):
    # f=integral function of x
    # a lower limit
    # b upper limit
    # steps:(steps-1) trapezoids
    L=b-a
    interval=L/(steps-1)
    cumulator=0
    for x in range(steps-1):
        f1 = f(x*interval)
        f2 = f((x+1)*interval)
        cumulator += 0.5 * (f1+f2) * interval
    return cumulator

 # How to use this function?
 # 1st, define parameters:
m=3 # parameter
N=1 # normalize factor
L=2 # length of 1-D box
 # 2nd, define integral equation


def Psi1(x):
    return N*(np.sin(np.pi*x/L) + np.sin(2*np.pi*x/L))

def Psi2(x):
    return N*np.sin(np.pi*x/L)
# density not change

def Psi3(x):
    return N*x**m*(L-x)
# different coefficients

def Psi3_sq(x):
    return Psi3(x)*Psi3(x)


#For Psi1
#(a) Nomalize these wavefunctions to get the normalized factor:
def Psi1_sq(x):
    return N*(np.sin(np.pi*x/L)+np.sin(2*np.pi*x/L)) * N*(np.sin(np.pi*x/L)+np.sin(2*np.pi*x/L))
a = integral(Psi1_sq,0,L,100)
print('The normalize factor of Psi1(x) is '+ str(1/np.sqrt(a)) + 
      '\n Now bring this factor back in to Psi1(i.e. change N and reload)')

#(b) Calculate the expansion coefficients c_n
def phi(n,x):
    return np.sqrt(2/L)*np.sin(n*np.pi*x/L)


# Time-dependent schordinger equation
Psi=Psi2 ############################################################  Change here
step=100 ############################################################  Integral Steps
M = 50   ############################################################  Change here M=50
List_cn_0=[]
List_cn_real=[]
List_cn_imaginary=[]

A=np.zeros([M,step])
interval=L/step

def cn_0(x):
    return phi(n,x)*Psi(x)

######################### plot #########################
# Plot the time dependent functions ΨR(x,τ),ΨI (x,τ),ρ(x,τ) for time τ =1.
# We can now plot the wave function pieces for different times, τ = 0.05, 0.1, 0.2, 0.3...
# we can create a dynamic video (gif) with different tau, and animate.
n_tau=50
seg=1/n_tau
tau_list=[]
for i in range(n_tau+1):
    tau_list.append(round(i*seg,2))

for tau in tau_list:
    List_cn_0=[]
    List_cn_real=[]
    List_cn_imaginary=[]

    for n in range(M):
        Cn_0 = integral(cn_0,0,L,step)
        Cn_real = np.cos(2*np.pi*n*n*tau)*Cn_0
        Cn_imaginary = (-1)*np.sin(2*np.pi*n*n*tau)*Cn_0
        List_cn_0.append(Cn_0)
        List_cn_real.append(Cn_real)
        List_cn_imaginary.append(Cn_imaginary)

        for x in range(step):
            A[n][x] = phi(n,x*interval)
            
    Psi_real = np.dot(List_cn_real, A)
    Psi_imaginary = np.dot(List_cn_imaginary, A)
    rho = Psi_real*Psi_real + Psi_imaginary*Psi_imaginary
    exact=np.dot(List_cn_0,A)
    
    # create list of x variable
    Lx=[]
    for x in range (step):
        Lx.append(x*interval)

    # save figures
    plt.figure(figsize=(6,8))
    plt.subplot(2,1,1)
    plt.title('Time Dependent Functions Psi(tau='+str(tau)+')')
    plt.xlim(0, 2)
    plt.ylim(-2, 2)
    #plt.plot(Lx,exact,label="Psi_((tau=0/1))")
    plt.plot(Lx,Psi_real,label="Psi_real")
    plt.plot(Lx,Psi_imaginary,label="Psi_imagin")
    plt.legend(loc=3)
    
    plt.subplot(2,1,2)
    plt.title('rho(tau='+str(tau)+')')
    plt.xlim(0, 2)
    plt.ylim(-2, 2)
    plt.plot(Lx,rho,label="rho")
    plt.legend(loc=3)
    plt.savefig('Psi_and_rho(tau='+str(tau)+').png')
os.chdir(path)


if str(Psi).split()[1]=='Psi3':
    filename = str(Psi).split()[1]+'_and_rho_m='+str(m)+'.gif'
else:
    filename = str(Psi).split()[1]+'_and_rho.gif'


with imageio.get_writer(filename, mode='I') as writer:
    imgs = os.listdir(path+'/mov')
    os.chdir(path+'/mov')
    for filename in imgs:
        image = imageio.imread(filename)
        writer.append_data(image)
os.chdir(path)

