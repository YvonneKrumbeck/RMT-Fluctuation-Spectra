#!/usr/bin/env python
# coding: utf-8

# # Computing the Power Spectrum for a Large Ecosystem

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# ### Functions

# In[2]:


def get_A(n, pc, p, A_ii, sigma):
    '''Generates a community matrix for `n` species. 
    
    Parameters
    ----------
    n : int
        Number of species. 
    pc : float
        Probability of connection between species i and j.
    p : float
        Proportion of predator-prey interactions.
    A_ii : float
        Diagonal elements.
    sigma : float
        The interaction coefficients are drawn from a bivariate normal-distribution 
        with variance :math:`\sigma^2`, and zero mean. The covariance is determined by 
        the propotion of predator-prey interactions.
    '''
    
    gamma = 1-2*p  # interaction symmetry parameter
    
    # Initialise diagonal matrix
    A = A_ii *np.eye(n)
    
    # Set off-diagonal elements
    for i in range(n):
        for j in range(i+1,n):
            
            # Check if there is a connection between species i and j
            rand = np.random.rand()
            if rand < pc:
                    
                # Draw random interaction coefficients
                cov = [[sigma**2, sigma**2*gamma], [sigma**2*gamma, sigma**2]]
                A[i,j], A[j,i] = np.random.multivariate_normal([0,0],cov)
                    
    return A

def get_B(n, B_ii, A):
    '''Generates the noise matrix for a given community matrix `A`. 
    
    Parameters
    ----------
    n : int
        Number of species. 
    B_ii : float
        Diagonal elements.
    A : ndarray
        The `n` by `n` community matrix.        
    '''
    
    # Initialise diagonal matrix
    B = B_ii *np.eye(n)
    
    # Set B_ij=-|A_ij|, if predator-prey interaction
    AA = A*A.transpose()
    B[AA<0] += -np.abs(A[AA<0])
    
    return B    

def get_expectations(A, B):
    '''Returns necessary quantities from the community and noise matrices. 
    
    Parameters
    ----------
    A : ndarray
        Community matrix.
    B : ndarray
        Noise matrix.   
    '''
    
    # Compute expectations over non-zero off-diagonal elements
    A0 = np.copy(A)  # avoids in-place issues
    B0 = np.copy(B)
    np.fill_diagonal(A0, 0)
    np.fill_diagonal(B0, 0)
    
    EA2 = (A0**2)[A0!=0].mean()
    EAA = (A0*A0.transpose())[A0!=0].mean()
    EAB = (A0*B0)[A0!=0].mean()
    
    return EA2, EAA, EAB    

#################################################

def integrate(A, B, nt=5*10**5, dt=0.001, v=1000):
    '''Returns a timeseries of an Ornstein-Uhlenbeck process.
    
    Parameters
    ----------
    A : array_like
        Community matrix.
    B : array_like
        Noise matrix.
    nt : int, optional
        Number of time steps.
    dt : float, optional
        Time step size.
    v : int, optional
        System size (scales the noise).
    '''
    
    M = np.linalg.cholesky(B)
    
    # Initialise empty timeseries
    n = len(A)  # number of species
    x = np.zeros((nt,n))

    # Integration of linearised Ornstein-Uhlenbeck
    for ti in range(1,nt):
        determinitic = np.dot(A,x[ti-1]) *dt
        stochastic = 1/np.sqrt(v) *np.dot(M,np.random.normal(size=n)) *np.sqrt(dt)
        x[ti] = x[ti-1] +determinitic +stochastic

        # Terminate if population diverges
        if (x[ti] > 10**5).any():
            break

    return x 

def sample_timeseries(x, sample_rate, nt0=5*10**5, dt0=0.001):
    '''Returns a sampled timeseries and new number of time steps `nt`, step size `dt`.
    
    Parameters
    ----------
    x : array_like
        Original timeseries.
    sample_rate : int
        Step size for slicing the array `x`. 
    nt0 : int, optional
        Original number of time steps.
    dt : float, optional
        Original time step size.
    '''
    
    # Sample timeseries (discarding transient)
    x = x[::sample_rate]
    
    # Compute new number of time steps and step size
    nt = int(nt0/sample_rate)
    dt = dt0 *sample_rate
    
    return x, nt, dt

############################################

def powerspec_simulation(x, nt, dt):
    '''Returns the power spectral density computed from a timeseries `x`.
    
    Parameters
    ----------
    x : ndarray
        Timeseries of shape `(nt,n)`, where `n` is the number of species.
    nt : int
        Number of time steps. 
    dt : int
        Time step size.
    '''
    
    xf = x.transpose()  # .shape=(n,nt)
    
    # Center timeseries at zero
    xf = xf -xf.mean()
    
    # Fast Fourier Transform
    xf = np.fft.fft(xf, axis=-1)
    
    # Compute power spectrum
    ps = np.abs(xf)**2 *dt/nt
    
    # Store only values for positive frequencies
    ps = ps[:,:nt//2]
    
    # Compute corresponding frequencies
    w = np.fft.fftfreq(nt, dt)[:nt//2] *2*np.pi
    
    return w, ps.transpose()

def powerspec_direct(w_lst, A, B):
    '''Returns the exact solution for the power spectral density.
    
    Parameters
    ----------
    w_lst : array_like
        A list of the frequencies where the power spectrum is evaluated.
    A : ndarray
        Community matrix. 
    B : ndarray
        Noise matrix.
    '''
    
    n = len(A)  # number of species
    
    # Initialise empty list for power spectrum
    ps = []
    
    # Compute power spectrum
    for w in w_lst:
        A1 = A -1j*w *np.eye(n)
        A1 = np.linalg.inv(A1)
        A2 = A1.transpose().conjugate()
        
        ps_tmp = np.matmul(np.matmul(A1,B),A2).real
        
        # Store diagonal elements
        ps.append(np.diag(ps_tmp))
        
    return np.array(ps)

def powerspec_meanfield(w, A_ii, B_ii, EA2, EAA, EAB):
    '''Returns the mean-field approximation for the power spectral density.
    
    Parameters
    ----------
    w : ndarray
        Array of frequencies.
    A_ii : float
        Diagonal elements of the communtiy matrix.
    B_ii : float
        Diagonal elements of the noise matrix.
    EA2 : float
        Expected value of :math:`A_ij^2`.
    EAA : float
        Expected value of :math:`A_ijA_ji`.
    EAB : float
        Expacted value of :math:`A_ijB_ij`.    
    '''
    
    # Compute resolvent
    tmp = -A_ii +1j*w    
    if EAA != 0:
        r = 1/2/c/EAA *(tmp -np.sqrt(tmp**2 -4*c*EAA))
    else:
        r = 1/tmp
    
    # Compute power spectrum
    ps = (r.real**2 +r.imag**2) *(2*r.real*c*EAB +B_ii) /(1 -(r.real**2 +r.imag**2)*c*EA2)
    
    return r, ps    

def powerspec_singledefect(w_lst, A, B, r_mf, ps_mf):
    '''Returns the single defect approximation for the power spectral density.
    
    Parameters
    ----------
    w_lst : array_like
        List of frequencies.
    A : ndarray
        Community matrix.
    B : ndarray
        Noise matrix.
    r_mf : array_like
        Resolvent from mean-field approximation.
    ps_mf : array_like
        Power spectral density from mean-field approximation.  
    '''
    
    A_ii = A[0,0]
    B_ii = B[0,0]
    
    # Compute relevant sums over off-diagonal elements
    A0 = np.copy(A)  # avoids in-place issues
    B0 = np.copy(B)
    np.fill_diagonal(A0, 0)
    np.fill_diagonal(B0, 0)
    sumA2 = (A0*A0).sum(axis=1)
    sumAA = (A0*A0.transpose()).sum(axis=1)
    sumAB = (A0*B0).sum(axis=1)
    
    # Initialise array for power spectrum
    ps = np.zeros((len(w_lst), n))
    
    # Compute power spectrum
    for iw in range(len(w_lst)):
        w = w_lst[iw]
        denom = 1/np.abs(-A_ii +1j*w -r_mf[iw] *sumAA)**2
        ps[iw] = (ps_mf[iw]*sumA2 +2*r_mf[iw].real*sumAB +B_ii) *denom
    
    return ps


# ### Parameters

# In[6]:


n = 250
pc = 0.1
c = n*pc

p = 1
gamma = 1-2*p

b = 0.5
sigma = 0.25
mu = sigma *np.sqrt(2/np.pi)

A_ii = -b
B_ii = 2*b +c*mu
B_ii


# ### Generate community and noise matrices

# In[7]:


A = get_A(n, pc, p, A_ii, sigma)
B = get_B(n, B_ii, A)

EA2, EAA, EAB = get_expectations(A, B)
(EA2, sigma**2), (EAA, gamma*sigma**2), (EAB, 0)


# ### Simulate time series

# In[8]:


x = integrate(A, B)


# ### Compute power spectra

# In[9]:


w, ps_sim = powerspec_simulation(*sample_timeseries(x, sample_rate=500))
ps_dir = powerspec_direct(w, A, B)
r_mf, ps_mf = powerspec_meanfield(w, A_ii, B_ii, EA2, EAA, EAB)
ps_sda = powerspec_singledefect(w, A, B, r_mf, ps_mf)


# ### Plot

# In[10]:


fig, ax = plt.subplots(1,2, figsize=(12,4))
ax[0].set_title(r'Mean-Field Approximation')
ax[1].set_title(r'Single Defect Approximation')

#ax[0].plot(w, 1000*ps_sim, ':', color='gray', alpha=0.1)
ax[0].plot(w, 1000*ps_sim.mean(axis=1), ':', color='black', label=r'Simulation')
ax[0].plot(w, ps_dir.mean(axis=1), '.', color='orange', label=r'Direct');
ax[0].plot(w, ps_mf, '-', color='blue', label=r'Theory');

ax[1].plot(w, ps_dir, '.', color='orange', alpha=0.2);
ax[1].plot(w, ps_dir[:,0], '.', color='orange', alpha=1, label=r'Direct');
ax[1].plot(w, ps_sda, '-', color='blue', alpha=0.2);
ax[1].plot(w, ps_sda[:,0], '-', color='blue', alpha=1, label=r'Theory');

for a in ax:
    a.set_xlabel(r'frequency $\omega$')
    a.set_ylabel(r'power spectrum $\phi$')
    a.legend()
    #a.set_yscale('log')
    #a.set_ylim(0,10)

