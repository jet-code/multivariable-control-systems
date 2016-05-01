
# coding: utf-8

# In[1]:

# Alexander Hebert
# ECE 6390
# Computer Project #2
# Method 2/3


# In[2]:

# Tested using Python v3.4 and IPython v2


##### Import libraries and functions

# In[3]:

import numpy as np


# In[4]:

from PlaneRotationFn import planeRotation1
from PlaneRotationFn import planeRotation2


# In[5]:

import scipy


# In[6]:

import sympy


# In[7]:

from IPython.display import display


# In[8]:

from sympy.interactive import printing


# In[9]:

np.set_printoptions(precision=6)


# In[10]:

#np.set_printoptions(suppress=True)


##### Original system:

# In[11]:

A = np.loadtxt('A_ex1.txt')


# In[12]:

A


# In[13]:

n,nc = A.shape


# In[14]:

B = np.loadtxt('B_ex1.txt')


# In[15]:

B


# In[16]:

nr,m = B.shape


##### Compute eigenvalues/poles of A to determine system stability:

# In[17]:

A_eigvals, M = np.linalg.eig(A)


# In[18]:

A_eigvals


# In[19]:

# Two poles lie in the RHP and are unstable.


# In[20]:

A_eigvals_desired = np.array([-0.2,-0.5,A_eigvals[2],A_eigvals[3]])


# In[21]:

A_eigvals_desired


# In[22]:

Lambda = np.diag(A_eigvals_desired)


# In[23]:

Lambda


##### Pole Assignment Algorithm from journal paper

# In[24]:

# Step A: Decomposition of B using SVD
# B = U*S*V.H


# In[25]:

U, s, VH = np.linalg.svd(B)


# In[26]:

U


# In[27]:

s


# In[28]:

S = np.zeros((4, 2))
S[:2, :2] = np.diag(s)


# In[29]:

S


# In[30]:

VH


# In[31]:

# Extract U_0 and U_1 from matrix U = [U_0,U_1]


# In[32]:

U_0 = U[:n,:m]


# In[33]:

U_0


# In[34]:

U_1 = U[:n,m:]


# In[35]:

U_1


# In[36]:

# B = [U_0,U_1][Z,0].T
# Compute Z from SVD of B


# In[37]:

Z = np.diag(s).dot(VH)


# In[38]:

Z


# In[39]:

# U_1.T *(A - lambda_j*I)
# Compute S_hat_j and S_j

for j in range(len(A_eigvals_desired)):
    
    lambda_j = A_eigvals_desired[j]
    
    # M_j is a temp matrix
    exec("M_%d = np.dot(U_1.T,(A - lambda_j*np.identity(n)))" %(j+1))
    
    # U_1.T *(A - lambda_j*I) = T_j *[Gamma_j,0]*[S_j_hat,S_j].T
    exec("T_%d, gamma_%d, SH_%d = np.linalg.svd(M_%d)" %(j+1,j+1,j+1,j+1))
    
    exec("S_hat_%d = SH_%d[:m,:].T" %(j+1,j+1))
    exec("S_%d = SH_%d[m:,:].T" %(j+1,j+1))
    


# In[40]:

# Initial eigenvectors in X_tilde
X_tilde = np.eye(n)
X_tilde


# In[41]:

# Initial eigenvectors in X
X = np.zeros((n,n))
X


# In[42]:

# Step X with Method 2/3

maxiter = 1
v4prev = 0
v4current = 0

for r in range(maxiter):
    
    for j in range(n):
        
        xt_j = X_tilde[:,j].reshape((n,1))
        
        for k in range(j+1,n,1):
            
            xt_k = X_tilde[:,k].reshape((n,1))
            
            exec("phi_j = np.linalg.norm(np.dot(S_hat_%d.T,xt_j))" %(j+1))
            exec("phi_k = np.linalg.norm(np.dot(S_hat_%d.T,xt_k))" %(k+1))
            
            if (phi_j < phi_k):
                
                sin_theta = phi_j
                cos_theta = np.sqrt(1 - phi_j**2)
                protation = planeRotation2(n,cos_theta,sin_theta,j,k)
                v4current = v4current + phi_j**2
                
            else:
                
                sin_theta = phi_k
                cos_theta = np.sqrt(1 - phi_k**2)
                protation = planeRotation2(n,cos_theta,sin_theta,j,k)
                v4current = v4current + phi_k**2
            
            X_tilde = np.dot(protation,X_tilde)
            
    v4current = np.sqrt(v4current)
    print(v4current - v4prev)
    v4prev = v4current
    


# In[43]:

# Compute eigenvectors x_j in X from X_tilde

for j in range(n):
    
    xt_j = X_tilde[:,j].reshape((n,1)) 
    exec("x_j = np.dot(np.dot(S_%d,S_%d.T),xt_j) / np.linalg.norm(np.dot(S_%d.T,xt_j))" %(j+1,j+1,j+1))
    X[:,j] = x_j[:,0]
    


# In[44]:

X


# In[45]:

np.linalg.matrix_rank(X)


# In[46]:

X_inv = np.linalg.inv(X)


# In[47]:

X_inv


# In[48]:

# M defined as A + BF
M = X.dot(Lambda).dot(X_inv)


# In[49]:

M


# In[50]:

# Compute feedback matrix F
F = np.dot(np.linalg.inv(Z),np.dot(U_0.T,(M - A)))


# In[51]:

F


# In[52]:



