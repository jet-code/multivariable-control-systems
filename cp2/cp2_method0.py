
# coding: utf-8

# In[1]:

# Alexander Hebert
# ECE 6390
# Computer Project #2


# In[2]:

# Tested using Python v3.4 and IPython v2


##### Import libraries

# In[3]:

import numpy as np


# In[4]:

import scipy


# In[5]:

import sympy


# In[6]:

from IPython.display import display


# In[7]:

from sympy.interactive import printing


# In[8]:

np.set_printoptions(precision=6)


# In[9]:

#np.set_printoptions(suppress=True)


##### Original system:

# In[10]:

A = np.loadtxt('A_ex1.txt')


# In[11]:

A


# In[12]:

n,nc = A.shape


# In[13]:

B = np.loadtxt('B_ex1.txt')


# In[14]:

B


# In[15]:

nr,m = B.shape


##### Compute eigenvalues/poles of A to determine system stability:

# In[16]:

A_eigvals, M = np.linalg.eig(A)


# In[17]:

A_eigvals


# In[18]:

# Two poles lie in the RHP and are unstable.


# In[19]:

A_eigvals_desired = np.array([-0.2,-0.5,A_eigvals[2],A_eigvals[3]])


# In[20]:

A_eigvals_desired


# In[21]:

Lambda = np.diag(A_eigvals_desired)


# In[22]:

Lambda


##### Pole Assignment Algorithm from journal paper

# In[23]:

# Step A: Decomposition of B using SVD
# B = U*S*V.H


# In[24]:

U, s, VH = np.linalg.svd(B)


# In[25]:

U


# In[26]:

s


# In[27]:

S = np.zeros((4, 2))
S[:2, :2] = np.diag(s)


# In[28]:

S


# In[29]:

VH


# In[30]:

# Extract U_0 and U_1 from matrix U = [U_0,U_1]


# In[31]:

U_0 = U[:n,:m]


# In[32]:

U_0


# In[33]:

U_1 = U[:n,m:]


# In[34]:

U_1


# In[35]:

# B = [U_0,U_1][Z,0].T
# Compute Z from SVD of B


# In[36]:

Z = np.diag(s).dot(VH)


# In[37]:

Z


# In[38]:

# Compute the nullspace of U_1.T *(A - lambda_j*I)
# for initial eigenvectors in X

X = np.zeros((n,n))

for j in range(len(A_eigvals_desired)):
    
    lambda_j = A_eigvals_desired[j]
    
    # M_j is a temp matrix
    exec("M_%d = np.dot(U_1.T,(A - lambda_j*np.identity(n)))" %(j+1))
    
    # U_1.T *(A - lambda_j*I) = T_j *[Gamma_j,0]*[S_j_hat,S_j].T
    exec("T_%d, gamma_%d, SH_%d = np.linalg.svd(M_%d)" %(j+1,j+1,j+1,j+1))
    
    exec("X[:,j] = SH_%d[-2,:]" %(j+1))
    # no transpose in SH_j due to 1-d vector
    
    exec("S_hat_%d = SH_%d[:m,:].T" %(j+1,j+1))
    exec("S_%d = SH_%d[m:,:].T" %(j+1,j+1))
    


# In[39]:

# Initial eigenvectors in X
X


# In[40]:

# Test X for full rank
X_rank = np.linalg.matrix_rank(X)


# In[41]:

all((X_rank,n))


# In[42]:

# Step X with Method 0

maxiter = 2
v2current = 0
v2prev = np.linalg.cond(X)
eps = 10e-5
flag = 0
X_j = np.zeros((n,n-1))
cond_num = np.zeros((n,1))

for r in range(maxiter):
    
    for j in range(n):
        
        X_j = np.delete(X,j,1)
        Q,R = np.linalg.qr(X_j,mode='complete')
        y_j = Q[:,-1].reshape((4,1))
        
        exec("S_j = S_%d" %(j+1))
        x_j = (S_j.dot(S_j.T).dot(y_j) / np.linalg.norm(np.dot(S_j.T,y_j)))
        X[:,j] = x_j[:,0]
        
        cond_num[j,0] = 1 / np.abs(np.dot(y_j.T,x_j))
        
    v2current = np.linalg.cond(X)
    
    if ((v2current - v2prev) < eps):
        print("Tolerance met")
        print("v2 = %.3f" %v2current)
        flag = 1
        
    else:
        v2prev = v2current
    
if (flag == 0):
    print("Tolerance not met")
    print("v2 = %.3f" %v2current)
    


# In[43]:

X


# In[44]:

np.linalg.matrix_rank(X)


# In[45]:

X_inv = np.linalg.inv(X)


# In[46]:

X_inv


# In[47]:

# M defined as A + BF
M = X.dot(Lambda).dot(X_inv)


# In[48]:

M


# In[49]:

# Eigenvalues of controlled system
M_eigvals, H = np.linalg.eig(M)
M_eigvals


# In[50]:

# Compute feedback matrix F
F = np.dot(np.linalg.inv(Z),np.dot(U_0.T,(M - A)))


# In[51]:

F


# In[52]:

np.linalg.norm(F)


# In[53]:

# Compute condition number norms


# In[54]:

# Inf norm
np.linalg.norm(cond_num,np.inf)


# In[55]:

# 2 norm
np.linalg.norm(cond_num)


# In[55]:



