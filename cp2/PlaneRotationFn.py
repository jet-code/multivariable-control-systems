import numpy as np

def planeRotation1(size,theta,j,k):
    """Compute the rotation matrix 
       in the plane of columns j an k (j < k)
       by angle theta (in radians)."""
    
    protation = np.identity(size)
    protation[j,j] = np.cos(theta)
    protation[k,j] = np.sin(theta)
    protation[k,k] = np.cos(theta)
    protation[j,k] = -1*np.sin(theta)
    
    return protation
    
def planeRotation2(size,cos_theta,sin_theta,j,k):
    """Compute the rotation matrix 
       in the plane of columns j an k (j < k)
       by cos_theta and sin_theta."""
    
    protation = np.identity(size)
    protation[j,j] = cos_theta
    protation[k,j] = sin_theta
    protation[k,k] = cos_theta
    protation[j,k] = -1*sin_theta
    
    return protation