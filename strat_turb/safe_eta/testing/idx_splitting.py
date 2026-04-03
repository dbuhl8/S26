import numpy as np
import matplotlib.pyplot as plt

k, n, m = (2,10,10)
A = np.zeros((k,n,m))

i = 0 

A[i,2,1] = 1
A[i,1,5] = 1
A[i,3,7] = 1
A[i,6,8] = 1

idx = np.where(A[i,:,:] == 2)
print(idx)
print(A[i,*idx].mean())
print(len(idx[0])/100)


