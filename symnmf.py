import sys
import pandas as pd
import numpy as np
import mysymnmfsp as msm
np.random.seed(0)

def defX(filename):
    return pd.read_csv(filename, header=None, delim_whitespace=True).astype(float).values.tolist()

# Hinitialization
'''
def initial_H(W, n, k):
    np.random.seed(1234) 
    W_np = np.array(W)
    m = np.mean(W_np)
    upper_bound = 2 * np.sqrt(m / k)
    H0_np = np.random.uniform(0, upper_bound, size=(n, k))
    H0 = H0_np.tolist()
    return H0

'''
# Use np.random.uniform() for random selection

def initial_H(W, n, k):
    W_np = np.array(W)
    m = np.mean(W_np)
    upper_bound = 2 * np.sqrt(m / k)
    H0_np = np.random.uniform(low=0, high=upper_bound, size=(n, k))
    H0 = H0_np.tolist()
    return H0





def printMatrix(M):
    for vector in M:
        ret_val = ""
        for cord in vector:
            ret_val = ret_val + str("{:.4f}".format(cord)) +","
        ret_val = ret_val[:len(ret_val)-1]
        print(ret_val)
    return None


if len(sys.argv) == 3:
    goal = sys.argv[1]
    inputfile = sys.argv[2]
else:
    k = int(sys.argv[1])
    goal = sys.argv[2]
    inputfile = sys.argv[3]

X = defX(inputfile)
dim = len(X[0])
N = len(X)

if goal == "symnmf":
    W = msm.norm(X, N, dim)
    H0 = initial_H(W, N, k)
    H = msm.symnmf(W, H0, k, N, dim)
    printMatrix(H)

elif goal == "sym":
    A = msm.sym(X, N, dim)
    printMatrix(A)

elif goal == "ddg":
    D = msm.ddg(X, N, dim)
    printMatrix(D)

else:  # goal == "norm"
    W = msm.norm(X, N, dim)
    printMatrix(W)
