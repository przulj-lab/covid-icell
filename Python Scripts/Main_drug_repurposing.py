# -*- coding: UTF-8 -*-

""" Non-negative matrix tri-factorization (numpy)"""
# Author: NMD

import Network_Matrices as nm
import Matrix_Factorization as mf
import numpy as np
import math

#Load matrix factor G from iCell
case='Patient_Infected'

G, genes, genes_features  = nm.Load_Matrix_Factor("./Data/G matrix factors/%s_B.csv"%(case))
gene2ind = {}
for i in range(len(genes)):
	gene2ind[genes[i]] = i

n,k1=np.shape(G)
features1 = [str(i) for i in range(k1)]

#Load Chemical similarity, and make it's Laplacian
Chem, drugs, drug2ind = nm.Load_X_to_X('./Data/Drug related data/Drugs_similarities_top5percent.csv', ',')
Lchem = nm.Make_Laplacian_Matrix(Chem)

#Load DTI
DTI = nm.Load_DTI('./Data/Drug related data/Drug_Target_Interactions.csv', gene2ind, drug2ind)
_,m = np.shape(DTI)
k2 = int(math.sqrt(m/2.))

features2 = [str(i) for i in range(k2)]

#influence of graph regularization, should be = 1., 0 for testing
gamma = 1.

#Run PB like data integration to repurpose drugs
Solver = mf.PSBLike_Drug_NMTF(max_iter=500, verbose = 10)
#D, S = Solver.Solve_MUR(DTI, Lchem, G, k1, k2, gamma, init="rand")
D, S = Solver.Solve_MUR(DTI, Lchem, G, k1, k2, gamma, init="SVD")

nm.Save_Matrix_Factor(D, "%s_D.csv"%(case), drugs, features2)
nm.ave_Matrix_Factor(S, "%s_S.csv"%(case), features1, features2)
