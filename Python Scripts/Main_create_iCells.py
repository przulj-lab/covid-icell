# -*- coding: UTF-8 -*-

""" Non-negative matrix tri-factorization (numpy)"""
# Author: NMD
import os
os.environ["OMP_NUM_THREADS"] = "4" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "4" # export OPENBLAS_NUM_THREADS=4
os.environ["MKL_NUM_THREADS"] = "4" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "4" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "4" # export NUMEXPR_NUM_THREADS=6

import sys
sys.path.append("./Python Scripts")

import networkx as nx
import numpy as np
import Network_Matrices as nm
import Matrix_Factorization as mf

print "\n\x1b[1;37;44m ############################################# \x1b[0m"
print "\x1b[1;37;44m #                                           # \x1b[0m"
print "\x1b[1;37;44m # ISF Framework                             # \x1b[0m"
print "\x1b[1;37;44m #                                           # \x1b[0m"
print "\x1b[1;37;44m ############################################# \x1b[0m"


cases = ['A549_Control', 'A549_Infected', 'NHBE_Control', 'NHBE_Infected','CALU_Control', 'CALU_Infected', 'Patient_Control', 'Patient_Infected']
tissue = cases[int(sys.argv[1])]

#Loading and preparing all the data

print "---- Loading PPI network"
PPI, nodes, node2ind = nm.Load_Network("./Data/Covid_Networks/%s_PPI.net"%(tissue))
PPI_mat = nm.Make_Adj(PPI, nodes, node2ind)


print "---- Loading COEX network"
COEX = nx.read_edgelist("./Data/Covid_Networks/%s_COEX.net"%(tissue))
COEX_mat = nm.Make_Adj(COEX, nodes, node2ind)

print "---- Loading GI network"
GI = nx.read_edgelist("./Data/Covid_Networks/%s_GI.net"%(tissue))
GI_mat = nm.Make_Adj(GI, nodes, node2ind)

nb_genes = len(nodes)

d = int( math.sqrt(nb_genes/2.) )
features = [str(i) for i in range(d)]

# Dataset in pairs, PPI + COEX
Solver = mf.iCell_SSNMTF(max_iter=200, verbose = 10)
G, S = Solver.Solve_MUR([PPI_mat, COEX_mat, GI_mat], d, init="SVD")
nm.Save_Matrix_Factor(G, "./Data/G matrix factors/%s_B.csv"%(tissue, nodes, features)
nm.Save_Matrix_Factor(S[0], "./Data/S matrix factors/.%s_S_PPI.csv"%tissue, features, features)
nm.Save_Matrix_Factor(S[1], "./Data/S matrix factors/.%s_S_COEX.csv"%tissue, features, features)
nm.Save_Matrix_Factor(S[2], "./Data/S matrix factors/.%s_S_GI.csv"%tissue, features, features)

nb_gene = len(nodes)
GT = np.transpose(G)
Raw_Net = np.dot(G,GT)
ind = int( float(nb_gene)*0.99 )

Int_Net = nx.Graph()
for i in range(nb_gene):
    sorted_row = np.sort(Raw_Net[i], kind='mergesort')
    val  = sorted_row[ind]
    genei = nodes[i]
    for j in range(nb_gene):
        genej = Glist[j]
        if (genei != genej) and Raw_Net[i][j]>= val:
            if not(Int_Net.has_edge(genei, genej)):
                Int_Net.add_edge(genei, genej, weigth=Raw_Net[i][j])
print("%s: Integrated net has %i nodes and %i edges"%(tissue, Int_Net.number_of_nodes(), Int_Net.number_of_edges()))
nx.write_edgelist(Int_Net, f'./Data/iCells/{tissue}_iCell.edgelist', data=['weight'])
