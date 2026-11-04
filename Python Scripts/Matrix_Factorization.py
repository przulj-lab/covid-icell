# -*- coding: UTF-8 -*-

""" Non-negative matrix tri-factorization (numpy)"""
# Author: NMD

import numpy as np
import scipy.sparse.linalg as ssl

epsilon = 1e-5

def Get_Hard_Clusters(M, nodes):
	n,k = np.shape(M)
	Clusters = [[] for i in range(k)]
	for i in range(n):
		idx = np.argmax(M[i])
		Clusters[idx].append(nodes[i])
	return Clusters


class iCell_SSNMTF:
	#
	#  Mol[i] = G.S[i].G^T
	#  subject to G >=0

	def __init__(self, max_iter=1000, verbose = 10):
		self.max_iter = max_iter
		self.verbose = verbose


	def Score(self, MOLs, G, S, norm_Mol):
		Square_Error = 0.
		norm_Err_MOLs = []
		rel_MOLs = []
		GT = np.transpose(G)
		for i in range(len(MOLs)):
			Err_MOLs = MOLs[i] - np.matmul(G,  np.matmul(S[i], GT))
			norm_Err_MOLs.append( np.linalg.norm(Err_MOLs, ord='fro'))
			rel_MOLs.append( norm_Err_MOLs[i] / norm_Mol[i] )
			Square_Error += norm_Err_MOLs[i]

		return Square_Error, rel_MOLs



	def Solve_MUR(self, MOLs, k, init="rand"):
		n,_ = np.shape(MOLs[0])
		nb_mols = len(MOLs)

		norm_MOLs = [np.linalg.norm(MOLs[i], ord='fro') for i in range(nb_mols)]

		if init == "rand":
			print " * Using random initialization"
			G = np.random.rand(n, k)
			#S = [np.random.rand(k_g,k_g) for i in range(nb_mols)]

		elif init == 'SVD':
			print " * Using SVD based initialization"

			G = np.zeros((n, k)) + epsilon
			#S = [np.zeros((k_g,k_g)) + epsilon for i in range(nb_mols)]

			nbf = float(nb_mols)
			SMOL = np.zeros((n,n))
			for i in range(nb_mols):
				SMOL += MOLs[i]
			print " * -- Eig decomposition on MOL to get Gi.Si.B^T"
			M, N = ssl.eigsh(SMOL, k)
			for i in range(n):
				for j in range(k):
					if N[i][j] > 0.:
						G[i][j] = N[i][j]



		else :
			print "Unknown initializer: %s"%(init)
			exit(0)

		# Computing init S
		S = []
		GT = np.transpose(G)
		GTG = np.matmul( GT, G)
		#GTG_inv = np.linalg.pinv(GTG, rcond=1e-5)
		GTG_inv = np.linalg.inv(GTG)
		G_GTG_inv = np.matmul(G, GTG_inv)
		GTG_inv_GT = np.matmul(GTG_inv, GT)
		for i in range(nb_mols):
			S.append( np.matmul(GTG_inv_GT, np.matmul(MOLs[i], G_GTG_inv)) )

		OBJ, RELs = self.Score(MOLs, G, S, norm_MOLs)
		print " - Init:\t OBJ:%.4f\t REL:"%(OBJ), RELs


		#Begining M.U.R.
		MolT = [np.transpose(Mol) for Mol in MOLs]

		for it in range(1,self.max_iter+1):

			### First, update B, which depends on multiple decompositions

			# Update rule for G, from all MOL[i] = G.Si.G^T
			Ge = np.zeros((n,k))
			Gd = np.zeros((n,k))
			for i in range(nb_mols):
				MOLi_G_Si = np.matmul(MOLs[i], np.matmul(G, S[i]))
				MOLi_G_Si_pos = (np.absolute(MOLi_G_Si)+MOLi_G_Si)/2.0
				MOLi_G_Si_neg = (np.absolute(MOLi_G_Si)-MOLi_G_Si)/2.0

				Si_GTG_Si = np.matmul(S[i], np.matmul(GTG, S[i]))
				Si_GTG_Si_pos = (np.absolute(Si_GTG_Si)+Si_GTG_Si)/2
				Si_GTG_Si_neg = (np.absolute(Si_GTG_Si)-Si_GTG_Si)/2

				G_Si_GTG_Si_pos = np.matmul(G,Si_GTG_Si_pos)
				G_Si_GTG_Si_neg = np.matmul(G,Si_GTG_Si_neg)

				Ge = Ge + MOLi_G_Si_pos + G_Si_GTG_Si_neg
				Gd = Gd + MOLi_G_Si_neg + G_Si_GTG_Si_pos
			G_mult = np.sqrt(np.divide(Ge, Gd + epsilon))
			G = np.multiply(G, G_mult)

			# Computing all S[i] using closed formula
			GT = np.transpose(G)
			GTG = np.matmul( GT, G)
			#GTG_inv = np.linalg.pinv(GTG, rcond=1e-5)
			GTG_inv = np.linalg.inv(GTG)
			G_GTG_inv = np.matmul(G, GTG_inv)
			GTG_inv_GT = np.matmul(GTG_inv, GT)
			for i in range(nb_mols):
				S[i] = np.matmul(GTG_inv_GT, np.matmul(MOLs[i], G_GTG_inv))

			if (it%self.verbose == 0) or (it==1):
				OBJ, RELs = self.Score(MOLs, G, S, norm_MOLs)
				print " - It %i:\t OBJ:%.4f\t REL:"%(it, OBJ), RELs
		return G, S

class PSBLike_Drug_NMTF:

	#
	#  Min  || DTI - G.S.D^T ||^2_F + gamma D^T.Lchem.D
	#  subject to D >=0

	def __init__(self, max_iter=1000, verbose = 10):

		self.max_iter = max_iter
		self.verbose = verbose


	def Score(self, DTI, G, S, D, LChem, gamma, norm_DTI):

		DT = np.transpose(D)
		Err = DTI - np.matmul(G,  np.matmul(S, DT))
		norm_Err = np.linalg.norm(Err, ord='fro')**2
		rel = norm_Err / norm_DTI
		Err_Cst = gamma*np.trace(np.matmul(DT, np.matmul(LChem, D)))

		return norm_Err+Err_Cst, rel

	def Solve_MUR(self, DTI, Lchem, G, k1, k2, gamma, init="rand"):

		n,m = np.shape(DTI)
		norm_DTI = np.linalg.norm(DTI, ord='fro')

		if init == "rand":
			print " * Using random initialization"
			D = np.random.rand(m, k2)
		elif init == 'SVD':
			print " * Using SVD based initialization"
			D = np.zeros((m, k2)) + epsilon
			J,K,L = np.linalg.svd(DTI)
			for i in range(m):
				for j in range(k2):
					if L[j][i] > 0.:
						D[i][j] = L[j][i]
		else :
			print "Unknown initializer: %s"%(init)
			exit(0)



		# Computing init S

		GT = np.transpose(G)
		GTG = np.matmul( GT, G)
		GTG_inv = np.linalg.inv(GTG)
		GTG_inv_GT = np.matmul(GTG_inv, GT)

		DT = np.transpose(D)
		DTD = np.matmul( DT, D)
		DTD_inv = np.linalg.inv(DTD)

		S = np.matmul(GTG_inv_GT, np.matmul(DTI, np.matmul(D, DTD_inv)))
		OBJ, REL = self.Score(DTI, G, S, D, Lchem, gamma, norm_DTI)

		print " - Init:\t OBJ:%.4f\t REL:%.4f"%(OBJ, REL)

		#Begining M.U.R.

		DTIT = np.transpose(DTI)

		for it in range(1,self.max_iter+1):

			# Update rule for D

			DTIT_G_S = np.matmul(DTIT, np.matmul(G, S))
			DTIT_G_S_pos = (np.absolute(DTIT_G_S) + DTIT_G_S)/2.0
			DTIT_G_S_neg = (np.absolute(DTIT_G_S) - DTIT_G_S)/2.0

			ST_GTG_S = np.matmul(np.transpose(S), np.matmul(GTG, S))
			D_ST_GTG_S_pos = np.matmul(D, ((np.absolute(ST_GTG_S) + ST_GTG_S)/2.0))
			D_ST_GTG_S_neg = np.matmul(D, ((np.absolute(ST_GTG_S) - ST_GTG_S)/2.0))

			DNum = DTIT_G_S_pos + D_ST_GTG_S_neg
			DDenom = DTIT_G_S_neg + D_ST_GTG_S_pos

			# Contrain on G2, contribution from G2T.L2.G2

			DConst_pos = np.matmul(((np.absolute(Lchem) + Lchem)/2.0), D)
			DConst_neg = np.matmul(((np.absolute(Lchem) - Lchem)/2.0), D)

			#Check that we are note dividing by zero

			Dupdate = np.sqrt(np.divide(DNum + gamma*DConst_neg, DDenom + gamma*DConst_pos + epsilon))
			D = np.multiply(D, Dupdate) + epsilon

			# Computing S using closed formula

			DT = np.transpose(D)
			DTD = np.matmul( DT, D)
			DTD_inv = np.linalg.inv(DTD)

			S = np.matmul(GTG_inv_GT, np.matmul(DTI, np.matmul(D, DTD_inv)))

			if (it%self.verbose == 0) or (it==1):
				OBJ, REL = self.Score(DTI, G, S, D, Lchem, gamma, norm_DTI)
				print " - It %i:\t OBJ:%.4f\t REL:%.4f"%(it, OBJ, REL)
		return D, S
