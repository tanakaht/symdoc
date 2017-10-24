import argparse
import sys
import sympy as sp
import typing

from symdoc import Markdown, doit, symfunc, gsym

cmdline_parser = argparse.ArgumentParser()
Markdown.add_parser_option(cmdline_parser)

M = Markdown('t3d', title='T3D (Transform 3D) for Python')
markdown, cmdline = M.markdown, M.cmdline


######################################################################
# 大域的な環境設定

import numpy as np
import sympy as sp

######################################################################
# ユーティリティ関数
dim = 2
ni = 2




def kk2(Pos,Lengthes,Cons,dim=dim,ni=ni,eps = 0.000001):
	P = sp.IndexedBase('P',shape=(dim,ni)) #ni個のPositions
	L = sp.IndexedBase('L') #L[i,j] is Natural Lengthes between pi and pj
	K = sp.IndexedBase('K') #K[i,j] is ij間のばね定数
	p = sp.IndexedBase('p') #Position of one of P

	i, j, d, x = [sp.Idx(*spec) for spec in [('i',ni),('j', ni), ('d', dim), ('x',dim)]]
	i_range, j_range, d_range = [(idx, idx.lower, idx.upper) for idx in [i, j, d]]

	distance_ij = sp.sqrt(sp.Sum((P[d,i] - P[d,j]) ** 2, d_range)).doit()
	Energy_ij = (K[i,j]*(distance_ij-L[i,j])*(distance_ij-L[i,j])/2)
	Potential_i = sp.Sum(Energy_ij,j_range).doit()
	total_energy = (sp.Sum(Potential_i,i_range).doit())/2 #2回足してるので2で割る

	delta_fomula_i = sp.sqrt(sp.Sum(sp.diff(Potential_i,P[d,i])**2,d_range)).doit()
	delta_formula_i_lambda = sp.lambdify(i,
								sp.lambdify((P,L,K),delta_fomula_i)(Pos,Lengthes,Cons))

	#変位
	__diffx__ =sp.diff(Potential_i,P[x,i])
	delta_x_lambda = sp.lambdify((x,p),
				 sp.lambdify((P,L,K),sp.Sum(sp.diff(__diffx__,P[d,i])*p[d],d_range)+__diffx__)(Pos,Lengthes,Cons))
	max_delta = 0
	max_delta_number = 0

	#serchmax
	for i in range(ni):
		if(delta_fomula_i_lambda(i)>max_delta):
			max_delta = delta_fomula_i_lambda(i)
			max_delta_number = i
	#P[d,mux_delta_number]を最適化
	while max_delta>eps:
		print('{0}'.format(max_delta_number))
		while delta_fomula_i_lambda(max_delta_number) > eps :
			s = sp.solve([delta_x_lambda(x,p) for x in range(dim)],[p[x] for x in range(x)]).doit()
			for d in range(dim):
				Pos[d,i] -= s[p[d]]

		#serch max
		max_delta = 0
		max_delta_number = 0
		for i in range(ni):
			if(delta_fomula_i_lambda(i,Pos,Lengthes,Cons)>max_delta):
				max_delta = delta_fomula_i_lambda(i,Pos,Lengthes,Cons)
				max_delta_number = i

	return Pos

#KK2を呼ぶんでみる
ni = 3
dim = 2
Pos= np.array([[1,0],[0,1],[-1,1]]).T
Lengthes = np.array([[0,1,1],
					 [1,0,1],
					 [1,1,0]])
Cons = np.array([[0,1,1],
			     [1,0,1],
				 [1,1,0]])
kk2(Pos,Lengthes,Cons,dim=dim,ni=ni)


