import argparse
import sys
import sympy as sp
import typing

from symdoc import Markdown, doit, symfunc, gsym

cmdline_parser = argparse.ArgumentParser()
Markdown.add_parser_option(cmdline_parser)

M = Markdown('kk', title='kk')
markdown, cmdline = M.markdown, M.cmdline


######################################################################
# 大域的な環境設定

import numpy as np
import sympy as sp

######################################################################
# ユーティリティ関数
dim = 2
n = 2

def KroneckerDelta(i,j):
	flag=0
	if i==j:
		flag=1
	return flag

mymodules =[{'KroneckerDelta':KroneckerDelta}, 'numpy']


def kk2(Pos,Lengthes,Cons,dim=dim,n=n,eps = 0.000001):
	P = sp.IndexedBase('P',shape=(dim,n)) #n個のPositions
	L = sp.IndexedBase('L') #L[i,j] is Natural Lengthes between pi and pj
	K = sp.IndexedBase('K') #K[i,j] is ij間のばね定数
	p = sp.IndexedBase('p') #Position of one of P

	i, j, d, x = [sp.Idx(*spec) for spec in [('i',n),('j', n), ('d', dim), ('x',dim)]]
	i_range, j_range, d_range = [(idx, idx.lower, idx.upper) for idx in [i, j, d]]



	distance_ij = sp.sqrt(sp.Sum((P[d,i] - P[d,j]) ** 2, d_range)).doit()
	Energy_ij = (K[i,j]*(distance_ij-L[i,j])*(distance_ij-L[i,j])/2)
	Potential_i = sp.Sum(Energy_ij,j_range)
	print(Potential_i)
	Potential_i_lambda = sp.lambdify((i,P,L,K),Potential_i)
#	print(Potential_i_lambda(1,Pos,Lengthes,Cons))
#	total_energy = (sp.Sum(Potential_i,i_range).doit())/2 #使わない

	delta_fomula_i = sp.sqrt(sp.Sum(sp.diff(Potential_i,P[d,i])**2,d_range))
	print('delta_fomula_i{}',delta_fomula_i)
	delta_fomula_i_lambda = lambda i,P,K,L:0
#	delta_fomula_i_lambda = sp.lambdify((i,P,L,K),delta_fomula_i,modules=mymodules)

	#変位
	__diffx__ =sp.diff(Potential_i,P[x,i])
	print({},__diffx__)
	delta_x_lambda = lambda x,p:1
#	delta_x_lambda = (lambda x,p:sp.lambdify((x,p,P,L,K),sp.Sum(sp.diff(__diffx__,P[d,i])*p[d],d_range)+__diffx__,modules=mymodules)(x,p,Pos,Lengthes,Cons))
	max_delta = 0
	max_delta_number = 0

	markdown(r'''
$$distanceij={distance_ij}$$
$$Energy_ij={Energy_ij}$$
$$Potential_ij={Potential_i}$$
$$delta_fomula_i={delta_fomula_i}$$
''',
			**locals())

	#serchmax
	for s in range(n):
		if(delta_fomula_i_lambda(s,Pos,Lengthes,Cons)>max_delta):
			max_delta = delta_fomula_i_lambda(s,Pos,Lengthes,Cons)
			max_delta_number = s
	#P[d,mux_delta_number]を最適化
	while max_delta>eps:
		print('{0}'.format(max_delta_number))
		while delta_fomula_i_lambda(max_delta_number,Pos,Lengthes,Cons) > eps :
			s = sp.solve([delta_x_lambda(x,p) for x in range(dim)],[p[x] for x in range(x)]).doit()
			for d in range(dim):
				Pos[d,i] -= s[p[d]]

		#serch max
		max_delta = 0
		max_delta_number = 0
		for i in range(n):
			if(delta_fomula_i_lambda(i,Pos,Lengthes,Cons)>max_delta):
				max_delta = delta_fomula_i_lambda(i,Pos,Lengthes,Cons)
				max_delta_number = i

	return Pos

#KK2を呼ぶんでみる
n = 3
dim = 2
Pos= np.array([[1,0],[0,1],[-1,1]]).T
Lengthes = np.array([[0,1,1],
					 [1,0,1],
					 [1,1,0]])
Cons = np.array([[0,1,1],
			     [1,0,1],
				 [1,1,0]])
kk2(Pos,Lengthes,Cons,dim=dim,n=n)


