import argparse
import sys
import sympy as sp
import typing
import math

from symdoc import Markdown, doit, symfunc, gsym

cmdline_parser = argparse.ArgumentParser()
Markdown.add_parser_option(cmdline_parser)

M = Markdown('kk2', title='kk2')
markdown, cmdline = M.markdown, M.cmdline


######################################################################
# 大域的な環境設定

import numpy as np
import sympy as sp

######################################################################

dim = 2
n = 3
#dim*n行列からi番目のを取り除く
def del_i(Matrix,i):
	Matrix2 = np.zeros((dim,n-1))
	for d in range(dim):
		Matrix2[d]=np.delete(Matrix[d],i)
	return Matrix2
#n*nから同様
def del_Matrix_i(Matrix,i):
	return np.delete(Matrix[i],i)

def KroneckerDelta(i,j):
	flag=0
	if i==j:
		flag=1
	return flag

mymodules =[{'KroneckerDelta':KroneckerDelta},'numpy']

def kk2(Pos,Lengthes,Cons,dim=dim,n=n,eps = 0.000001):
	P = sp.IndexedBase('P',shape=(dim,n)) #n個のPositions
	L = sp.IndexedBase('L') #L[i,j] is Natural Lengthes between pi and pj
	K = sp.IndexedBase('K') #K[i,j] is ij間のばね定数
	p = sp.IndexedBase('p') #Position of one of Pで動かす点
	q = sp.IndexedBase('q')#最適化の式の解

	i,j, d, x =[sp.Idx(*spec) for spec in [('i',n),('j', n),('d',dim),('x',dim)]]
	i_range, j_range = [(idx, idx.lower, idx.upper -1) for idx in [i, j]]
	d_range = (d,d.lower,d.upper)

	def exc_i(a):
		return 1*(a>=i)+a

	distance_j= sp.sqrt(sp.Sum((p[d] - P[d,exc_i(j)])*(p[d] -P[d,exc_i(j)]), d_range))
	Energy_j = (K[i,exc_i(j)]*(distance_j-L[i,exc_i(j)])*(distance_j-L[i,exc_i(j)])/2).doit()
	Potential = (sp.Sum(Energy_j,j_range).doit())
	Energy_j_diffx=sp.diff(Energy_j,p[1])

	diff_x=sp.simplify(sp.diff(Potential,p[d]))
	delta_fomula = sp.sqrt(sp.Sum(sp.diff(Potential,p[d])**2,d_range))
	delta_fomula_i_lambda = sp.lambdify((p,i,P,L,K),delta_fomula,dummify=False)


	delta_x=sp.Sum(sp.diff(sp.diff(Potential,p[x]),p[d])*q[d],d_range)+sp.diff(Potential,p[x])
	delta_x_lambda = sp.lambdify((x,q,p,i,P,L,K),delta_x)

	markdown(r'''
$$distance_j={distance_j}$$
$$Energy_j={Energy_j}$$
$$Enediff={Energy_j_diffx}$$
$$Potential={Potential}$$
$$diff_x={diff_x}$$
$$delta_fomula={delta_fomula}$$
''',
			**locals())

	max_delta = 0
	max_delta_number = 0
	#確認用
	def xmi(i,m):
		return Cons[i,m]*((Pos[0,m]-Pos[0,i])-Lengthes[i,m]*(Pos[0,m]-Pos[0,i])/
			math.sqrt((Pos[0,m]-Pos[0,i])**2+(Pos[1,m]-Pos[1,i])**2))
	def ymi(i,m):
		return Cons[i,m]*((Pos[1,m]-Pos[1,i])-Lengthes[i,m]*(Pos[1,m]-Pos[1,i])/
			math.sqrt((Pos[0,m]-Pos[0,i])**2+(Pos[1,m]-Pos[1,i])**2))
	def deltam(m):
		tmpx = 0
		tmpy = 0
		for i in range(n):
			if i!=m:
				tmpx+=xmi(i,m)
				tmpy+=ymi(i,m)
		return math.sqrt(tmpx**2+tmpy**2)

	for i_move in range(n):
		p_move=Pos.T[i_move]
		print(p_move)
#		print(delta_fomula_i_lambda(p_move,i_move,Pos,Lengthes,Cons))
		print(deltam(i_move))
#		print(delta_x_lambda(0,q2,p_move,i_move,Pos,Lengthes,Cons))
#		ds=sp.solve([delta2_x_lambda(x,q,p_move,i_move,Pos,Lengthes,Cons) for x in range(dim)],[q[x] for x in range(dim)])
		for d in range(dim):
#			Pos[d,i_move] -=ds[q[d]]
			d==d


	return Pos

#KK2を呼ぶんでみる
n = 3
dim = 2
Pos= np.array([[1,2],[3,1],[5,3]]).T
Lengthes = np.array([[0,1,1],
					 [1,0,1],
					 [1,1,0]])
Cons = np.array([[0,1,1],
			     [1,0,1],
				 [1,1,0]])
kk2(Pos,Lengthes,Cons,dim=dim,n=n)


