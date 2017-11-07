import argparse
import sys
import sympy as sp
import typing
import math

from symdoc import Markdown, doit, symfunc, gsym

cmdline_parser = argparse.ArgumentParser()
Markdown.add_parser_option(cmdline_parser)

M = Markdown('kk3', title='kk3')
markdown, cmdline = M.markdown, M.cmdline


######################################################################
# 大域的な環境設定

import numpy as np
import sympy as sp

import matplotlib.pyplot as plt
from time import time as current_time
######################################################################

dim = 2
n = 3

markdown(r'''
グラフを描画する際の手法としてKK法というものがあります。
グラフの各頂点間に仮想のばねを導入し、その系のエネルギーを最小にさせることにで、グラフのいい描画を得られます。

KK法を実装するためのいくつかの関数群を説明します。
''')

@symfunc
def convertToMatrix(nodeList,adj_list):
    '隣接リスト -> 隣接行列'
    n = len(nodeList)
    adj_matrix = np.zeros([n,n])

    for i in range(n):
        edges_of_i = len(adj_list[i])
        for j in range(edges_of_i):
            for k in range(n):
                if(adj_list[i][j] == nodeList[k]):
                    adj_matrix[i][k] = 1

    return adj_matrix

@doit
def __convertToMatrix__():
	markdown(r'''
#convertToMatrix(nodeList,adj_list)
隣接リストを隣接行列に変換します。
		''')

@symfunc
def warshall_floyd(adj_matrix):
    n = adj_matrix.shape[0]
    # generate distance_matrix
    distance_matrix = -n*adj_matrix + n+1 -(n+1)*np.identity(n)

    for k in range(n):
        for i in range(n):
            for j in range(n):
                if(distance_matrix[i][j] > distance_matrix[i][k] + distance_matrix[k][j] ):
                    distance_matrix[i][j] = distance_matrix[i][k] + distance_matrix[k][j]
    return distance_matrix

@doit
def __warshall_floyd__():
	i,j,k,d=sp.symbols('i j k d')
	dij,dik,dkj = ['d_{i,j}','d_{i,k}','d_{k,j}']
	markdown(r'''
#{fwarshall_floyd}(adj_matrix)
次元$n$の隣接行列:$A$から、warshall floyd法を用いて、距離行列:$D$を計算します。実際の処理は以下の通りです。

1. 距離行列:Dの初期値を以下のように設定する$$D=-n*A+(n+1)*J-(n+1)*I$$

2. $k$を$[0,n-1]$で動かしながら、以下を繰り返す

	1. $i$を$[0,n-1]$で動かしながら、以下を繰り返す

		1. $j$を$[0,n-1]$で動かしながら、以下を繰り返す

			1. ${dij}>{dik}+{dkj}$ならば、${dij}={dik}+{dkj}$とする

3. $D$を出力する

''',
	**locals())
	return 0


@symfunc
def draw_2Dgraph(adj_matrix,pos,name='sample'):
    figname = 'pics/' + name + '.png'
    n,dim = pos.shape
    for i in range(n-1):
        for j in range(i,n):
            if(adj_matrix[i,j] == 1):
                plt.plot([pos[i,0],pos[j,0]],[pos[i,1],pos[j,1]],'k-')

    plt.plot(pos.T[0],pos.T[1],'go')
    plt.title(name)
    plt.axes().set_aspect('equal', 'datalim')
    plt.savefig(figname)
    plt.clf()

@doit
def __draw_2Dgraph__():
	markdown(r'''
#draw_2Dgraph(adj_matrix,pos,name='sample')
グラフの各頂点の２次元上の配置:posと辺:adj_matrix、グラフ名:nameを受け取ってそのグラフを描画します。

		''')

@symfunc
def kk2(Pos,Lens,Cons,dim=dim,n=n,eps = 0.000001):
	start_time = current_time()
	ni=n-1
	P = sp.IndexedBase('P',shape=(ni,dim)) #ni個のPositions
	L = sp.IndexedBase('L') #L[j] is Natural Lens between p and P[j]
	K = sp.IndexedBase('K') #K[j] is p,P[j]間のばね定数
	p = sp.IndexedBase('p') #Position of one of Pで動かす点
	q = sp.IndexedBase('q')#最適化の式の解
	#i番目をスライスするよう
	Pos2=np.concatenate([Pos,Pos],axis=0).astype(np.float64) #Pos2[i+1:ni+i+1]
	Len2=np.concatenate([Lens,Lens],axis=0) #Len2[i+1:ni+i+1,i]
	Cons2=np.concatenate([Cons,Cons],axis=0) #Cons2[i+1:ni+i+1,i]


	i,j, d, x =[sp.Idx(*spec) for spec in [('i',ni),('j', ni),('d',dim),('x',dim)]]
	i_range, j_range = [(idx, idx.lower, idx.upper) for idx in [i, j]]
	d_range = (d,d.lower,d.upper)

	distance_j= sp.sqrt(sp.Sum((p[d] - P[j,d])*(p[d] -P[j,d]), d_range))
	Energy_j = (K[j]*(distance_j-L[j])*(distance_j-L[j])/2).doit()
	Energy = (sp.Sum(Energy_j,j_range).doit())

	variables = [p[d] for d in range(dim)]
	jac_sp,hess_sp = [sp.Matrix([Energy]).jacobian(variables),sp.hessian(Energy,variables)]
	jac,hess = [(sp.lambdify((K,L,P,p),f,dummify = False)) for f in [jac_sp,hess_sp]]

	#Sum(jac_sp[d]**2,d_range)の代わり
	abs_grad_E_sp=0
	for i in range(dim):
		abs_grad_E_sp += jac_sp[i]*jac_sp[i]
	abs_grad_E_sp=sp.sqrt(abs_grad_E_sp)

	abs_grad_E=sp.lambdify((K,L,P,p),abs_grad_E_sp,dummify=False)

	#ラムダ式の変数を略記する
	def vars_i(i):
		return Cons2[i+1:ni+i+1,i],Len2[i+1:ni+i+1,i],Pos2[i+1:ni+i+1],Pos2[i]

	max_delta = sp.oo
	max_delta_number = 0
	loop=0


	while max_delta>eps:
		#方程式を解いて一点を動かす
		ds=np.linalg.solve(hess(*vars_i(max_delta_number)),jac(*vars_i(max_delta_number)).flatten())
		Pos2[max_delta_number]-=ds
		Pos2[max_delta_number+n]-=ds

		#abs_grad_Eの最大となるところを探す
		max_delta=0
		for i_move in range(n):
			if(abs_grad_E(*vars_i(i_move))>max_delta):
				max_delta=abs_grad_E(*vars_i(i_move))
				max_delta_number=i_move
		loop+=1
	print('loop=%d'%loop)
	print("かかった時間: {}".format(current_time() - start_time))
	Pos_after=sp.Matrix(Pos2[:n])
	markdown(r'''
##kk2(Pos=${Pos}$,Lens=${Lens}$,Cons=${Cons}$)
ループを${loop}$回まわし、最適化した結果、以下が得られました。
$$Pos={Pos_after}$$
''',
			**locals())


	return Pos2[:n]




def kk(name,nodeList,adj_list,L0=1,K=1,dim=2,eps=0.000001):
	#kk2を呼ぶ

	n=len(nodeList)
	adj_matrix=convertToMatrix(nodeList,adj_list)
	distance_matrix=warshall_floyd(adj_matrix)
	Pos=np.zeros((n,dim))
	if dim==2:
		Pos=np.array([[L0*np.cos(a*2*np.pi/n)/2,L0*np.sin(a*2*np.pi/n)/2] for a in range(n)])
	else:
		Pos = np.array(L0*np.random.rand(n*dim).reshape(n,dim))
	L=L0/np.max(distance_matrix)
	Lens=L*distance_matrix
	Cons=K/((distance_matrix+np.identity(n))*(distance_matrix+np.identity(n)))-K*np.identity(n)
	print('{}を最適化'.format(name))
	Pos_sp,Lens_sp,Cons_sp=[sp.Matrix(M) for M in [Pos,Lens,Cons]]
	markdown(r'''
#{name}の最適化
まず、kkを呼びます。

##kk(name=${name}$,nodeList=${nodeList}$,adj_list=${adj_list}$)
距離行列:distance_matrixは以下のよう求まります。
$$distance\_matrix={distance_matrix}$$
kk2に以下の変数を渡します。
$$Lens={Lens_sp}$$
$$Cons={Cons_sp}$$
$$Pos={Pos_sp}$$
kk2を呼び、最適化後、{name}のグラフを描画します
		''',**locals())
	aftermove=kk2(Pos,Lens,Cons,dim=dim,n=n,eps=eps)
	draw_2Dgraph(adj_matrix,Pos,name=name+'_before')
	draw_2Dgraph(adj_matrix,aftermove,name=name+'_after')
	markdown(r'''
グラフ:{name}は以下のように描画されました。

![最適化前](pics/{name}_before.png)
![最適化後](pics/{name}_after.png)

''',**locals())


@doit
def __kk__():
	i,j,k,L,d,n,L0=sp.symbols('i j k L d n L0')
	lij, dij,kij= ['l_{i,j}','d_{i,j}','k_{i,j}']
	initPos_i=sp.Matrix([L0*sp.cos(i*2*sp.pi/n)/2,L0*sp.sin(i*2*sp.pi/n)/2])

	markdown(r'''
#kk(name,nodeList,adj_list,L0=1,K=1,dim=2,eps=0.000001)
kk2を呼び出すための関数です。

グラフの名前:name、その頂点のリスト:nodeList、グラフの隣接リスト:adj_listを受け取り、グラフに仮想のばねを導入し、そのばね定数と自然長をグラフ描画の適当な初期配置と合わせて、kk2に渡します。

1.{fconvertToMatrix}、{fwarshall_floyd}を使い、nodeList,adj_listから、距離行列:distance_matrixを求める。

2.以下のように$i,j$間の、自然長:{lij}、ばね定数:{kij}を定義し、$Lens=({lij})$,$Cons=({kij})$とする。
$${lij}=L0*{dij}/\left(\max_\left(i,j\in[0,n-1]\right){dij}\right)$$
$${kij}=K/({dij}^2)$$

3.各頂点の適当な初期配置:Posを与える。ここでは、描画する次元を2ならば、平面上に直径L0の正n角形を与える。すなわち、
$$Pos_i={initPos_i}$$
  次元が2でなければ、ランダムな値を与える。

4.kk2(P,L,K,dim=dim,n=n,eps=eps)をよび、最適化させたものを{fdraw_2Dgraph}で描画し、pics/(name).pngに保存する。


		''',**locals())


@doit
def __kk2__():
	ni=n-1
	P = sp.IndexedBase('P',shape=(ni,dim)) #ni個のPositions
	L = sp.IndexedBase('L') #L[j] is Natural Lens between p and P[j]
	K = sp.IndexedBase('K') #K[j] is p,P[j]間のばね定数
	p = sp.IndexedBase('p') #Position of one of Pで動かす点
	q = sp.IndexedBase('q')#最適化の式の解


	i,j, d, x =[sp.Idx(*spec) for spec in [('i',ni),('j', ni),('d',dim),('x',dim)]]
	i_range, j_range = [(idx, idx.lower, idx.upper) for idx in [i, j]]
	d_range = (d,d.lower,d.upper)

	distance_j= sp.sqrt(sp.Sum((p[d] - P[j,d])*(p[d] -P[j,d]), d_range))
	Energy_j = (K[j]*(distance_j-L[j])*(distance_j-L[j])/2)
	Energy = (sp.Sum(Energy_j,j_range))

	variables = [p[d] for d in range(dim)]
	jac_sp,hess_sp = [sp.Matrix([Energy]).jacobian(variables),sp.hessian(Energy,variables)]
	jac,hess = [(sp.lambdify((K,L,P,p),f,dummify = False)) for f in [jac_sp,hess_sp]]

	#Sum(jac_sp[d]**2,d_range)の代わり
	abs_grad_E_sp=0
	for i in range(dim):
		abs_grad_E_sp += jac_sp[i]*jac_sp[i]
	abs_grad_E_sp=sp.sqrt(abs_grad_E_sp)

	abs_grad_E=sp.lambdify((K,L,P,p),abs_grad_E_sp,dummify=False)

	markdown(r'''
#kk2(Pos,Lens,Cons,dim=2,eps=0.000001)
グラフ描画の初期配置:Pos、ばね定数:Cons、ばねの自然長:Lensをもらい、エネルギーが最小(初期配置によっては極小)の近似的な配置を返します。

1.動かす点:$p$を一つとり、Pos,Cons,Lensから点$p$を除いたものを$P,K,L$と置く。

2.動かす点$p$に対し、そのエネルギー:$E$,その勾配の絶対値:$\left|gradE\right|$、ヤコビアン:$J$、ヘシアン$H$を以下のよう定義する。
$$E={Energy}$$
$$\left|gradE\right|={abs_grad_E_sp}$$
$$J={jac_sp}$$
$$H={hess_sp}$$

3.動かす点:$p$を変えながら、$\left|gradE\right|$が最大となる点探す。

4.$p$を以下の手順で最適化する。

	1.$x$を変数として以下の連立方程式を解く。
$$H x=J$$

	2.その解を用いて、点$p$を更新する。
$$p=p-x$$

5.動かす点:$p$を変えながら、$\left|gradE\right|$が最大となる点探し、$\left|gradE\right|>eps$ならば4に戻る。


		''',**locals())


markdown(r'''
実際にいくつかのグラフで試してみましょう。
	''')

kk('triangle',[0,1,2],[[1,2],[0,2],[0,1]])
kk('double_triangle',[0,1,2,3,4,5],[[2,3],[4,5],[0,3,5],[0,2],[1,5],[1,2,4]])
kk('4d_cube',[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],[[1,2,4,8],[0,3,5,9],[0,3,6,10],[1,2,7,11],[0,5,6,12],
[1,4,7,13],[2,4,7,14],[3,5,6,15],[0,9,10,12],[1,8,11,13],[2,8,11,14],[3,9,10,15],[4,8,13,14],[5,9,12,15]
,[6,10,12,15],[7,11,13,14]])