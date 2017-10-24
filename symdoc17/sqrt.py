import argparse
import sys
import sympy as sp
import typing
import numpy as np

from symdoc import Markdown, doit, symfunc, gsym


cmdline_parser = argparse.ArgumentParser()
Markdown.add_parser_option(cmdline_parser)

M = Markdown('sqrt', title='How to calcurate sqrt using NewtonRapshon method')
markdown, cmdline = M.markdown, M.cmdline

x, a =sp.var('x a')



def NewtonRaphson(f,x,*args):
	F = sp.Function('F')
	newrap = x - F(x)/sp.diff(F(x),x)
	thisnewrap = sp.simplify(newrap.subs(F(x),f).doit())
	calnewrap = sp.lambdify((x,*args),thisnewrap)
	fanc = sp.simplify(f)
	def newtonraphson(*args2,xstart=1,eps = 10**(-20),solution='α'):
		xloop = xstart
		xback = xstart-1
		i = 0
		notfind = False
		while abs(xloop - xback) > eps:
			xback = xloop
			xloop = calnewrap(xloop,*args2)
			i = i+1
			if i>100:
				print('収束しませんでした')
				notfind = True
				break
		#startfordoc
		if M.cmdline.symdoc:
			fancargs = fanc.subs(*args,*args2)
			xn, n1 = sp.var('x_n {n+1}')
			thisnewrapxnargs = thisnewrap.subs(*args,*args2).subs(x,xn)
			markdown(r'''
$F(x)={fancargs}$とすると、漸化式は
$$x_{n1}={thisnewrapxnargs}$$
''',
			**locals())
			if notfind ==False:
				markdown(r'''
$x_0={xstart}$として、${solution}≈x_{i}={xloop}$が近似解として得られます。
''',
				**locals())
			if notfind:
				markdown(r'''
$x_0={xstart}$としたところ、
収束しませんでした。
''',
				**locals())
		#endfordoc
		return xloop

	#startfordoc
	if M.cmdline.symdoc:
		xn, xn1, n1 = sp.var('x_n x_{n+1} {n+1}')
		thisnewrapxn = thisnewrap.subs(x,xn)
		newrapxn = newrap.subs(x,xn)
		markdown(r'''
##Newton-Rapshon法による${F}(x):={f}=0$の求解
Newton-Rapshon法を用いて、
方程式${F}(x)={f}=0$の近似解を求めます。\n
以下の漸化式による数列は${F}(x)=0$の解$α$に収束します。
$$x_{n1}={newrapxn}={thisnewrapxn}$$
''',
		**locals())
	#endfordoc
	return newtonraphson

sqrt = NewtonRaphson(x**2-a,x,a)
sqrt(2,solution = sp.sqrt(2))
pow = NewtonRaphson(a**(x),x,a)
pow(2)

ff = NewtonRaphson(x/sp.tan(x)-1/a,x,a)
ff(1,xstart = 5)