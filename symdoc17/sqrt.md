---
title: How to calcurate sqrt using NewtonRapshon method
layout: page
---






##Newton-Rapshon法による$F(x):=- a + x^{2}=0$の求解
Newton-Rapshon法を用いて、
方程式$F()=- a + x^{2}=0$の近似解を求めます。\n
以下の漸化式による数列は$F(x)=0$の解$α$に収束します。
$$x_{n+1}=x_{n} - \frac{F{\left (x_{n} \right )}}{\frac{d}{d x_{n}} F{\left (x_{n} \right )}}=\frac{a + x_{n}^{2}}{2 x_{n}}$$



$F(x)=x^{2} - 2$とすると、漸化式は
$$x_{n+1}=\frac{x_{n}^{2} + 2}{2 x_{n}}$$



$x_0=1$として、$\sqrt{2}≈x_6=1.414213562373095$が近似解として得られます。



##Newton-Rapshon法による$F(x):=a^{x}=0$の求解
Newton-Rapshon法を用いて、
方程式$F()=a^{x}=0$の近似解を求めます。\n
以下の漸化式による数列は$F(x)=0$の解$α$に収束します。
$$x_{n+1}=x_{n} - \frac{F{\left (x_{n} \right )}}{\frac{d}{d x_{n}} F{\left (x_{n} \right )}}=x_{n} - \frac{1}{\log{\left (a \right )}}$$



$F(x)=2^{x}$とすると、漸化式は
$$x_{n+1}=x_{n} - 1.44269504088896$$



$x_0=1$としたところ、
NewtonRapshon法で収束しませんでした。


