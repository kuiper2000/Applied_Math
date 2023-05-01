(Laplace)=
# Week 13:  Laplace Equation  
## Problem Setup and Initial Conditions
While Laplace equation is nothing new to readers since we have introduced two related special functions in the chapter of {ref}`SpecialF`, Legendre polynomial and Bessel function. While those special functions form a complete set of solution of diffusion equation for the corresponding coordinate, there is an alternative way of solving diffusion problem when certain boundary conditions are satisfied. (i.e., Sturm-Liouville theorem is satisfied). 


Let's begin with a diffusion problem on a Retangle (figure below),


```{figure} Dirichlet_1.png
---
name: FIG12
scale: 20%
---
The Dirichlet problem on a retangle. 
```   

the figure above can be written as {eq}`eq213`

```{math}
:label: eq213
\begin{align}
\nabla ^2 u &= 0 \textrm{for } 0<x<L,0<y<K \\
u(x,0)      &= 0 \textrm{for } 0\leq x\leq L \\
u(0,y)      &= u(L,y) = 0 \textrm{for } 0\leq y\leq K \\
u(x,K)      &= f(x) \textrm{for } 0\leq x\leq L 
\end{align}
```
