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
\nabla ^2 u &= 0 ;\;\ \textrm{  for  } 0<x<L,0<y<K \\
u(x,0)      &= 0 ;\;\ \textrm{  for  } 0\leq x\leq L \\
u(0,y)      &= u(L,y) = 0 ;\;\ \textrm{  for  } 0\leq y\leq K \\
u(x,K)      &= f(x) ;\;\ \textrm{  for  } 0\leq x\leq L 
\end{align}
```

The retangle is defined as domain $\Omega$. In such problem, we want a function that is harmonic on $\Omega$, equals $f(x)$ on the upper side and zero on the lower side  and two vertical side. Using the same technique of _separation of variable_, we can assume the structure of $u$ follows 


```{math}
:label: eq214
\begin{align}
u(x,y) = X(x)Y(y)
\end{align}
```

we have 

```{math}
:label: eq215
\begin{align}
X^{''}Y+XY^{''}=0
\end{align}
```

Dividing the entire {eq}`eq215` by $XY$, we can get  


```{math}
:label: eq215
\begin{align}
\frac{X^{''}}{X}+\frac{Y^{''}}{Y}=0
\end{align}
```

Similar to how we deal with heat and wave equations, the first part is only related to $x$ and the second part is only related to $y$. 
