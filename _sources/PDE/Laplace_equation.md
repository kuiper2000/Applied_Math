(Laplace)=
# Week 13:  Laplace Equation  
## Problem Setup and Initial Conditions
While Laplace equation is nothing new to readers since we have introduced two related special functions in the chapter of {ref}`SpecialF`, Legendre polynomial and Bessel function. While those special functions form a complete set of solution of diffusion equation for the corresponding coordinate, there is an alternative way of solving diffusion problem when certain boundary conditions are satisfied. (i.e., Sturm-Liouville theorem is satisfied). 



## Dirichlet Problem for a Rectangle
Let's begin with a diffusion problem on a Rectangle (figure below),


```{figure} Dirichlet_1.png
---
name: FIG12
scale: 20%
---
The Dirichlet problem on a rectangle. 
```   

the figure above can be written as {eq}`eq213`

```{math}
:label: eq213
\begin{align}
\nabla ^2 u &= 0 \;\; \textrm{  for  } 0<x<L,0<y<K \\
u(x,0)      &= 0 \;\; \textrm{  for  } 0\leq x\leq L \\
u(0,y)      &= u(L,y) = 0 \;\; \textrm{  for  } 0\leq y\leq K \\
u(x,K)      &= f(x) \;\; \textrm{  for  } 0\leq x\leq L 
\end{align}
```

The rectangle is defined as domain $\Omega$. In such problem, we want a function that is harmonic on $\Omega$, equals $f(x)$ on the upper side and zero on the lower side  and two vertical side. Using the same technique of _separation of variable_, we can assume the structure of $u$ follows 


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
:label: eq216
\begin{align}
\frac{X^{''}}{X}+\frac{Y^{''}}{Y}=0
\end{align}
```

Similar to how we deal with heat and wave equations, the first part is only related to $x$ and the second part is only related to $y$. This indicates that both are the same constant with opposite sign. i.e., 

```{math}
:label: eq217
\begin{align}
\frac{X^{''}}{X}+\lambda & =0 \;\; \textrm{  for } X(0)=X(L)=0\\
\frac{Y^{''}}{Y}-\lambda & =0 \;\; \textrm{  for } Y(0)=0 \\
\end{align}
```


It's not hard to find that {eq}`eq217` has a solution of (also see solution )

```{math}
:label: eq218
\begin{align}
\lambda_n = \frac{n^2\pi^2}{L^2},X_n(x)=\sin(\frac{n\pi x}{L})
\end{align}
```

For solution in $y$ direction, 

```{math}
:label: eq219
\begin{align}
Y_n(y)=\sinh (\frac{n\pi y}{L})
\end{align}
```

The whole solution can therefore be written as 


```{math}
:label: eq220
\begin{align}
u(x,y) = \sum_{n=1}^{n=\infty} c_n \sin(\frac{n\pi x}{L})\sinh(\frac{n\pi y}{L})
\end{align}
```

For the boundary condition at $y=K$, we know 

```{math}
:label: eq221
\begin{align}
u(x,K) = \sum_{n=1}^{n=\infty} c_n \sin(\frac{n\pi x}{L})\sinh(\frac{n\pi K}{L})
\end{align}
```

which infers that 


```{math}
:label: eq222
\begin{align}
c_n \sinh(\frac{n\pi K}{L}) = \frac{2}{L}\int_{0}^{L}  f(\xi)\sin(\frac{n\pi \xi}{L})d\xi
\end{align}
```

Substitute {eq}`eq222` back to {eq}`eq220`, we have 

```{math}
:label: eq223
\begin{align}
u(x,y) = \sum_{n=1}^{\infty}  \frac{2}{L}\int_{0}^{L}  (f(\xi)\sin(\frac{n\pi \xi}{L})d\xi ) \sin(\frac{n\pi x}{L})\frac{\sinh(\frac{n\pi y }{L})}{\sinh (\frac{n\pi K }{L})}
\end{align}
```


While we only specify boundary condition at the upper end, we can approach the problem the same way if other edges have non-zero boundary conditions. (see figure below)


```{figure} Dirichlet_3.png
---
name: FIG13
scale: 20%
---
The Dirichlet problem on a rectangle with four specified boundaries. 
```   

The solution for the condition above can be decomposed into four cases and each case with only one boundary specified (The lower panel of {ref}`FIG13`).  The total solution is the linear combination of the four solution above. i.e., 

```{math}
:label: eq224
\begin{align}
u(x,y) = \sum_{i=1}^{4}  u_i(x,y)
\end{align}
```

## Dirichlet Problem for a Disk

If today, we setup a problem on a disk ($\Omega \; \in x^2+y^2<R^2$). It's a rare case when we use Cartesian coordinate but not so rare when we are using polar coordinate i.e., $x=r\cos(\theta)$ and $y=r\sin(\theta)$. One can consider Dirichlet's approach is an alternative way to solve the polar coordinate diffusion problem. 


With given boundary condition, 

```{math}
:label: eq225
\begin{align}
u(x,y) = f(x,y) \; \textrm{  for  } x^2+y^2 = R^2
\end{align}
```

We can first apply chain rule to find the diffusion in polar coordinate. (see {ref}`SpecialF`) 

```{math}
:label: eq226
\begin{align}
\nabla^2 U(r,\theta) = U_{rr}+\frac{1}{r}U_r+\frac{1}{r^2}U_{\theta\theta} \; \textrm{  for  } -\pi\leq \theta \leq \pi
\end{align}
```

For {eq}`eq226`, we know three basis, $1,r^n\cos(n\theta),r^n\sin(n\theta)$. form a completed solution space. i.e., 


```{math}
:label: eq227
\begin{align}
U(r,\theta) = \frac{1}{2}a_0+\sum_{n=1}^{\infty} (a_n r^{n}\cos(n\theta)+b_n r^{n}\sin(n\theta))
\end{align}
```

Using the boundary condition 


```{math}
:label: eq228
\begin{align}
U(R,\theta) = \frac{1}{2}a_0+\sum_{n=1}^{\infty} (a_n R^{n}\cos(n\theta)+b_n R^{n}\sin(n\theta))
\end{align}
```

$a_0$ and $a_n$ are simply the corresponding Fourier coefficients. 


```{math}
:label: eq229
\begin{align}
a_0 = \frac{1}{\pi} \int_{\-pi}^{\pi}f(\xi)d\xi
\end{align}
```

and 

```{math}
:label: eq230
\begin{align}
a_n &= \frac{1}{\pi R^n} \int_{\-pi}^{\pi}f(\xi)\cos(n\xi)d\xi \\
b_n &= \frac{1}{\pi R^n} \int_{\-pi}^{\pi}f(\xi)\sin(n\xi)d\xi \\
\end{align}
```
