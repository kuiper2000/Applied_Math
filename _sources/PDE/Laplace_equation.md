(Laplace)=
# Week 13:  Laplace Equation  
## Problem Setup and Initial Conditions
Laplace equation is nothing new to readers since we have introduced two related special functions in the chapter of {ref}`SpecialF`, Legendre polynomial and Bessel function. While those special functions form a complete set of solution of diffusion equation for the corresponding coordinate, there is an alternative way of solving diffusion problem when certain boundary conditions are satisfied. (i.e., Sturm-Liouville theorem is satisfied). 



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

If today, we setup a problem on a disk ( $\Omega \in x^2+y^2 \leq R^2$). It's a rare case when we use Cartesian coordinate but not so rare when we are using polar coordinate i.e., $x=r\cos(\theta)$ and $y=r\sin(\theta)$. One can consider Dirichlet's approach is an alternative way to solve the polar coordinate diffusion problem. 


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
a_0 = \frac{1}{\pi} \int_{-\pi}^{\pi}f(\xi)d\xi
\end{align}
```

and 

```{math}
:label: eq230
\begin{align}
a_n &= \frac{1}{\pi R^n} \int_{-\pi}^{\pi}f(\xi)\cos(n\xi)d\xi \\
b_n &= \frac{1}{\pi R^n} \int_{-\pi}^{\pi}f(\xi)\sin(n\xi)d\xi \\
\end{align}
```

Substituting {eq}`eq229` and {eq}`eq230` back to {eq}`eq227`, the complete equation can be written as 

```{math}
:label: eq231
\begin{align}
U(r,\theta) = & \frac{1}{2\pi}\int_{-\pi}^{\pi}f(\xi)d\xi + \\
& \frac{1}{\pi}\sum_{n=1}^{\infty} (\frac{r}{R})^n \int_{-\pi}^{\pi} (f(\xi)\cos(n\xi)\cos(n\theta)+f(\xi)\cos(n\xi)\sin(n\theta)) d\xi
\end{align}
```


:::{admonition} Example 1
Find the solution 

If the disk has radius $R=4$ and $U(4,\theta) = \theta^2$, then the solution by {eq}`eq231` is 

```{math}
\begin{align}
U(r,\theta) = & \frac{1}{2\pi}\int_{-\pi}^{\pi}\xi^2 d\xi + \\
& \frac{1}{\pi}\sum_{n=1}^{\infty} (\frac{r}{4})^n \int_{-\pi}^{\pi} (\xi^2\cos(n\xi)\cos(n\theta)+\xi^2\cos(n\xi)\sin(n\theta)) d\xi \\
= & \frac{1}{3}\pi^2 + \sum_{n=1}^{\infty}\frac{4(-1)^n}{n^2}
\end{align}
```

The integration used are 

```{math}
\frac{1}{\pi} \int_{-\pi}^{\pi} \xi^2 \cos(n\xi)d\xi = \frac{4(-1)^n}{n^2}
```

and 
```{math}
\frac{1}{\pi} \int_{-\pi}^{\pi} \xi^2 \sin(n\xi)d\xi = 0
```
:::



:::{admonition} Example 2
Solving the Dirichlet problem

```{math}
\begin{align}
& \nabla^2 u(x,y) = 0 \; \; \textrm{for x^2+y^2<9}
& u(x,y) = x^2y^2
\end{align}
```

We first convert this problem to polar coordinate, letting 

```{math}
\begin{align}
U(r,\theta) = u(r\cos(\theta),r\sin(\theta))
\end{align}
```

and the boundary condition is 

```{math}
\begin{align}
x=3\cos(\theta), y=3\sin(\theta)
\end{align}
```
so 

```{math}
\begin{align}
U(3,\theta) = 3^2\cos^2(\theta)\times 3^2\sin^2(\theta) = f(\theta)
\end{align}
```


Solve the problem using {eq}`eq231`
```{math}
\begin{align}
U(r,\theta) = & \frac{1}{2\pi}\int_{-\pi}^{\pi}81\cos^2(\xi)\sin^2(\xi) d\xi + \\
& \frac{1}{\pi}\sum_{n=1}^{\infty} (\frac{r}{4})^n \int_{-\pi}^{\pi} (81\cos^2(\xi)\sin^2(\xi) \cos(n\xi)\cos(n\theta)+81\cos^2(\xi)\sin^2(\xi)\cos(n\xi)\sin(n\theta)) d\xi \\
\end{align}
```

All we need is to find the integral of 

```{math}
\begin{align}
\int_{-\pi}^{\pi} 81\cos^{2}(\xi)\sin^{2}(\xi)d\xi = \frac{81\pi}{4} 
\end{align}
```

```{math}
\int^{\pi}_{-\pi} 81\cos^{2}(\xi)\sin^{2}(\xi)\cos(n\xi)d\xi = 
\begin{cases}
  0 & \text{if }n\neq 4 \\
  \frac{-81\pi}{8} & \text{if }n= 4
\end{cases}
```
,and 

```{math}
\int^{\pi}_{-\pi} 81\cos^{2}(\xi)\sin^{2}(\xi)\sin(n\xi)d\xi = 0 
```
:::


## Poisson Integral Formula

Poisson integral is a technique of writing a series solution into an integral form. 

In a special case, where $R=1$ and $U(1,\theta)=f(\theta)$

```{math}
:label: eq232
U(r,\theta) = \frac{1}{2\pi} \int_{-\pi}^{\pi}[1+2\sum_{n=1}^{\infty}r^n\cos(n(\xi-\theta))]f(\xi)d\xi
```

Define a _Poisson Kernel_

```{math}
:label: eq233
P(r,\xi) = \frac{1}{2\pi}[1+2\sum_{n=1}^{\infty} r^n \cos(n\xi)]
```

so 

```{math}
:label: eq234
U(r,\theta) =\int_{-\pi}^{\pi} P(r,\xi-\theta)f(\xi)d\xi
```

Thinking of a point inside the unit disk as a complex number and also having a polar coordinates $(r,\xi)$. Use the Euler's formula to write 

```{math}
:label: eq235
z = re^{i\xi} = r[\cos(\xi)+i\sin(\xi)]
```

then 

```{math}
:label: eq236
z^n = r^n e^{in\xi} = r^n [\cos(n\xi)+i\sin(n\xi)]
```

Substitute {eq}`eq236` back into {eq}`eq234`, we can find 

```{math}
:label: eq237
1+2\sum_{n=1}^{\infty}\cos(n\xi) & =\Re (1+2\sum_{n=1}^{\infty}z^n) \\
&= \Re(1+2\frac{z}{1-z}) \\
& =\Re(\frac{1+z}{1-z})  \\
& = \Re(\frac{1+re^{i\xi}}{1-re^{i\xi}})
```

To find the real part of $\frac{1+re^{i\xi}}{1-re^{i\xi}}$, we can apply the following algebraic manipulation, 

```{math}
:label: eq238
\frac{1+re^{i\xi}}{1-re^{i\xi}} & = \frac{1+re^{i\xi}}{1-re^{i\xi}}(\frac{1-re^{i\xi}}{1-re^{i\xi}}) \\
& = \frac{1-r^2+r(e^{i\xi}-e^{-i\xi})}{1+r^2-r(e^{i\xi}+e^{-i\xi})} \\
& = \frac{1-r^2+r(\cos(\xi)+i\sin(\xi)-\cos(\xi)+i\sin(\xi))}{1+r^2-r(\cos(\xi)+i\sin(\xi)+\cos(\xi)-i\sin(\xi))} \\
& = \frac{1-r^2+2ir\sin(\xi)}{1+r^2+2r\cos(\xi)} \\ 
& = \frac{1-r^2}{1+r^2-2r\cos(\xi)}+i \frac{2r\sin(\xi)}{1+r^2-2r\cos(\xi)}
``` 

and the real part is $\frac{1-r^2}{1+r^2-2r\cos(\xi)}$

This leads to the final form of _Poisson Integral_

```{math}
:label: eq239
U(r,\theta) = \frac{1}{2\pi}\int^{\pi}_{-\pi} \frac{1-r^2}{1+r^2-2r\cos(\xi-\theta)}f(\xi)d\xi
```


For a disk with radius R about the origin, a change of variables yields the Poisson solution. 

```{math}
:label: eq240
U(r,\theta) = \frac{1}{2\pi}\int^{\pi}_{-\pi} \frac{R^2-r^2}{R^2+r^2-2rR\cos(\xi-\theta)}f(\xi)d\xi
```

:::{admonition} Example 3

Find the Poisson Integral of example 1. 

With $R=4$ and $f(\theta) = \theta^2$

```{math}
U(r,\theta) = \frac{1}{2\pi}\int_{-\pi}^{\pi}\frac{16-r^2}{16+r^2-8r\cos(\xi-\theta)}\xi^2 d\xi = \frac{16-r^2}{2\pi}\int_{-\pi}^{\pi}
```
:::


## The Neumann Problem 
Different from Dirichlet problem where the boundary temperature (value) is directly given, Neumann problem describes the jump condition. Mathematically, a Dirichlet problem usually has a boundary condition like $u(s,t)=f(s)$ where $s$ is the boundary along a closed region $C$. While the Neumann usually has a form of 
$u_n(s,t)=g(s)$, where $u_n(s,t)$ is the normal derivative (perpendicular) on the boundary. (see figure below)

```{figure} Neumann.png
---
name: FIG14
scale: 20%
---
The Neumann boundary condition. 
```   

Here, we first define the unit normal vector and unit tangent vector. 

```{math}
:label: eq241
\frac{dx}{ds}\mathbf{i} +\frac{dy}{ds}\mathbf{j}
``` 

and the unit normal vector 

```{math}
:label: eq242
\mathbf{n}(s) = \frac{dy}{ds}\mathbf{i} -\frac{dx}{ds}\mathbf{j}
``` 

Then the normal derivative of a point on C can be written as 

```{math}
:label: eq243
\frac{\partial u}{\partial n} = \frac{\partial u}{\partial x}\frac{dy}{ds}-\frac{\partial u}{\partial y}\frac{dx}{ds}
``` 

which can be recognized as the dot product of the gradient u with the unit normal vector on C. 


Combining with what we learn in previous sector, a complete set of Neumann problem can be written as 

```{math}
:label: eq244
\nabla^2 u = 0 \; \; \textrm{for $(x,y)$ interior to D} \\ 
\frac{\partial u}{\partial n} = g(x,y) \; \; \textrm{for $(x,y)$ on C} \\ 
``` 


The set of equation above says, if we don't consider the heat (or some tracers) transport normal to the boundary, the heat flux integrated over the domain should equal 0 (which physically makes sense because the amount of heat loss for a given grid point is equivalent to the amount of heat received by other grid point).  



```{math}
:label: eq244
& \nabla^2 u = 0 \; \; \textrm{for $(x,y)$ interior to D} \\ 
& \frac{\partial u}{\partial n} = g(x,y)
``` 

With {eq}`eq244`, there is a simple condition must be satisfied. 

```{math}
:label: eq245
\oint_{C} g(x,y) ds &= \oint_{C} \frac{\partial u}{\partial n} ds \\
&=  \oint_{C} [\frac{\partial u}{\partial x}\frac{dy}{ds}-\frac{\partial u}{\partial y}\frac{dx}{ds}]ds \\
&=  \oint_{C} -u_y dx + u_x dy \\
&=  \iint_D(u_xx+u_yy) dA \textrm{by Green's theorem} 
```

if $u$ is harmonic on D (i.e., satisfies Sturm-Liouville problem), then $\nabla^2 u = 0$. Then we can conclude that $\oint_{C} g(x,y) ds=0$. 


:::{admonition} Example 4
Find if the following set of equation satisfies a Neumann problem

```{math}
\begin{align}
\nabla^2 u & = 0 \textrm{ for $0<x<1,0<y<1$} \\ 
u_n(x,y) & = 
\begin{cases}
0 \textrm{ on the lower, upper, and the left sides} \\ 
y^2 \textrm{ on the right}
\end{cases}
\end{align}
```
This indicates that 

```{math}
u_n(x,0) &= u_n(x,1) =u_n(0,y) = 0 \\
u_n(1,y) &= y^2\\
```

To test if the problem above is an Neumann problem, we can take the line integral along the boundary, which leads to 

```{math}
\oint_{C} g(x,y) ds = \int_{0}^{1}y^2 dy = \frac{1}{3} \neq 0 
``` 

Immediately, this problem is not a Neumann problem. 
:::



:::{admonition} Example 5
A very import characteristic of Neumann problem is...there is no unique solution. We will use this example to illustrate the result. 

```{math}
\begin{cases}
\nabla^2 u(x,y) &= 0 \textrm{for $0<x<a,0<y<b$} \\
\frac{\partial u}{\partial y}(x,0) &=      \;\;\; \frac{\partial u}{\partial y} (x,b)=0 \textrm{for $0<x<a$} \\
\frac{\partial u}{\partial x}(0,y) &= 0    \;\;\; \textrm{for $0<y<a$} \\
\frac{\partial u}{\partial x}(x,y) &= g(y) \;\;\; \textrm{for $0<y<a$} \\
\end{cases}
```

From the very beginning of this chapter, we know the above equations have solution of 

```{math}
u(x,y) = c+\sum_{n=1}^{\infty} c_n \cosh(\frac{n\pi x}{b})\cos(\frac{n\pi y}{b})
```

also, based on the derivative boundary 

```{math}
u_n(0,y) = \sum_{n=1}^{\infty} \frac{n\pi}{b} c_n \sinh(\frac{n\pi 0}{b})\cos(\frac{n\pi y}{b}) = g(y)
```

indicating $\frac{n\pi}{b} c_n \sinh(\frac{n\pi 0}{b})$ is the Fourier $\cos$ coefficient of the solution in y structure. i.e., 

```{math}
\frac{n\pi}{b}c_n\sinh(\frac{n\pi a}{b}) = \frac{2}{b}\int_{0}^{b} g(\xi) \cos(\frac{n\pi\xi}{b}) d\xi \; \; \; \textrm{or...}\\ 
u(x,y) = c+\sum_{n=1}^{\infty} c_n\cosh(\frac{n\pi x}{b})\cos(\frac{n\pi y}{b})
``` 

However, to satisfy Neumann condition, we also know, 

```{math}
c=\frac{1}{b} \int_{0}^{b} g(y)dy = 0 
```

i.e., the line integral of $g(y)$ along the boundary equals 0. One should notice that, $c$ is an arbitrary constant and it's not a Fourier $\cos$ coefficient for wave number 0 (otherwise, it will violate the Neumann conditions). This also implies that if $u$ is a solution, $u+c$ will be a constant as well (where $c$ is an arbitrary constant). 

:::


:::{admonition} Example 6
Here is an application of Neumann problem to a disk. 
Suppose $D$ is a disk with a radius $R$ about the origins. The boundary of $D$ is the circle $C$. In polar coordinate, the Neumann problem for $D$ is 

```{math}
\begin{cases}
\nabla^2 u(r,\theta) = 0 \;\;\; \textrm{for $\leq r < R,-\pi\leq\theta<\pi$} \\
\frac{\partial u}{\partial r}(R,\theta)=f(\theta) \;\;\; \textrm{for $-\pi\leq \theta\leq \pi$}
\end{cases}
```

A necessary condition for the existence of a solution is that the line integral of the normal derivative over the boundary is 0. 

```{math}
\int_{-\pi}^{\pi}f(\theta)d\theta = 0
```


:::