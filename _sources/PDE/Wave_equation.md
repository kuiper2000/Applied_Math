(Wave)=
# Week 12:  Wave Equation  
## Problem Setup and Initial Conditions
Waves are ubiquitous and especially central to atmospheric sciences. Vibration on any material with restoring forces can generate waves. The existence of restoring force implies the existence of 2nd-order derivatives in time. One of the most well-known wave equation is Hooke's Law


```{math}
:label: eq176
\begin{align}
\frac{d^2x}{dt^2} = -kx
\end{align}
```

where the restoring force is equivalent to how much the material is stretched except with an opposite sign. This suggesting the force is always going against to the stretched direction. It's not hard to find the solutions share a form of 



```{math}
:label: eq177
\begin{align}
x = a_n \cos(\sqrt{k}t)+ b_n \sin(\sqrt{k}t)
\end{align}
```

The coefficients in {eq}`eq175` will be determined by the given initial values of $x$ and $x_t$. This is the general setup of wave solutions, where only the temporal structure is considered.  


## Wave Solution with Space and Time Structures
Here we will consider a slightly more complicated case where the propagating over a space (i.e., we are observing wave at both space and time.). Therefore, the original equation can be rewritten as 


```{math}
:label: eq178
\begin{align}
y_{tt} = c^2 y_{xx} ;\ \textrm{for   } 0<x<L,t>0
\end{align}
```

with initial and boundary conditions 

```{math}
:label: eq179
\begin{align}
y(x,0) & = f(x), ;\ y_t(x,0) = g(x) \\ 
y(0,t) & = y(L,t) = 0 \\ 
\end{align}
```


Using separation of variable, one can expect that the solution shares a form of (leave this practice to readers)

```{math}
:label: eq180
\begin{align}
y_n(t) = X_n(x)T_n(t) = [a_n\cos(\frac{n\pi ct}{L})+b_n\sin(\frac{n\pi ct}{L})]\sin(\frac{n\pi x}{L})
\end{align}
```

In most cases, we can satisfy the initial condition with finite sum of $a_n$ and $b_n$. We therefore attempt a solutions 

```{math}
:label: eq181
\begin{align}
y(x,t) &= \sum_{n=1}^{\infty} y_n(x,t) \\
       &= \sum_{n=1}^{\infty}[a_n\cos(\frac{n\pi ct}{L})+b_n\sin(\frac{n\pi ct}{L})]\sin(\frac{n\pi x}{L})
\end{align}
```

using initial condition, we know 

```{math}
:label: eq182
\begin{align}
y(x,0) &= \sum_{n=1}^{\infty}a_n \sin(\frac{n\pi x}{L}) = f(x)
\end{align}
```

indicating 

```{math}
:label: eq183
\begin{align}
a_n = \frac{2}{L} \int_{0}^{L} f(\xi)\sin(\frac{n\pi \xi}{L})d\xi
\end{align}
```

and 


```{math}
:label: eq183
\begin{align}
\frac{n\pi c}{L}b_n = \frac{2}{L} \int_{0}^{L} g(\xi)\sin(\frac{n\pi \xi}{L})d\xi
\end{align}
```



:::{admonition} Example 1

Suppose a string with fixed ends at $x=0$ and $x=\pi$, with initial conditions of 

```{math}
f(x) = \begin{cases} 
x     \textrm{   for } 0\leq x \leq \frac{\pi}{2}\\
\pi-x \textrm{   for } \frac{\pi}{2}< x \leq \pi\\
\end{cases} 
```

and 

```{math}
g(x) = x (1-\cos(x))
```

We also assume $c=1$, which represents the propagating speed of waves 


From {eq}`eq181`, we know the solution has a form of 

```{math}
\begin{align}
y(x,t) = \sum_{n=1}^{\infty}[a_n\cos(nt)+b_n\sin(nt)]\sin(nx)
\end{align}
```

where 

```{math}
\begin{align}
a_n = \frac{2}{\pi} \int_{0}^{\pi} f(\xi)\sin(\frac{n\pi \xi}{L})d\xi  = \frac{4\sin(\frac{n\pi}{2})}{n^2\pi}
\end{align}
```

and 

```{math}
\begin{align}
b_n &= \frac{2}{n\pi c}\int_{0}^{\pi} g(\xi)\sin(\frac{n\pi \xi}{L})d\xi \\
    &= \frac{-2}{n^2}(-1)^{n} + \frac{2(-1)^{n}}{n^2-1} \\
    &= \frac{2}{n^2(n^2-1)}(-1)^{n}  
\end{align}
```
:::


Now considering a more complicated case, where we have external forcing. i.e.,

```{math}
:label: eq184
\begin{align} 
y_{tt} &= c^2y_{xx}+x;\ \textrm{for } 0 <x< L, t>0 \\
y(0,t) &=y(L,t) = 0, \\
y(x,0) &=x(L-x) = f(x) \\
y_t(x,0) &=x(1+\cos(\frac{\pi x}{L})) = g(x)
\end{align} 
```

Similar to how we deal with heat equation with external forcing, we can assume the solution shares a form of 


```{math}
:label: eq185
y(x,t) = Y(x,t) + \psi(x)
```

Substitute {eq}`eq185` back to {eq}`eq184`, we can have a homogeneous equation for $Y(x,t)$. 

```{math}
:label: eq186
Y_{tt} = c^2Y_{xx} 
```

and an 2nd-order ODE for the forcing term i.e., 

```{math}
:label: eq187
c^2\psi^{''} = -x 
```

Integrating {eq}`eq187` twice, we have 

```{math}
:label: eq188
\psi(x) = -\frac{x^3}{6c^2}+\alpha x +\beta
```

Now, look at boundary conditions, 

```{math}
:label: eq189
y(0,t) = 0 = c^2(Y(0,t)+\psi(0)) = c^2 Y(0,t)+c^\beta
```

To make {eq}`eq189` a homogeneous problem, we choose $beta=0$. 
For the other boundary condition, 

```{math}
:label: eq190
y(L,t) = 0 = c^2(Y(L,t)+\psi(L)) = c^2 Y(L,t)-\frac{L^3}{6c^2}+\alpha L
```


Let 

```{math}
:label: eq191
\alpha  = \frac{L^2}{6c^2}
```

we can have the same homogeneous equation at the second boundary. 


{eq}`eq191` indicates that 

```{math}
:label: eq192
\psi  = -\frac{x^3}{6c^2}+\frac{L^2}{6c^2}x = \frac{x}{6c^2}(L^2-x^2)
```

With the result above, we can rewrite the entire set of wave equation to 

```{math}
:label: eq184
\begin{align} 
Y_{tt} &= c^2Y_{xx};\ \textrm{for } 0 <x< L, t>0 \\
Y(0,t) &=Y(L,t) = 0, \\
Y(x,0) &=x(L-x)-\psi = x(L-x)- \frac{x}{6c^2}(L^2-x^2)\\
y_t(x,0) &=x(1+\cos(\frac{\pi x}{L})) = g(x)
\end{align} 
```

## D'Alembert's Solutions, Characteristic Lines, and Dispersion Relationship 