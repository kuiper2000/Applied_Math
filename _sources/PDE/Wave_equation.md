(Wave)=
# Week 12:  Wave Equation  
## Problem Setup and Initial Conditions
Waves are ubiquitous and especially central to atmospheric sciences. Vibration on any material with restoring forces can generate waves. The existence of restoring force implies the existence of 2nd-order derivatives in time. One of the most well-known wave equation is Hooke's Law


```{math}
:label: eq175
\begin{align}
\frac{d^2x}{dt^2} = -kx
\end{align}
```

where the restoring force is equivalent to how much the material is stretched except with an opposite sign. This suggesting the force is always going against to the stretched direction. It's not hard to find the solutions share a form of 



```{math}
:label: eq176
\begin{align}
x = a_n \cos(\sqrt{k}t)+ b_n \sin(\sqrt{k}t)
\end{align}
```

The coefficients in {eq}`eq175` will be determined by the given initial values of $x$ and $x_t$. This is the general setup of wave solutions, where only the temporal structure is considered.  


## D'Alembert's Solutions, Characteristic Lines, and Dispersion Relationship 
Here we will consider a slightly more complicated case where the propagating over a space (i.e., we are observing wave at both space and time.). Therefore, the original equation can be rewritten as 


```{math}
:label: eq177
\begin{align}
y_{tt} = c^2 y_{xx} ;\ \textrm{for   } 0<x<L,t>0
\end{align}
```

with initial and boundary conditions 

```{math}
:label: eq178
\begin{align}
y(x,0) & = f(x), ;\ y_t(x,0) = g(x) \\ 
y(0,t) & = y(L,t) = 0 \\ 
\end{align}
```


Using separation of variable, one can expect that the solution shares a form of (leave this practice to readers)

```{math}
:label: eq179
\begin{align}
y_n(t) = X_n(x)T_n(t) = [a_n\cos(\frac{n\pi ct}{L})+b_n\sin(\frac{n\pi ct}{L})]\sin(\frac{n\pi x}{L})
\end{align}
```

In most cases, we can satisfy the initial condition with finite sum of $a_n$ and $b_n$. We therefore attempt a solutions 

```{math}
:label: eq180
\begin{align}
y(x,t) &= \sum_{n=1}^{\infty} y_n(x,t) \\
       &= \sum_{n=1}^{\infty}[a_n\cos(\frac{n\pi ct}{L})+b_n\sin(\frac{n\pi ct}{L})]\sin(\frac{n\pi x}{L})
\end{align}
```

using initial condition, we know 

```{math}
:label: eq181
\begin{align}
y(x,0) &= \sum_{n=1}^{\infty}a_n \sin(\frac{n\pi x}{L}) = f(x)
\end{align}
```

indicating 

```{math}
:label: eq182
\begin{align}
a_n = \frac{2}{L} \int_{0}^{L} f(\xi)\sin(\frac{n\pi \xi}{L})d\xi
\end{align}
```

and 


```{math}
:label: eq182
\begin{align}
\frac{n\pi c}{L}b_n = \frac{2}{L} \int_{0}^{L} g(\xi)\sin(\frac{n\pi \xi}{L})d\xi
\end{align}
```



:::{admonition} Example 1

Suppose a string with fixed ends at $x=0$ and $x=\pi$, with initial conditions of 

```{math}
:label: eq183
f(x) = \begin{cases} 
x     \textrm{   for } 0\leq x \leq \frac{\pi}{2}\\
\pi-x \textrm{   for } \frac{\pi}{2}< x \leq \pi\\
\end{cases} 
```
:::
