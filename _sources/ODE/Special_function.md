(SpecialF)=
# Week 9 and 10:  Special Functions 
## Special Functions, Series Solutions and Recurrence Relations


Following the discussion of Sturm-Liouville form in the last chapter, we will go a step further to check the cases of  $r(x)\neq 1$,  $q(x)\neq 0$, and $p(x)\neq 0$. Indicating that we have inhomogeneous grid spacing. Some cases such as cylindrical coordinate, spherical harmonics and parabolic cylinder functions, all of their solutions correspond to a set of Sturm-Liouville equation with $r(x)\neq 1$,  $q(x)\neq 0$, and $p(x)\neq 0$. The reason why these special functions are important is that they have been widely used in geophysical fluid dynamics, quantum physics and electromagnetism. To solve these special functions, series solutions play the most important role where we will leverage it for finding the recurrence relation of each special function.

## Bessel Function 

Bessel function is the solution of Laplace equation on cylindrical coordinate. (see figure below). Considering a diffusion equation on a cylindrical coordinate,  


```{math}
:label: eq94
\nabla^2 \psi = \frac{1}{x}\frac{\partial}{\partial x} (x \frac{\partial \psi}{\partial x}) + \frac{1}{x^2} \frac{\partial^2 \psi}{\partial \phi^2}+\frac{\partial^2 \psi}{\partial z^2}
```


```{figure} cylindrical.png
---
name: FIG9
scale: 30%
---
Cylindrical coordinate. 
```

{eq}`eq94` is a partial differential equation. While we will cover more details in the later chapter, let us assume that {eq}`eq94` can be simplified by introducing separation of variables. i.e., 

```{math}
:label: eq95
\psi = X(x)\Phi(\phi)Z(z)
```
where $X(x)$, $\Phi(\phi)$ and $Z(z)$ are radial, tangential and vertical structure of the solution. 

By substituting {eq}`eq95` into equation {eq}`eq94`, we can have 

```{math}
:label: eq96
\begin{align}
\frac{d^2 }{dx^2}X+\frac{1}{x}\frac{d }{dx}X+(k^2-\frac{\mu^2}{x^2})X &= 0 \\
\frac{d^2 }{d\phi^2}\Phi+\mu^2\Phi &= 0 \\
\frac{d^2 }{dz^2}Z-k^2 Z &= 0 \\
\end{align}
``` 

where $k$ and $\mu$ are eigen values for vertical and tangential solutions respectively. Here, we focus on the radial structure given that the vertical and tangential solutions can be solved using Fourier series. 

By rearranging {eq}`eq96`, we can have 
```{math}
:label: eq97
\begin{align}
\frac{d^2 }{dx^2}X+\frac{1}{x}\frac{d }{dx}X+(1-\frac{\mu^2}{k^2 x^2})X &= 0 \\
\textrm{or} &  \\
x^2\frac{d^2 }{dx^2}X+x\frac{d }{dx}X+(x^2-\nu^2)X &= 0 \\
\end{align}
``` 

Observing {eq}`eq97`, one can find that this problem can be solved by combining Euler form with series solutions (we usually assume the solution has a form of $e^{At}$). By assuming the solution has a form of 


```{math}
:label: eq98
y= \sum_{n=0}^{\infty} c_n x^{n+r} 
```

which gives us 

```{math}
:label: eq99
y^{'}=\sum_{n=0}^{\infty} (n+r)c_n x^{n+r-1} 
```

and 

```{math}
:label: eq100
y^{''}=\sum_{n=0}^{\infty} (n+r)(n+r-1)c_n x^{n+r-2} 
```

Substituting {eq}`eq98`, {eq}`eq99`, and {eq}`eq100` into {eq}`eq97`, we can find 

```{math}
:label: eq101
\sum_{n=0}^{\infty} (n+r)(n+r-1)c_n x^{n+r} + \sum_{n=0}^{\infty} (n+r)c_n x^{n+r} + \sum_{n=0}^{\infty} c_n x^{n+r+2} - \sum_{n=0}^{\infty} \nu^2 c_n x^{n+r}=0
```

or 


```{math}
:label: eq102
\begin{align}
& c_0[r(r-1)+r-\nu^2]x^{r}+ \\ 
& c_1[(r+1)(r)+(r+1)-\nu^2]x^{r+1}+ \\
& \sum_{n=2}^{\infty} [[(n+r)(n+r-1)+(n+r)-\nu^2]c_n+c_{n-2}]x^{n+r}=0
\end{align}
```



