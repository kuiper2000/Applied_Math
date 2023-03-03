(SpecialF)=
# Week 9 and 10:  Special Functions 
## Special Functions, Series Solutions and Recurrence Relations


Following the discussion of Sturm-Liouville form in the last chapter, we will go a step further to check the cases of  $r(x)\neq 1$,  $q(x)\neq 0$, and $p(x)\neq 0$. Indicating that we have inhomogeneous grid spacing. Some cases such as cylindrical coordinate, spherical harmonics and parabolic cylinder functions, all of their solutions correspond to a set of Sturm-Liouville equation with $r(x)\neq 1$,  $q(x)\neq 0$, and $p(x)\neq 0$. The reason why these special functions are important is that they have been widely used in geophysical fluid dynamics, quantum physics and electromagnetism. To solve these special functions, series solutions play the most important role where we will leverage it for finding the recurrence relation of each special function.

## Bessel Function 
### The origins and solutions 
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
J= \sum_{n=0}^{\infty} c_n x^{n+r} 
```

which gives us 

```{math}
:label: eq99
J^{'}=\sum_{n=0}^{\infty} (n+r)c_n x^{n+r-1} 
```

and 

```{math}
:label: eq100
J^{''}=\sum_{n=0}^{\infty} (n+r)(n+r-1)c_n x^{n+r-2} 
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

To get a non-trivial solution, we need $[r(r-1)+r-\nu^2]=0$, and $[[(n+r)(n+r-1)+(n+r)-\nu^2]c_n+c_{n-2}]=0$. (reader can think why we don't give the same constraint to $[(r+1)(r)+(r+1)-\nu^2]$). This forms a complete set a recurrence relationship. i.e., 

```{math}
:label: eq103
\begin{align}
& r^2-\nu^2=0 \\
& (2\nu+1)c_1 = 0 \\
& [(n+r)(n+r+1)+(n+r)-\nu^2]c_n+c_{n-2}=0
\end{align}
```

from the first equation, we know $r\pm \nu$. This also gives us $c_1=0$ according to the second equation of {eq}`eq103`. Similarly, solving for $c_n$, we have a recurrence for $n\geq2$  

```{math}
:label: eq104
c_n=-\frac{1}{n(n+2\nu)}c_{n-2}
```

On the other hand, because $c_1=0$, we have 

```{math}
:label: eq105
c_3=c_5=c_7=\cdots=0
```


according to {eq}`eq105`, we can find the entire equation is an even-indexed function (any odd-indexed term is dropped). This also indicates we can define a new index, $n=2N$, to rewrite {eq}`eq104`. 

```{math}
:label: eq106
c_{2N}=-\frac{1}{2N(2N+2\nu)}c_{2(N-1)} = (-1)\frac{1}{2^2 [N(N+\nu)]}c_{2(N-1)} 
```

based on {eq}`eq106`, we can find the relation between $2N$th-order term with the $0$th-order term. i.e., 



```{math}
:label: eq107
\begin{align}
c_{2N} & =-\frac{1}{2N(2N+2\nu)}c_{2(N-1)} \\
       & = (-1)^2\frac{1}{2^4 [N(N+\nu)][(N-1)(N-1+\nu)]}c_{2(N-2)} \\
       & = \cdots \\
       & =  (-1)^N\frac{1}{2^{2N} [N(N+\nu)][(N-1)(N-1+\nu)]\cdots[(1)(1+\nu)]}c_{2(N-N)} 
\end{align}
```

Therefore, the solution of Bessel function can be written as 


```{math}
:label: eq108
\begin{align}
J = c_0\sum_{N=0}^{\infty} (-1)^N\frac{1}{2^{2N} N!(1+\nu)(2+\nu)(3+\nu)\cdots(N+\nu)} x^{2N+\nu}
\end{align}
```


This is a nontrivial solution of Bessel function of order $\nu$ for any nonzero choice of constant $c_0$. This solution is usually expressed expressed in terms of the gamma function, which is defined as 


```{math}
:label: eq109
\begin{align}
\Gamma = \int_{0}^{\infty} t^{x-1}e^{-t}dt
\end{align}
```

One can find that it looks very similar to the Laplace Transform of $t^{x-1}$ indicating that we can solve the problem using integrate by part. 


```{math}
:label: eq110
\begin{align}
\Gamma(x) = \int_{0}^{\infty} t^{x-1}e^{-t}dt = - t^{x-1}e^{-t}|_{t=0}^{t=\infty} + (x-1) \int_{0}^{\infty} t^{x-2}e^{-t} dt 
\end{align}
```

The term $-t^{x-1}e^{-t}|_{t=0}^{t=\infty}$ is apparently 0, suggesting 

```{math}
:label: eq111
\begin{align}
\Gamma(x) = (x-1) \int_{0}^{\infty} t^{x-2}e^{-t} dt =(x-1) \Gamma (x-1)
\end{align}
```

If we go all the way to $\Gamma(1)$

```{math}
:label: eq112
\begin{align}
\Gamma(x) = (x-1)!\Gamma(1)
\end{align}
```

and given the fact $\Gamma(1)=1$, we find $\Gamma(x)=(x-1)!$. Similarly, one can prove $\Gamma(x+\nu+1) = (x+\nu)\cdots(3+\nu)(2+\nu)(1+\nu)\Gamma(\nu+1)$. Substitute this relation back to {eq}`eq108`, we have 


```{math}
:label: eq113
\begin{align}
J = c_0\sum_{N=0}^{\infty} (-1)^N\frac{\Gamma(\nu+1)}{2^{2N} N!\Gamma(n+\nu+1)} x^{2N+\nu}
\end{align}
```

for homogeneous problem, $c_0$ can be an arbitrary (but nonzero) constant. If we choose, $c_0=\frac{1}{2^{\nu}\Gamma(1+\nu)}$, then the entire equation can be written as 

```{math}
:label: eq114
\begin{align}
J = \sum_{N=0}^{\infty} (-1)^N\frac{1}{N!\Gamma(n+\nu+1)} (\frac{x}{2})^{2N+\nu}
\end{align}
```


{eq}`eq114` is so-called _Bessel function of the first kind of order $\nu$_

But don't forget, the Bessel function is a 2nd-order ODE. This indicates we have another set of solution which is linearly independent of {eq}`eq114`. One such solution is the _Bessel function of the second kind_ of order $\nu$ (also known as Weber function), which is denoted $Y_{\nu}$. It's related to the first kind as follows: 


```{math}
:label: eq115
\begin{align}
Y_{\nu} = \frac{J_{\nu}\cos(\nu \pi)-J_{-\nu}}{\sin(\pi \nu)}
\end{align}
```

where $J_{-n}=(-1)^{n}J_{n}$. An alternative form of {eq}`eq115` to represent the second kind of Bessel function solution with $\nu=0$ is 

```{math}
:label: eq116
\begin{align}
Y_{0} = \frac{2}{\pi} [J_0(x)\ln(x)+\sum_{n=1}^{\infty}\frac{(-1)^{n+1}}{2^{2n}(n!)^2}\phi(n)x^{2n}] + \frac{2}{\pi}(\gamma-\ln(2))J_0(x)
\end{align}
``` 

where $\phi(n)=1+\frac{1}{2}+\frac{1}{3}+\cdots+\frac{1}{n}$ and $\gamma=\lim_{n\rightarrow\infty}(\phi-\ln(x))\sim 0.57721566\cdots$


Recall that when we deal with Euler differential equation with repeated roots, we multiply the first root by $\ln (x)$ to find the second root. From {eq}`eq116`, we can find similar approach has been applied to find the solution of second kind of Bessel function solution. 

By inspecting {eq}`eq115`, we know it has a singularity at $x\rightarrow 0$, which is a rare condition in geophysical dynamics. 

