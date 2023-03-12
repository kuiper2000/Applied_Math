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
name: FIG10
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

where $J_{-N}=(-1)^{N}J_{N}$ (readers can try to prove this equation by assuming $N_{\textrm{new}}=N-\nu$). An alternative form of {eq}`eq115` to represent the second kind of Bessel function solution with $\nu=0$ is 

```{math}
:label: eq116
\begin{align}
Y_{0} = \frac{2}{\pi} [J_0(x)\ln(x)+\sum_{n=1}^{\infty}\frac{(-1)^{n+1}}{2^{2n}(n!)^2}\phi(n)x^{2n}] + \frac{2}{\pi}(\gamma-\ln(2))J_0(x)
\end{align}
``` 

where $\phi(n)=1+\frac{1}{2}+\frac{1}{3}+\cdots+\frac{1}{n}$ and $\gamma=\lim_{n\rightarrow\infty}(\phi-\ln(x))\sim 0.57721566\cdots$


Back to {eq}`eq115`, when $\nu\rightarrow N$ (N is an integer), readers can find the following relation according to L'Hôpital's rule (I will leave this practice to readers). 

```{math}
:label: eq117
\begin{align}
Y_{\nu} = \frac{J_{\nu}\cos(\nu \pi)-J_{-\nu}}{\sin(\pi \nu)} = \lim _{\nu\rightarrow N} [\frac{\partial J_{\nu}}{\partial \nu}-(-1)^N\frac{\partial J_{-\nu}}{\partial \nu}]
\end{align}
```


:::{admonition} Notes
To prove that $Y_{\pm\nu}=\lim _{\nu\rightarrow N} [\frac{\partial J_{\nu}}{\partial \nu}-(-1)^N\frac{\partial J_{-\nu}}{\partial \nu}]$  is a solution of Bessel function, we can take the derivative of the second equation of {eq}`eq97` with respect to $\nu$

```{math}
\frac{\partial}{\partial \nu}[x^2\frac{d^2 }{dx^2}X+x\frac{d }{dx}X+(x^2-\nu^2)X]= 0 \\
```

we have 

```{math}
x^2\frac{d^2 }{dx^2}\frac{\partial}{\partial \nu}X+x\frac{d }{dx}\frac{\partial}{\partial \nu}X+(x^2-\nu^2)\frac{\partial}{\partial \nu}X= 2\nu X \\
```

if we let $X=J_{\nu}-(-1)^N J_{-\nu}$, we can find $\frac{\partial}{\partial \nu}J_{\pm\nu}$ is the solution of Bessel function if and only if $2\nu X \rightarrow 0$. It's not hard to find $J_{\pm\nu}=0$, given $(-1)^NJ_{-\nu}= (-1)^{2N}J_{\nu}=J_{\nu}$. 

:::

Based on {eq}`eq114`, we can plot the solution of Bessel function of the first kind with different $\nu$. 

```julia
using Plots
using SpecialFunctions

x  = range(0.001, 5π, length=1000)
bessel_output = zeros(Float64, 5, 1000)
for nu  in 0:4
    sum_x = zero.(x)
    for j in 0:100
        sum_x = (-1).^j*1/(gamma(j+1)*gamma(j+nu+1))*(x/2).^(2*j+nu)+sum_x
    end
    bessel_output[nu+1,:] = sum_x
end

```


```{figure} Bessel_function.png
---
name: FIG11
---
An visualization of Bessel function solution of the first kind with $\nu=0-4$
```

The readers can also try to plot Bessel function solution of the second kind and visualize the result. 


### Modified Bessel function

There is one special kind of Bessel function where the coordinate is defined on imaginary coordinate i.e., $x=ix$. If we use $ix$ to replace $x$ (or $k$=$ik$) in the first equation of {eq}`eq96`, we can have so-called modified Bessel function. 

```{math}
:label: eq118
x^2\frac{d^2 }{dx^2}X+x\frac{d }{dx}X-(x^2+\nu^2)X = 0 
``` 

Modified Bessel function is a more widely used form in atmospheric science, where the wavy solution in radial direction (i.e., along $x$) becomes an exponent-like function where the solution is one-side bounded. i.e., $X=0$ when $x=0$ or $x=\infty$. On the other hand, the exponent-like function in $z$ direction becomes a Fourier-series problem (i.e., Sturm-Liouville with $p(x)=1$). In atmospheric science, it is called _vertical normal mode decomposition_. 

The first kind of solution of the modified Bessel function is related to the original Bessel function as follows: 

```{math}
:label: eq119
I_{\nu} \equiv (i)^{-\nu}J_{\nu}(ix)
``` 

and so is the second kind 

```{math}
:label: eq120
K_{\nu} \equiv \frac{\pi}{2}\frac{I_{-\nu}(x)-I_{\nu}(x)}{\sin(\nu \pi)};\ \textrm{with $\nu\rightarrow N$} 
``` 

The figure below visualizes the first kind solution of the modified Bessel function. 
```{figure} Modified_Bessel_function.png
---
name: FIG12
---
An visualization of modified Bessel function solution of the first kind with $\nu=0-4$
```

One can easily find the solution is bounded on one end but unbounded on the other end. 

A famous example of Bessel function is its application in solving the transverse circulation of tropical cyclone. For more details, readers can refer to paper Schubert and Hack 1982. 


### Generating Function and Recurrence Relations 
For all of these special function, we can find some kinds of generating functions which enable us to calculate the polynomial in a more efficient. For Bessel function, it's generating function is related to 

```{math}
:label: eq121
e^{\frac{x(t-1/t)}{2}}
``` 

is expanded in an infinite series in $t$, then the coefficient coefficient of $t^n$ is $J_{n}(x)$. i.e., 

```{math}
:label: eq122
e^{\frac{x(t-1/t)}{2}}=\sum_{-\infty}^{\infty}J_n(x)t^n
``` 

For this reason, this is called the generating function for integer order Bessel functions of the first kind. 
One can further rewrite {eq}`eq122` into 

```{math}
:label: eq123
\begin{align}
e^{\frac{x(t-1/t)}{2}} & =e^{\frac{xt}{2}}e^{\frac{-x}{2t}} \\
                       & =(\sum_{m=0}^{\infty}\frac{1}{m!}(\frac{xt}{2})^m)(\sum_{k=0}^{\infty}\frac{1}{k!}(-1)^{k}(\frac{x}{2t})^{k}) \\ 
                       & =(1+\frac{xt}{2}+\frac{1}{2!}\frac{x^2t^2}{2^2}+\frac{1}{3!}\frac{x^3t^3}{2^3}+\cdots)(1-\frac{x}{2t}+\frac{1}{2!}\frac{x^2}{2^2t^2}-\frac{1}{3!}\frac{x^3}{2^3t^3}+\cdots)
\end{align}
``` 

For each n, determine the coefficient of $t^n$, we can find the corresponding Bessel function. To illustrate, look for the coefficient of $t^4$ in this product. To retrieve the $t^4$ term, it can be the multiplication of $\frac{1}{4!}\frac{x^4t^4}{2^4}$ on the left with 1 on the right, or $\frac{1}{5!}\frac{x^5t^5}{2^5}$ on the left with $\frac{x}{2t}$ on the right...so on and so forth. 

The coefficient of $t^4$ in the product of the two series can be written as 


```{math}
:label: eq124
\begin{align}
& \frac{1}{2^4 4!}x^4-\frac{1}{2^6 5!}x^5+\frac{1}{2^8 6!}x^6+\cdots...
& =\sum_{n=0}^{\infty}\frac{(-1)^n}{2^{2n+4}n!(n+4)!}x^{2n+4}=J_4(x)
\end{align}
``` 


In addition to the generative function of Bessel function solutions, the Bessel functions with different order ex:$\nu$, $\nu-1$, and $\nu+1$ are related to each other in several important. 


```{math}
:label: eq125
\begin{align}
&\frac{d}{dx}(x^{\nu}J_{\nu}(x))  =x^{\nu}J_{\nu-1}(x) \\
&\frac{d}{dx}(x^{-\nu}J_{\nu}(x)) =-x^{-\nu}J_{\nu+1}(x) 
\end{align}
``` 

and 

```{math}
:label: eq126
\begin{align}
\frac{2\nu}{x}J_{\nu}(x)=J_{\nu+1}(x)+J_{\nu-1}(x)
\end{align}
``` 

I will leave the proof to the readers. From {eq}`eq126`, one can easily use lower order Bessel function to retrieve higher order Bessel function. 


## Legendre Polynomial 
Legendre is another very important and widely used function in atmospheric science, ex: the development of dynamical core, studying polar amplification, and solving Hadley circulation on a sphere. One can find that all of them are in a spherical coordinate. Indeed, one of the simplest way to decompose the phenomena on a sphere is using a spherical harmonics where the Fourier decomposition is applied to the zonal direction and the Legendre polynomial is applied to the meridional decomposition. 

### The origins and solutions 
The Legendre polynomial arises from the spherical harmonics.  

```{figure} Spherical_map.png
---
name: FIG13
scale: 30%
---
Spherical Coordinate
```

For Laplace equation on a spherical coordinate, it can be written. 
```{math}
:label: eq127
\begin{align}
\nabla^2 \psi = \frac{1}{r^2}\frac{\partial }{\partial r}(r^2 \frac{\partial}{\partial r} \psi)+\frac{1}{r^2\cos\theta}\frac{\partial}{\partial \theta}(\cos\theta \frac{\partial}{\partial \theta} \psi)+\frac{1}{r^2 \cos^2\theta} \frac{\partial^2}{\partial \phi^2}\psi
\end{align}
``` 

Again, if we applied separation of variable to {eq}`eq127` and let $\psi=R(r)\Theta(\theta)\Phi(\phi)$, we will have a set of ODEs and they are related to each others through eigenvalues $\lambda$ and $m$. 

```{math}
:label: eq128
\begin{align}
&\frac{1}{R}\frac{d}{r}r^2\frac{dR}{dr} = \lambda \\
&\frac{1}{\phi}\frac{d^2\Phi}{d\phi^2}  = -m^2 \\
&\lambda \cos^2\theta+\frac{\cos\theta}{\Theta}\frac{d}{d\theta}(\cos\theta\frac{d\Theta}{d\theta}) = m^2
\end{align}
``` 

Here we take the zonal average (average over $\phi$), the second equation becomes 0 and the third equation can be written as, 

```{math}
:label: eq129
\begin{align}
&\lambda \cos^2\theta+\frac{\cos\theta}{\Theta}\frac{d}{d\theta}(\cos\theta\frac{d\Theta}{d\theta}) = 0
\end{align}
``` 
From the equation above, we know the weighting function is $\cos\theta$. Since $d\sin\theta=\frac{1}{\cos\theta}d\theta$, we can further rewrite the above equation into 

```{math}
:label: eq130
\begin{align}
&\lambda \cos^2\theta \Theta+\cos^2\theta\frac{d}{d\sin\theta}((1-\sin^2\theta)\frac{d\Theta}{d\sin\theta}) = 0
\end{align}
``` 

Then subtract $\cos^2\theta$ from both terms and let $\sin\theta=x$, we will have  

```{math}
:label: eq131
\begin{align}
&\lambda \Theta+\frac{d}{dx}((1-x^2)\frac{d\Theta}{dx}) = (1-x^2)\frac{d^2}{dx^2}\Theta -2x\Theta +\lambda \Theta 
\end{align}
``` 

which satisfies the Sturm-Liouville form given that the solution is bounded at both poles. 

Following the typical procedure of solving Sturm-Liouville problem, assume the solution has a form of 

```{math}
:label: eq132
\begin{align}
\Theta=\sum_{n=0}^{\infty} a_nx^n
\end{align}
``` 

which gives us 

```{math}
:label: eq133
\begin{align}
\Theta^{'}&=\sum_{n=1}^{\infty} n a_nx^{n-1} \\
\Theta^{''}&=\sum_{n=2}^{\infty} n(n-1) a_nx^{n-2} \\
\end{align}
``` 

Substitute {eq}`eq132` and {eq}`eq133` into, 

```{math}
:label: eq134
\begin{align}
& \sum_{n=2}^{\infty} n(n-1) a_nx^{n-2} -\sum_{n=2}^{\infty} n(n-1) a_nx^{n}-2\sum_{n=1}^{\infty} n a_nx^{n}+\lambda \sum_{n=0}^{\infty} a_nx^n = 0\\
& \rightarrow \sum_{n=0}^{\infty} (n+2)(n+1) a_{n+2}x^{n}  -\sum_{n=2}^{\infty} n(n-1) a_nx^{n}-2\sum_{n=1}^{\infty} n a_nx^{n} +\lambda \sum_{n=0}^{\infty} a_nx^n = 0 
\end{align}
```


Examine each order, we have 
```{math}
:label: eq135
\begin{align}
\textrm{0th order: } & [(2)(1)a_{2}+\lambda a_0] x^0\\
\textrm{1st order: } & [(3)(2)a_{3}-2(1)a_1+\lambda a_1] x^1\\
\textrm{nth order (n$\geq$2): } & [(n+2)(n+1)a_{n+2}-n(n-1) a_n-2(n)a_n+\lambda a_n] x^n\\
\end{align}
```

we obtain
```{math}
:label: eq136
\begin{align}
\textrm{0th order: } & a_{2} =-\frac{\lambda}{2} a_0\\
\textrm{1st order: } & a_{3} =\frac{2-\lambda}{6} a_1\\
\textrm{nth order (n$\geq$2): } & a_{n+2}=\frac{n(n-1)-\lambda}{(n+2)(n+1)} a_n \\
\end{align}
``` 

One can find, as long as we determine $a_0$ and $a_1$, we determine all of the higher order terms. i.e., 
```{math}
:label: eq137
\begin{align}
\Theta =&\sum_{n=0}^{\infty} a_n x^n \\
       =&a_0[1+\frac{-\lambda}{2}x^2+\frac{-\lambda}{2}\frac{2(2-1)-\lambda}{(2+2)(2+1)}x^4+\cdots] + \\
        &a_1[x+\frac{2-\lambda}{6}x^3+\frac{2-\lambda}{6}\frac{3(3-1)-\lambda}{(3+2)(3+1)}x^5+\cdots]
\end{align}
```

We can find something interesting from the second and third lines of {eq}`eq137`. The entire equation is truncated when $\lambda = n(n-1)$. If we further set one of $a_0$ or $a_1$ equals 0, all of the even/odd terms vanish. This will form so-called _Legendre Polynomials_ where, 


```{math}
:label: eq136
\begin{align}
P_0(x) & =1  \\
P_1(x) & =x  \\
P_2(x) & =\frac{1}{2}(3x^2-1) \\ 
P_3(x) & =\frac{1}{2}(5x^3-3x) \\
P_4(x) & =\frac{1}{8}(35x^4-30x^2+3) \\
& \cdots 
\end{align}
``` 