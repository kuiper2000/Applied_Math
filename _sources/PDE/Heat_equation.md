(Heat)=
# Week 11:  Heat Equation  
Heat diffusion equation has important applications in engineering and Earth science. Especially when it comes to studying the Earth energy balance. In this chapter, we will consider a very simple 1D thermal diffusion problem. We will also demonstrate that we can always go to higher-dimension problem using _separation of variable_. 


## History and Formula 
The heat equation was first developed by Joseph Fourier in 1822 for the purpose of modeling how heat diffuses over a certain material. (later on, the readers will see why the solution of heat diffusion can be approached with Fourier transform.) The heat equation follows a very simple energy conservation law. 

```{math}
:label: Heat
\textrm{heat change rate} = \textrm{heat flux through boundary to neighbor} + \textrm{forcing}
```

if we write it down in a mathematical form 

```{math}
:label: Heatmath
\frac{\partial c(x)\rho(x)u(x,t)}{\partial t} = -\frac{\partial q}{\partial x} + Q(x,t)
```

where $c(x)$ is the specific heat, $\rho(x)$ is the density of the stick and $q$ is the heat flux through the lateral (if other places are insulated) and $Q(x,t)$ is the local heating source, which can come from molecular process (i.e., friction between atom or radiation). The equation above can be visualized as follow. 

```{figure} Heat.png
---
name: FIG10
scale: 30%
---
Heat diffusion on a stick
```

One can notice that the heat flux is proportional to the temperature difference between the stick and its neighbor, which is the key process in determining the energy redistribution. When there is no temperature gradient (and no external heat source), there is no temperature change. Our goal is to understand how the initial temperature, $u(x,t)$, evolve as a function of $x$ and $t$. In some special cases, we can add convective process (advection) and radiation to the problem but we will keep it simple for now. 


{eq}`Heatmath` can be rewritten as  

```{math}
:label: eq143
\frac{\partial u}{\partial t} = k \frac{\partial ^2 u}{\partial x^2}
```

where $Q(x,t)$ is dropped for simplification and $k=\frac{\kappa}{c\rho}$, so-called diffusion coefficient.  

Considering the solutions are bounded in a stick with $0\leq x \leq L$. If we closely observe {eq}`eq143`, we can find it is the combination of a 1st-order ODE in time and a 2nd-order ODE in space. Therefore, we need at least 1 initial condition and 2 boundary conditions, which usually have form of 

```{math}
:label: eq144
\begin{align}
& u(0,t) = T_1, \; u(L,t) = T_2  \\
& u(x,0) = f(x)\; \textrm{for }  \; 0\leq x \leq L
\end{align}
```

The first two equations of {eq}`eq144` indicates that the temperature at both ends of the stick equals to their nearby environment. 
We can also use energy constraint where the temperature at both ends are not necessary the same as their environment. Therefore, the temperature difference will cause the heat radiating to its neighbor if the stick has higher temperature than its environment. This can be written as 

```{math}
:label: eq145
\begin{align}
& u_x(0,t) = A[u(0,t)-T], \; u_x(L,t) = -A[u(x,t)-T]  \\
& u(x,0) = f(x) \; \textrm{for }  \; 0\leq x \leq L
\end{align}
```   

where the $\pm$ sign simply represents the direction of radiation. One can easily find these two kinds of boundary condition fall within the category of type I Sturm-Liouville problem. It's also possible to have a mixed boundary condition, where 

```{math}
:label: eq146
\begin{align}
& u_x(0,t) = T_1, \; u_x(L,t) = -A[u(x,t)-T_2]  \\
& u(x,0) = f(x) \; \textrm{for }  \; 0\leq x \leq L
\end{align}
```   

or vice versa. 

Now let's take a look of a few cases. 



:::{admonition} Example 1
First consider the problem 

```{math}
\begin{align}
& u_t = ku_{xx} \; \textrm{for }  \; 0\leq x \leq L \\
& u(0,t) = u(L,t) = 0 \\ 
& u(x,0) = f(x) 
\end{align}
```   

Using _separation of variable_ of as what we did in previous chapter, consists of attempting a solution of the form 

```{math}
\begin{align}
u(x,t) = X(x)T(t)
\end{align}
```   

Substitute this into the differential equation to get 

```{math}
\begin{align}
XT^{'} = kX^{''}T
\end{align}
```   

which is equivalent to 

```{math}
\begin{align}
\frac{T^{'}}{kT} = \frac{X^{''}}{X}
\end{align}
```   

Observing the equation above, we can find the left hand side only depends on $t$ and the right hand side only depends on $x$. We also know that $X(x)$ and $T(t)$ vary independently. Therefore, only possibility exists: the ration of $T^{'}$ and $kT$ is a constant and so do $X^{''}$ and $X$. According to this, we can write down 


```{math}
\begin{align}
\frac{T^{'}}{kT} = \frac{X^{''}}{X} = -\lambda 
\end{align}
```  

One should notice that we can choose $lambda$ on the right in stead. However, to satisfy the Fourier solution is space, we can only have $-\lambda $. (readers can think about why?). 

Now we have two differential equations. For the spatial structure equation, 

```{math}
\begin{align}
& X^{''} + \lambda X = 0  \; \textrm{with }  X(0) = X(L) = 0 
\end{align}
```  

The first equation has a solution of Fourier $\sin$ function, 

```{math}
\begin{align}
& X_n(x) = \sin(\frac{n\pi x}{L}) \; \textrm{for } n=1,2,\cdots \\
& \textrm{where} \lambda =\frac{n^2\pi^2}{L^2}
\end{align}
```  

For the temporal structure equation, 
```{math}
\begin{align}
& T^{'} + k\frac{n^2\pi^2}{L^2} T = 0 
\end{align}
```  

This implies 
```{math}
\begin{align}
T = e^{-k\frac{n^2\pi^2}{L^2}t}
\end{align}
```  

Put two solutions together, we have 

```{math}
\begin{align}
u(x,t) = \sum_{n=1}^{\infty} b_n \sin(\frac{n\pi x}{L})e^{-k\frac{n^2\pi^2 k}{L^2}t}
\end{align}
```  

where 


```{math}
\begin{align}
u(x,0) = \sum_{n=1}^{\infty} b_n \sin(\frac{n\pi x}{L}) = f(x)
\end{align}
```  

$b_n$ is the Fourier $\sin$ coefficients for each Fourier mode, which can be derived by projecting $f(x)$ onto different Fourier basis

```{math}
b_n = \frac{2}{L}\int_{0}^{L}f(x)\sin(\frac{n\pi x}{L}) dx
```  

and the final solution is 

```{math}
\begin{align}
u(x,t) = \sum_{n=1}^{\infty} (\frac{2}{L}\int_{0}^{L}f(\xi)\sin(\frac{n\pi \xi}{L}) d\xi) \sin(\frac{n\pi x}{L})e^{-k\frac{n^2\pi^2 k}{L^2}t} 
\end{align}
```  

From the solution above, one can find it will gradually approach 0 when $t\rightarrow\infty$ because of $e^{-k\frac{n^2\pi^2}{L^2}t}$. This is consistent with physical intuition. Therefore, for a heat diffusion without the existence of external forcing, the solution will gradually be smoothed out. 
:::


:::{admonition} Example 2: Diffusion on an Insulated Stick
Now considering a case where we have temperature gradient on a insulated stick. The first half has temperature of $T$ and the second half has temperature $0$. 

```{math}
f(x) = \begin{cases} 
T \textrm{   for } 0\leq x \leq \frac{L}{2}\\
0 \textrm{   for } \frac{L}{2}< x \leq L\\
\end{cases}
```

Given the entire stick is insulated, we can expect that the equilibrium temperature will be $\frac{T}{2}$ i.e., half of the heat moves from the left to the right. In addition, the insulated stick will have Fourier $\cos$ function as solution given 0 radiation boundary condition. i.e., 

```{math}
u_x(0,t) = u_x(L,t) = 0
```

The equation above can be considered as that the environment and stick always have the same temperature which leads to 0 heat exchange between sticks and outside environment. For such condition, the solution has a form of Fourier $\cos$ function. Therefore, the solution can be written as

```{math}
\begin{align}
u(x,t) = \sum_{n=0}^{\infty} a_n \cos(\frac{n\pi x}{L})e^{-k\frac{n^2\pi^2k}{L^2}t}
\end{align}
```  

To get the Fourier coefficient, 

```{math}
a_n = \frac{2}{L} \int^{\frac{L}{2}}_{0} T \cos(\frac{n\pi  x}{L}) dx = \frac{2T}{n\pi}\sin(\frac{n\pi}{2})
```   

where 

```{math}
a_0 = \frac{T}{2}
```   

The solution can be written as 


```{math}
\begin{align}
u(x,t) = \frac{T}{2}+\frac{2T}{\pi}\sum_{n=1}^{\infty} \frac{1}{n}\sin(\frac{n\pi}{2}) \cos(\frac{n\pi x}{L})e^{-k\frac{n^2\pi^2k}{L^2}t}
\end{align}
```  

One can see that $k$ only appears in $e^{-k\frac{n^2\pi^2k}{L^2}t}$. If we choose a $k$ big enough, the signal will flatten out very quickly. 
:::


From two cases above, the readers can extend to some more complicated cases such as one side has a constant temperature and the other side is radiation boundary condition. I will leave the practice to the readers.   

:::{admonition} Example 3: Nonhomogeneous Problem
Solve the following case, 

```{math}
\begin{align}
& u_t = ku_{xx} \; \textrm{for }  \; 0\leq x \leq L \\
& u(0,t) = T_1 \\
& u(L,t) = T_2 \\ 
& u(x,0) = f(x) 
\end{align}
```
where either $T_1$ or $T_2$ is not 0 (and $T_1\neq T_2$). Such case can be considered as a nonhomogeneous case. Because the temperature on both ends are different, there is always heat moving from one side to the other side. Thus, this problem is equivalent to a problem with a constant forcing. The simplest way of doing that is assuming $u(x,t)$ has the following form 

```{math}
\begin{align}
u(x,t) = v(x,t) + \psi(x)
\end{align}
```  

where $\psi(x)$ can take care of the nonhomogeneous temperature of the stick and $v(x,t)$ will to lead a homogeneous problem. i.e., 

```{math}
\begin{align}
& v_t = kv_{xx}+ \psi_{xx}(x); \textrm{for }  \; 0\leq x \leq L \\
& v(0,t) = 0 \\
& v(L,t) = 0 \\ 
& v(x,0) = f(x) 
\end{align}
```

to make the set of equation above homogeneous, we have $\psi_{xx}(x)=0$. This implies $\psi(x) = Ax+B$. Using the boundary condition,  $\psi(0) = B = T_1$ and $\psi(L) = AL+T_1 = T_2$, we know $\psi(0) = (T_2-T_1)x/L+T_1$. To solve $v$, we can follow the same steps in example 1. 

:::


:::{admonition} Example 4: Inclusion of Convection and other Processes 
We can include convection and other processes, 
```{math}
\begin{align}
& u_t = (ku_{xx}+Au_x+Bu) \; \textrm{for }  \; 0\leq x \leq L \\
& u(0,t) = u(L,t) = 0 \\
& u(x,0) = f(x) 
\end{align}
```

We can first convert it to a form which we are familiar with


```{math}
\begin{align}
& v_t    = kv_{xx} \; \textrm{for }  \; 0\leq x \leq L \\
& v(0,t) = v(L,t) = 0 \\
& v(x,0) = g(x) 
\end{align}
```


Recall that solving the 1st/2nd-order ODE, we usually guess a solution of $e^{\alpha x+\beta t}$ to get a characteristic equation. However, it doesn't satisfy the boundary condition where $u(0,t) = u(L,t) = 0$. Therefore, we can use the same approach of variation of parameters and assume 

```{math}
\begin{align}
u(x,t) = v(x,t)e^{\alpha x+\beta t}
\end{align}
```

which will lead to 
```{math}
\begin{align}
& u_t(x,t)  = \beta  v(x,t)e^{\alpha x+\beta t}+v_t(x,t)e^{\alpha x+\beta t} \\
& u_x(x,t)  = \alpha v(x,t)e^{\alpha x+\beta t}+v_x(x,t)e^{\alpha x+\beta t} \\
& u_xx(x,t) = \alpha^2 v(x,t)e^{\alpha x+\beta t}+2\alpha v_x(x,t)e^{\alpha x+\beta t} + v_xx(x,t)e^{\alpha x+\beta t} \\
\end{align}
```

and 

```{math}
\begin{align}
& \beta  v(x,t)  = (k\alpha^2+kA\alpha-\beta+kB)v+(2k\alpha+kA)v_x+kv_{xx}
\end{align}
```

The equation above implies 

```{math}
\begin{align}
& k\alpha^2+kA\alpha-\beta+kB = 0 \\
& 2k\alpha+kA                 = 0
\end{align}
```

or 


```{math}
\begin{align}
& \alpha = -\frac{A}{2} \\
& \beta  = k(B-\frac{A^2}{4})
\end{align}
```

With the chosen $\alpha$ and $\beta$, we can reorganize the original equation to a solvable form. 
:::




## Forced Solutions
While {eq}`eq143` takes heat flux as the only process for redistributing heat, we can have additional heat sources, $F(x,t)$. Therefore, {eq}`eq143` is written as 


```{math}
:label: eq149
\begin{align}
& \frac{\partial u}{\partial t} = k \frac{\partial ^2 u}{\partial x^2}+F(x,t) \; \textrm{for }  \; 0\leq x \leq L \\
& u(0,t) = u(L,t) = 0 \\
& u(x,0) = f(x) 
\end{align}
```

It is easy to verify that _separation of variables_ can be applied to the questions above when $F(x,t)$ exist. However, we can do some scale analysis first... One can find that when $F(x,t)$ is small, we can expect the solution has a form of Fourier $\sin$ function. i.e., 

```{math}
:label: eq150
\begin{align}
u(x,t) = \sum_{n=0}^{\infty} b_n \sin(\frac{n\pi x}{L})e^{-n^2\pi^2 kt/L^2}
\end{align}
``` 

with the numbers the Fourier $\sin$ coefficients of $f(x)$. This suggests that, for the problem with the term of $F(x,t)$, attempt a solution 


```{math}
:label: eq151
\begin{align}
u(x,t) = \sum_{n=0}^{\infty} T_n(t) \sin(\frac{n\pi x}{L})
\end{align}
``` 

because the spatial structure is constrained by both ends, which ensures the solution is the linear combination of Fourier $\sin$ functions. According to {eq}`eq151`, one can find that $T_n(t)$ is simply the coefficients of Fourier $\sin$ functions. i.e.,  


```{math}
:label: eq152
\begin{align}
T_n(t)  = \frac{2}{L}\int_{0}^{L} u(\xi,t) \sin(\frac{n\pi \xi}{L})d\xi
\end{align}
``` 

Following similar veins, we take Fourier transform of the entire heat diffusion equation, which leads to 

```{math}
:label: eq153
\begin{align}
T_n^{'}(t)  & = \frac{2k}{L}\int_{0}^{L} u_{xx}(\xi,t) \sin(\frac{n\pi \xi}{L})d\xi + \frac{2}{L}\int_{0}^{L} u(\xi,t) \sin(\frac{n\pi \xi}{L})d\xi \\
& = \frac{2k}{L}\int_{0}^{L} u_{xx}(\xi,t) \sin(\frac{n\pi \xi}{L})d\xi + B_n(t)
\end{align}
``` 

where 

```{math}
:label: eq154
\begin{align}
B_n(t) = \frac{2}{L}\int_{0}^{L} F(\xi,t) \sin(\frac{n\pi \xi}{L})d\xi
\end{align}
``` 


Evaluate the last integral in {eq}`eq153` by carrying out two integrate by parts. 

```{math}
:label: eq155
\begin{align}
\frac{2k}{L}\int_{0}^{L} u_{xx}(\xi,t) \sin(\frac{n\pi \xi}{L})d\xi &= \frac{2k}{L}([u_{x}(\xi,t)\sin(\frac{n\pi \xi}{L})]_{0}^{L}-\frac{n\pi}{L}\int_{0}^{L}u_{x}\cos(\frac{n\pi \xi}{L})d\xi) \\
&=0-\frac{2k n\pi}{L^2}[[u\cos(\frac{n\pi \xi}{L})]_{0}^{L}+\frac{n\pi}{L}\int_{0}^{L}u\sin(\frac{n\pi \xi}{L})d\xi]
\end{align}
```





