(Heat)=
# Week 11:  Heat Equation  
Heat diffusion equation has important applications in engineering and Earth science. Especially when it comes to studying the Earth energy balance. In this chapter, we will consider a very simple 1D thermal diffusion problem. We will also demonstrate that we can always go to higher-dimension problem using _separation of variable_. 


## History and Formula 
The heat equation was first developed by Joseph Fourier in 1822 for the purpose of modeling how heat diffuses over a certain material. Let's start from a 1D heat equation, which can be written as 

```{math}
:label: eq143
\frac{\partial u}{\partial t} = k \frac{\partial ^2 u}{\partial x^2}
```

Considering the solutions are bounded in a stick with $0\leq x \leq L$ and $k$ is so-called diffusion coefficient. If we closely observe {eq}`eq143`, we can find it is the combination of a 1st-order ODE in time and a 2nd-order ODE in space. Therefore, we need at least 1 initial condition and 2 boundary conditions, which usually have form of 

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

where the $\pm$ sign simply represents the direction of radiation. One can easily find these two kinds of boundary condition fall within the category of type I Sturm-Liouville problem. Therefore, it's also possible to have a mixed boundary condition, where 

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

Using _separation of variable_ of what we did in previous chapter, consists of attempting a solution of the form 

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

:::
