(2ndODE)=
# Week 3 and 4: 2nd-order ODE
The 2nd-order ODE is defined as an ordinary differential equation involving the second order derivative, perhaps the first but no third or higher order derivatives. In general, the 2nd-order ODE shares the following form, 

```{math}
:label: eq27
y"+p(x)y'+q(x)y=f(x)
```

where $p(x)$, $q(x)$ and $f(x)$ can be 0. When $f(x)=0$, we call it _homogeneous problem_. Different from the homogeneous problem in 1st-order ODE, it has another meaning here. In a differential equation, $f(x)$ can usually be considered as a forcing term. A lot of time, the forcing can be localized. Now you can imagine why $f(x)=0$ is called homogeneous. For example, if we have a warm SST associated with ENSO, what kind atmospheric response we will get? In this case, the forcing, $f(x)=0$, is the warm SST and the atmospheric response is $y$. Another example is, if today we have a wok with a heat source (e.g., stove fire) sitting right below it, what kind temperature response will we have? In this case, the fire is the external forcing $f(x)$ and the temperature is the response $y(x)$. 

Similar to the 1st-order, in a 2nd-order ODE problem, we need at least 2 initial conditions (boundary conditions, i.e., {eq}`eq28`) for $y$ and $y^{'}$ to find the particular solution. Otherwise, we can only find general solutions. 

```{math}
:label: eq28
\begin{align}
y(x_0) &= A \\
y'(x_0) &= B \\  
\end{align}
```
For the examples above, the boundary condition of the Earth atmosphere is how much heat was taken by the ocean or emitted to the space through radiation. Of course, all of these problems are way more complicated (we will talk a bit when at the chapter of special function).  

In the following section, we will walk through some fundamental characteristics/theorems of 2nd-order ODE. 


## Some Basic Theorems of 2nd-order ODE
### Theorem 1: The Existence of Solutions of an Initial Value Problem

For a 2nd-order ODE, like {eq}`eq27`, which has $p(x)$, $q(x)$ and $f(x)$ continuous on an open interval $I$ (the initial/boundary conditions are included in this interval), the solutions have some very important properties: 

:::{admonition} Properties
(1) The combinations of two solutions is still an solution. \
(2) A constant multiple of the solution is still an solution. (and is considered as the same solution) 
:::

One can easily find that it shares some similarity with linear function. Indeed, when {eq}`eq27` is linear, both properties above are also hold. 


### Theorem 2: The Uniqueness of Solutions of Initial Value Problems and Wronskian Test for Independence 
However, for a 2nd-order ODE, at least two solutions are required to completely solve the problem. i.e., we can use a set of solution like $y(x)=$ and $y'(x)=$, to represent the given system or using $y_1(x)=$ and $y_2(x)=$ as an alternative set of solution. One classic example is the harmonic oscillation. We can use location （displacement), x, and velocity, $dx/dt$, to describe the current state of a particle. An alternative way to describe harmonic oscillation is using circular motion, where we have two linear independent coordinates $\mathbf{x}=[x_i\hat{i},x_j\hat{j}]$

One key point is: $y_1$ and $y_2$ should be linearly independent. According to linear algebra, we know we can calculate the _cross product_ to test if two solutions are linearly independent. For $y_1$ we can consider its first and second coordinates as $y_1$ and $y^{'}_1$ and so does $y_2$. 

```{math}
:label: eq29
W[y_1,y_2](x) = \begin{vmatrix}
y_1(x)   & y_1(x) \\
y^{'}_1(x) & y^{'}_2(x)   
\end{vmatrix}
```
$W[y_1,y_2](x)$ is the _Wronskian_ of two functions, $y_1$ and $y_2$. If two solutions are independent, their cross product or Wronskian, i.e., the area spanned by two solutions, should be a non-zero value. We will walk through the details about Wronskian later, including how does it server as an additional constrain, which enables us to find the solution. 

### Theorem 3: General and Particular Solutions 
The third theorem is related to finding particular solutions in an ODE problem. 
Given that $p(x)$ and $q(x)$ are continuous in an open interval $I$ and $y_1$ and $y_2$ are two independent solutions. Then any general solution can be written as the linear combinations of the two solutions. i.e., 

```{math}
:label: eq30
y(x) = c_1y_1(x)+c_2y_2(x)
```

also, we have the initial/boundary values for $y$ 
```{math}
:label: eq31
\begin{align}
y(x_0) & = c_1y_1(x_0)+c_2y_2(x_0) = A \\
y^{'}(x_0) & = c_1y^{'}_1(x_0)+c_2y^{'}_2(x_0) = B
\end{align}
```

Using Cramer's rule, one can find 
```{math}
:label: eq32
\begin{align}
c_1 &= \frac{Ay'_2(x_0)-By_2(x_0)}{W[y_1,y_2](x_0)} \\ 
c_2 &= \frac{By_1(x_0)-Ay'_1(x_0)}{W[y_1,y_2](x_0)} \\ 
\end{align}
```

With the choice of $y_1$ and $y_2$, we can find a corresponding $y$. One should notice that as long as we decide what $y_1$ and $y_2$ are, we can find a unique solution of $y$. One can also find that if $W[y_1,y_2](x_0)$ equals 0, $Ay'_2(x_0)-By_2(x_0)=0$ and $By_1(x_0)-Ay'_1(x_0)=0$ are also required. At the end, this will lead to trivial solution. Therefore, a non-zero Wronskian is very important for finding an unique solution.


:::{admonition} Example 1
We know $e^x$ and $e^{2x}$ is a solution of 
```{math}
y"-3y'+2y=0
```
Let the general solution to be the linear combination of $e^x$ and $e^{2x}$

```{math}
y=c_1e^{x}+c_2e^{2x}
```

If today we have an initial condition, $y(0)=-1$ and $y'(0)=3$. We also know the Wronskian

```{math}
W[y_1,y_2](x) = \begin{vmatrix}
e^x & e^{2x} \\
e^x & 2e^{2x}   
\end{vmatrix} = e^{3x} = 1
```

Using {eq}`eq32`, we can find c_1 and c_2 


```{math}
\begin{align}
c_1 &= \frac{Ay'_2(x_0)-By_2(x_0)}{W[y_1,y_2](x_0)}=\frac{-2*2-3*1}{1}=-7 \\ 
c_2 &= \frac{By_1(x_0)-Ay'_1(x_0)}{W[y_1,y_2](x_0)}=\frac{3*1+2}{1}=5 \\ 
\end{align}
```

so the final solution will be $y=-7e^{x}+5e^{2x}$
:::


From the discussion above, we know how to find general solution along with Wronskian and initial conditions. Now, let's consider a nonhomogeneous case. Let $y_p$ to be any solution of the non-homogeneous equation {eq}`eq27`. We also have the initial condition of $y(x_0)=A$ and $y'(x_0)=B$. Then every particular solution can be written with the following form. 


```{math}
:label: eq33
Y(x) = c_1y_1(x)+c_2y_2(x)+y_p(x)
```

One can easily find that it's the combination of the homogeneous solution and the given particular solution. This can be proved by defining $Y(x)=y(x)+y_p(x)$, where $y(x)$ is the homogeneous solution we mentioned above, and substituting $Y(x)$ into {eq}`eq27`. We can find 


```{math}
:label: eq34
(y^{"}+y^{"}_p)+p(x)(y^{'}+y^{'}_p)+q(x)(y+y_p)=f(x)
```

and reorganize the equation 

```{math}
:label: eq35
(y^{"}+p(x)y^{'}+q(x)y)+(y^{"}_p+p(x)y^{'}_p+q(x)y_p)=f(x)
```

given that $y$ is the solution for homogeneous equation and $y_p$ is the solution for the nonhomogeneous equation, the entire equation of {eq}`eq34` also equals 0 indicating $Y$ is also an particular solution for {eq}`eq27`. 


:::{admonition} Example 2
For a given equation $y^{"}+4y=8x$, we know there is one particular solution, $2x$. Using theorem 3, we can find all possible particular solutions. 

First, we try to find the general solution of $y^{"}+4y$. We can assume the solution shares a form of $ce^{Ax}$ (we will talk a bit more in the section of solving the constant-homogeneous ODE). By substituting $ce^{Ax}$ into equation $y^{"}+4y=0$, we can find $A=\pm2i$, which gives as an general solution of $ce^{i2x}$

Then any particular solution can be written as
```{math}
y_p = ce^{i2x}+2x
```
:::


## Categories of 2nd-order ODE 
In this section, we will introduce some widely seen 2nd-order ODEs and the corresponding way to approach it, including (1) Constant Coefficient Homogeneous equation, (2)non-homogeneous ODE, (3) Euler Form, and (4) Series Solution 


:::{admonition} Forms of ODE 
-- Constant Coefficient Homogeneous equation: $y^{''}+Ay^{'}+By=0$ ($A$ and $B$ are constant)\
-- Nonhomogeneous ODE: $y^{''}+p(x)y'+q(x)y=f(x)$ \
-- Euler Form: $x^2y^{''}+xAy'+By=0$, where $A$ and $B$ are constant\
-- Series Solution : $y^{''}+p(x)y'+q(x)y=f(x)$, $y=\sum a_i x^{i}$
:::


### Constant Coefficient Homogeneous equation

For constant coefficient homogeneous equation, all of the equations share the same form: 

```{math}
:label: eq36
y^{''}+Ay^{'}+By=0
```
where $A$ and $B$ are constant.
One simple way to approach it is assuming that the solution shares a form of $Ce^{kx}$ where c and k are constant. By substituting the solution into {eq}`eq36` and solve for k, we can find $k=\frac{-A\pm\sqrt{A^2-4B}}{2}$. According to the root of $k$, we can further categorize the solutions into three different types (1) different but real roots, (2) repeated roots, and (3) complex roots. 

#### Case I: Different but Real Roots
This is the simplest case, where the solution is written as 

```{math}
:label: eq37
c_1 e^{\frac{-A+\sqrt{A^2-4B}}{2}x}+c_2 e^{\frac{-A-\sqrt{A^2-4B}}{2}x}
```
 To test if two solutions are linearly independent, we can use Wronskian.  

```{math}
W[y_1,y_2](x) = \begin{vmatrix}
e^{\frac{-A+\sqrt{A^2-4B}}{2}x}   & e^{\frac{-A-\sqrt{A^2-4B}}{2}x} \\
\frac{-A+\sqrt{A^2-4B}}{2}e^{\frac{-A+\sqrt{A^2-4B}}{2}x} & \frac{-A-\sqrt{A^2-4B}}{2}e^{\frac{-A-\sqrt{A^2-4B}}{2}x}   
\end{vmatrix} = -\sqrt{A^2-4B}e^{-Ax}
```

#### Case II: Repeated Roots 
For the case of repeated roots, if we apply the same Wronskian test, we will find that the two solutions are linearly dependent, suggesting we need to find another solution. A simple way to do that is assuming the second solution has a form of 

```{math}
:label: eq38
c_1 e^{\frac{-A}{2}x}+c_2 xe^{\frac{-A}{2}x}
```

Readers can use {eq}`eq38` and walk through Wronskian all over again and see if the answer changes this time. 


#### Case III: Complex Roots
The third case is similar to the first case where we have two different roots but both of them are complex (and they are complex conjugates). In this case, we have the same set of solution as the first cases. 

```{math}
:label: eq39
c_1 e^{(\alpha+i\beta)x}+c_2 e^{(\alpha-i\beta)x}
```

where $\alpha=-\frac{A}{2}$ and $\beta = \frac{\sqrt{4B-A^2}}{2}$. 
We can go one step further using _Euler identify_ to rewrite {eq}`eq39`

```{math}
:label: eq40
e^{i\beta x} = \cos(\beta x)+i\sin(\beta x)
```
this will lead to 
```{math}
:label: eq41
y = c_3e^{\alpha x}\cos(\beta x)+c_4e^{\alpha x}\sin(\beta x)
```
where $c_3=c_1+c_2$ and $c_4=i(c_1-c_2)$


### Nonhomogeneous Equation
To solve a nonhomogeneous 2nd-order ODE (i.e., $f(x)$ term exists), we can leverage two methods (1) _principle of superposition_ and (2) _Variation of Parameters_. 
Before we introduce _Principle of superposition_, we can first observe if $f(x)$ in {eq}`eq27` can be decomposed into the superposition of a few different functions. i.e., 

```{math}
:label: eq42
f(x) = f_1(x)+f_2(x)+\cdots+f_n(x)
```

and $Y_i$ is the particular solution of the equation $y^{''}+p(x)y^{'}+q(x)y=f_i(x)$. 

_principle of superposition_ then tells us, if today we have a function $y_p$ which equals to $\sum_i Y_i$, $y_p$ is also a particular solution of {eq}`eq27`. Based on this rule, we can decomposition $f(x)$ into different pieces which is easier for us to approach. 


The second method, _variations of parameters_, which is a method in combination of the way we solved homogeneous problem and Wronskian. (we will cover more details about how to solve more complicated homogeneous problem than a constant coefficient case but just leave it for convenience at the current stage) For example, if $y_1$ and $y_2$ are the solution of the homogeneous form of {eq}`eq27`. We can assume that the particular solution has a form 

```{math}
:label: eq43
y_p(x) = u_1(x)y_1(x) + u_2(x)y_2(x)
``` 

This gives us 

```{math}
:label: eq44
y^{'}_p(x) = u_1^{'}(x)y_1(x)+u_1(x)y^{'}_1(x) + u^{'}_2(x)y_2(x) + u_2(x)y^{'}_2(x)
``` 

For convenience, we will first assume we can find certain $u_1^{'}$ and $u^{'}_2(x)$ which makes 

```{math}
:label: eq45
u_1^{'}(x)y_1(x)+u^{'}_2(x)y_2(x) = 0 
``` 

Then we take another derivative of {eq}`eq44`, which gives us
```{math}
:label: eq46
y^{''}_p(x) = u_1^{'}(x)y^{'}_1(x)+u_1(x)y^{''}_1(x) + u^{'}_2(x)y^{'}_2(x) + u_2(x)y^{''}_2(x)
``` 

Then substitute {eq}`eq46` and {eq}`eq44` into {eq}`eq27`, we can have 
```{math}
:label: eq47
\begin{align}
&u_1^{'}(x)y^{'}_1(x)+u_1(x)y^{''}_1(x) + u^{'}_2(x)y^{'}_2(x) + u_2(x)y^{''}_2(x) \\
&+p(x)(u_1(x)y^{'}_1(x)+u_2(x)y^{'}_2(x)) \\ 
&+q(x)(u_1(x)y_1(x) + u_2(x)y_2(x)) =f(x)
\end{align}
``` 

Reorganize equation, we can find 

```{math}
:label: eq48
u_1^{'}(x)y^{'}_1(x)+u^{'}_2(x)y^{'}_2(x) = f(x)
``` 

Since both {eq}`eq45` and {eq}`eq48` should hold, this will lead to 

```{math}
:label: eq49
\begin{align}
u^{'}_1(x)&=-\frac{y_2(x)f(x)}{W(x)} \\
u^{'}_2(x)&=\frac{y_1(x)f(x)}{W(x)} \\
\end{align}
``` 
where $W(x)$ is the Wronskian of $y_1$ and $y_2$. 

Comparing {eq}`eq49` with {eq}`eq32`, it is evident that the _variation of parameters_ is using theorem 3 that we introduced above except that the case in theorem 3 has a constant coefficient. 


:::{admonition} Example 3
Find the general solution of the following equation

```{math}
y^{''}+4y=sec(x)
``` 
Before we solve the differential equation, we can first observe what type of differential equations it is? 
First, without $\sec(x)$, it's a constant coefficient of 2nd-order ODE, which has solution of $c_1e^{i2x}+c_2e^{i2x}$ (or can be written as $c_1 \cos(2x)+c_2 \sin(2x)$). 

Now, we can use _variation of parameters_ to find $c_1$ and $c_2$ as a function x. First, the Wronskian is 
```{math}
W[y_1,y_2](x) = \begin{vmatrix}
\cos(2x)  & \sin(2x) \\
-2\sin(2x) & 2\cos(2x)  
\end{vmatrix} = 2
```

using {eq}`eq49`, we can find 


```{math}
\begin{align}
u^{'}_1(x)&=-\frac{\sin(2x)\sec(x)}{2} \\
u^{'}_2(x)&=\frac{\cos(2x)\sec(x)}{2} \\
\end{align}
``` 

Given that $\sin(2x)=2\sin(x)\cos(x)$ and $\cos(2x)=2\cos(x)^2-1$, we have 
```{math}
\begin{align}
u^{'}_1(x)&=-2\sin(x) \\
u^{'}_2(x)&=\cos(x)-\frac{1}{2}\sec(x) \\
\end{align}
``` 

and integrate both sides of the equation above,

```{math}
\begin{align}
u_1(x)&=2\cos(x) \\
u_2(x)&=\sin(x)-\frac{1}{2}\ln |\sec(x)-tan(x)| \\
\end{align}
```

substitute $u_1$ and $u_2$ back to the general solution, we have 

```{math}
y_p = 2\cos(x)\cos(2x)+\sin(x)\sin(2x)-\sin(2x)\ln|\sec(x)-\tan(x)|
```

To find all possible general solution, we use theorem 3: the combination of solutions for homogeneous problem and nonhomogeneous problem will be the general solution for any nonhomogeneous problem. Therefore, we have the final solution 

```{math}
y = c_1\cos(2x)+c_2\sin(2x)+2\cos(x)\cos(2x)+y_p
```

:::


There is one special technique, _undetermined coefficients_, which is usually used along with _principle of superposition_ and theorem 3 (i.e., $y=y_g+y_p$, where $y_g$ is the general solution and $y_p$ is the particular solution). However, different from _variation of parameter_, it can only be applied to a constant coefficient ODE. The concept behind it is _guessing the $y_p$ based on $f(x)$_. 


For example, if today $f(x)$ evolves a trigonometric function, we can simply guess a particular like $c_1\cos(kx)+c_2\sin(kx)$. If today $f(x)$ shares a polynomial form, we can guess a solution of $Ax^2+B^x+C$. By substitute these particular solution, we can determine the coefficients such as $c_1$, $c_2$ or $A$, $B$ and $C$. The readers can see the similarity between _undetermined coefficients_ and _variation of parameters_. The main difference is that all of the coefficients here is _constant_ instead of a function of x. More details will be covered in chapter of _Laplace Transform_.  

:::{admonition} Example 4
Using undermined coefficient to solve the following ODE.
```{math}
y^{''}+4y=7e^{3x}
```
According to _undermined coefficients_, we can guess a particular solution of $Ae^{3x}$. Then substitute $Ae^{3x}$ into the original equation, we have

```{math}
9Ae^{3x}+4Ae^{3x}=7e^{3x}
``` 

indicating $A=\frac{7}{13}$. Now, we have one particular solution, $\frac{7}{13}e^{3x}$. 
Also, by inspecting the left hand side of the ODE, we know the solution of the homogeneous problem can be written as 

```{math}
y = c_1\cos(2x)+c_2\sin(2x)
``` 

Therefore, the completed solution is 
```{math}
y = c_1\cos(2x)+c_2\sin(2x)+y_p
``` 
:::


### Euler Form 
Euler form is a special case of 2nd-order ODE, where the 2nd-order derivative is multiplied by $x^2$. i.e.,

```{math}
:label: eq50
x^2y^{''}+xAy'+By=0
``` 
where $A$, and $B$ are constant. In such case, we can't apply the same approach that we used for solving a constant coefficient ODE. To transform the equation back to an constant coefficient ODE. We can first assume the solution shares a form of $y=x^{r}$. This results in 

```{math}
:label: eq51
\begin{align}
y^{''} & = r(r-1)x^{r-2} \\
y^{'} &  = rx^{r-1}
\end{align}
``` 

substitute {eq}`eq51` into {eq}`eq50`, we will have 

```{math}
:label: eq52
r^2+(A-1)r+B = 0
``` 

Solving for $r$ in {eq}`eq52`, we can find the answer for {eq}`eq50`. 
Similar to a constant coefficient ODE, we will have three different scenarios for $r$: (1) distinct but real roots, (2) repeated roots, and (3) complex roots. 

#### Case I: Distinct but real roots

This is the simplest case, where the solution can be directly written as 

```{math}
:label: eq53
y=c_1x^{r_1} +c_2x^{r_2}
``` 
where $c_1$ and $c_2$ are constant. $r_1$ and $r_2$ are two roots for {eq}`eq52`

#### Case II: Repeated Roots 
In the case where we have repeated roots for {eq}`eq52`, we can assume the second root has a form of 

```{math}
:label: eq54
y_2 = \ln(x)x^{r_1}
``` 

One can find the similarity in repeated-root case between Euler equation and constant coefficient ODE (see {eq}`eq38`).  This can be proofed by using _variation of parameters_, where we assume the second solution is 

```{math}
:label: eq55
y_2 = u(x)x^{r_1}
``` 

#### Case III: Distinct but Complex Roots  
When we have two complex root, the solution can be written as 
```{math}
:label: eq56
y=c_1x^{r_1} +c_2x^{r_2}
``` 
as those shown in case I. Here, we can go one step further given that $e^{a\ln(x)}=x^a$. This leads to 

```{math}
:label: eq57
y=c_1e^{r_1 \ln(x)} +c_2e^{r_2 \ln(x)}
``` 

Recall $r_1$ and $r_2$ need to be each others' complex conjugate. Therefore, we can rewrite them into $r_1=a+ib$ and $r_2=a-ib$, which results in the final solution of 

```{math}
:label: eq58
y=c_1e^{a}e^{i\ln(x)} +c_2e^{a}e^{-i\ln(x)}
``` 

or 

```{math}
:label: eq59
y=c_3e^{a}\cos(\ln(x)) +c_4e^{a}\sin(\ln(x))
``` 

### Series Solution
In the section, we will introduce one of the most powerful and general approaches for solving the ODE, _Series Solutions_.  
Let's start from what we learn above. For all of the previous ODE categories, we can find many of them share a solution form of $c_1e^{r_1x}+c_2e^{r_2x}$. Recall the Taylor expansion of an exponent (center around 0, i.e., $x_0$=0) can be written as 

```{math}
:label: eq60
e^{x} = 1+\frac{x^{1}}{1!}+\frac{x^2}{2!}+\cdots
``` 

This implies that most of the ODEs that we introduced above can be represented with a polynomial. Indeed, when we have _enough_ terms in the polynomial, we might be able to find a certain combination of these terms to mimic the solution of an ODE (see figure below). On the other hand, by decomposing the solution into different polynomial terms help us better interpret the results. For example, $y(x) =4e^{-x^3/3}+e^{-x^3/3}\int^{x}_0e^{\xi^3/3}d\xi$ is the solution of the function $y^{'}+x^2y=1$, where a closed form doesn't exist. It's apparent that understanding this solution is more challenging than analyzing a power series. 

```{figure} Exp_series.gif
---
name: FIG2
---
A schematic figure shows how the increase of power series to make a better prediction of a given function. [Wikipedia](https://en.wikipedia.org/wiki/Taylor_series#:~:text=The%20Taylor%20series%20of%20any%20polynomial%20is%20the%20polynomial%20itself.&text=The%20above%20expansion%20holds%20because,term%20in%20the%20infinite%20sum.) 
```


To find the series solution of an ODE, we can first assume it has a form of 

```{math}
:label: eq61
y=\sum a_n x^{n}
``` 

By substituting {eq}`eq61` into {eq}`eq27`, we can find the corresponding coefficient, $a_n$,of each term in the polynomial. 

:::{admonition} Example 5
Find the polynomial solution of the following 1st-order ODE. 

```{math}
y^{'}+2xy=\frac{1}{1-x}
``` 

First, assume the solution has a form of {eq}`eq61`. That also implies 

```{math}
y^{'} = \sum_{n=1} a_n (n)x^{n-1} 
``` 
Inspect the right hand side of ODE, it can be written as $\sum x^n$ (if the series converges). 
Reorganize the entire equation, we have 


```{math}
\sum_{n=1} a_n (n)x^{n-1} +  \sum_{n=0} 2a_n x^{n+1} = \sum_{n=0} x^n
``` 
or

```{math}
\sum_{n=0} a_{n+1} (n+1)x^{n} +  \sum_{n=1} 2a_{n-1} x^{n} = \sum_{n=0} x^n
``` 

For the term of $x^0$, we have 
```{math}
a_1 x^{0}  = x^{0}
``` 

For other higher order terms, we have a _recurrence relation_ 

```{math}
a_{n+1} (n+1)x^{n}+ 2a_{n-1}x^{n} = x^{n}
``` 

and drop the trivial solution $x=0$, we have a sets of final solution 

```{math}
\begin{align}
(n=0) &\; a_1-1=0 \\
(n=1) &\; a_2=\frac{1}{2}(1-2a_0) \\
...    &
\end{align}
``` 


:::
