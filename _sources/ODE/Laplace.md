(Laplace)=
# Week 5 and 6: Laplace Transform
## What is a Laplace Transform and how it works? (s space and t or x space)
From the first two chapters, we know a lot of ODEs share a similar form of solution : i.e., $c_1e^{s_1t}+c_2e^{s_2t}$. When we substitute the solution back to the ODEs, we can usually subtract the exponential term and make the entire equation into an _algebra problem_. This is exactly the goal of _Laplace Transform_, which can be simplified as the flowchart below. 

```{math}
\begin{align}
\textrm{initial value problem (at t space)} & \Longrightarrow \textrm{algebra problem (at s space)} \\
& \Longrightarrow \textrm{solution of the algebra problem (at s space)} \\ 
& \Longrightarrow \textrm{solution of the initial value problem (at t space)}
\end{align}
```  

From the flowchart, we can find the role of _Laplace Transform_ is converting the entire equation from _t space_ to the _s space_ and simplify the problem. Mathematically, it can written as 

```{math}
:label: eq62
\mathcal{L}[f] = \int_0^{\infty} e^{-st}f(t)dt=F(s)
```  

One can consider {eq}`eq62` as a projection or mapping of $f(t)$ onto the s space. 


:::{admonition} Example 6
Let $f(t)=e^{at}$, with a non-zero a. The Laplace Transform of $f(t)$ is 

```{math}
\begin{align}
\mathcal{L}[f] & = \int^{\infty}_{0} e^{-st} e^{at} dt \\
               & = \int^{\infty}_{0} e^{(a-s)t} dt \\
               & = \frac{1}{(a-s)}e^{(a-s)t}|^{t=\infty}_{t=0} \\
\end{align}
```  

The Laplace Transform exist when $(a-s)<0$, which will lead to 

```{math}
F(s) = \frac{1}{s-a}
```

One should notice that $(a-s)<0$ is the necessary criteria for the existence of Laplace Transform of a function. That also means, if any function which has a growth rate $a$ greater than $s$, then there is no Laplace Transform of that given function. A very famous example is the function $e^{t^2}$, which doesn't have a Laplace Transform. 

:::

Table 1 are some of the widely used Laplace Transforms. Readers can follow the example above to find the Laplace Transform of the corresponding function. 

:::{admonition} Table 1
```{math}
\begin{align}
f(t) \Longrightarrow & F(s) \\
1    \Longrightarrow & \frac{1}{s} \\
t    \Longrightarrow & \frac{1}{s^2} \\
t^n  \Longrightarrow & \frac{n!}{s^{n+1}} \\
\frac{1}{\sqrt{t}}  \Longrightarrow & \sqrt{\frac{\pi}{s}} \\
e^{at} \Longrightarrow & \frac{1}{s-a} \\ 
te^{at} \Longrightarrow & \frac{1}{(s-a)^2} \\
t^{n}e^{at} \Longrightarrow & \frac{n!}{(s-a)^{n+1}} \\
\frac{1}{a-b}(e^{at}-e^{bt}) \Longrightarrow & \frac{1}{(s-a)(s-b)} \\
\sin(at) \Longrightarrow & \frac{a}{s^2+a^2} \\
\cos(at) \Longrightarrow & \frac{s}{s^2+a^2} \\
t\sin(at) \Longrightarrow & \frac{2as}{(s^2+a^2)^{2}} \\
t\cos(at) \Longrightarrow & \frac{s^2-a^2}{(s^2+a^2)^{2}} \\
e^{at}\sin(bt) \Longrightarrow &  \frac{b}{(s-a)^2+b^2} \\
e^{at}\cos(bt) \Longrightarrow &  \frac{(s-a)}{(s-a)^2+b^2} \\
\sinh(at) \Longrightarrow & \frac{a}{s^2-a^2} \\
\cosh(at) \Longrightarrow & \frac{s}{s^2-a^2} \\
\end{align}
```  
:::


One can find that Laplace Transform is linear, indicating any rules applied to a linear operator is also applied to Laplace Transform. i.e., $\mathcal{L}(f(x)+g(x))=\mathcal{L}(f)+\mathcal{L}(g)$ and $\mathcal{L}(af(x))=a\mathcal{L}(f(x))$. 


For the inverse transform, we can setup a problem, where we would like to find a function leads to 
```{math}
:label: eq63
1 = \int_0^{\infty} e^{-st} \mathcal{L}^{-1} dt
```

One should notice that $\mathcal{L}^{-1}$ needs to be a function of $s$ since we would like to transform $F(s)$ back to $t$ space. 
By observing {eq}`eq63`, we can find it with a form of

```{math}
:label: eq64
\mathcal{L}^{-1}[F(s)]= \frac{1}{2\pi i}\int_{\alpha-i\infty}^{\alpha+i\infty} e^{st} F(s)ds 
```

If readers still recall the complex roots of a constant coefficient ODE, it usually shares a solution of $e^{\alpha t}\cos(\beta t)+e^{\alpha t}\sin(\beta t)$. The $\alpha$ in the upper/lower bound of the integral is identical to the meaning of $e^{\alpha t}$, where we specify a basis for projection should have a certain growth rate $c$. On the other hand, $-i\infty$ and $i\infty$ indicates we loop over all possible wave number $\beta$, because sometimes we can have multiple $\beta$. Readers can find the proof of inverse Laplace Transform [here](https://math.stackexchange.com/questions/2926685/proof-of-inverse-laplace-transform). We will cover more details in the chapter of Fourier Transform. ({eq}`eq63` and {eq}`eq64` are actually a more general definition of Fourier and inverse Fourier Transforms. 


## Initial Value Problem 
To deal with an initial value problem, we first need to know the Laplace transform of a derivative. i.e., 

```{math}
:label: eq65
\mathcal{L}[f'](s) = \int_0^{\infty} e^{-st} f' dt 
```

Using integrate by part, we can find  

```{math}
:label: eq66
\int_0^{\infty} e^{-st} f' dt = -f e^{-st}|^{t=\infty}_{t=0} +s\int^{\infty}_{0} e^{-st}f(t)dt = sF(s)-f(t=0)
```

Similarly, we can go to higher order derivative

```{math}
:label: eq67
\begin{align}
\int_0^{\infty} e^{-st} f^{''} dt & = f^{'}e^{-st}|^{t=\infty}_{t=0} + s\int^{\infty}_{0} e^{-st}f^{'}(t)dt \\
                                  & = -f^{'}(t=0) + s(sF(s)-f(t=0)) \\
                                  & = s^2F(s)-f^{'}(t=0)-sf(t=0)
\end{align}
```

It is evident that if we would like to calculate the Laplace Transform of the 2nd-order derivative, we need the initial condition from the 0th-order and the 1st-order, which is consistent with what we learn in week 1. 


The readers can try to go to higher order derivative. It can be proofed that 
```{math}
:label: eq68
\begin{align}
\mathcal{L}[f^{(n)}](s) = s^{(n)}F(s)-s^{(n-1)}f(0)-s^{(n-2)}f^{(1)}(0)-\cdots-s^{0}f^{(n-1)}(0)
\end{align}
```
:::{admonition} Example 1
Solve the following ODE with Laplace Transform 

```{math}
y'-4y=1; \; y(0)=1
```
As we know $\mathcal{L}[y'] = sY(s)-y(0)$, $\mathcal{L}[y] = Y(s)$ and $\mathcal{L}[1] = \frac{1}{s}$. This leads to 

```{math}
(s-4)Y(s)-y(0)=\frac{1}{s}; \; y(0)=1
```

Rearrange the above equation, 

```{math}
Y(s)=\frac{1}{s(s-4)}+\frac{1}{s-4} 
```

From the table above, we know 

```{math}
y(t)=\mathcal{L}^{-1}[\frac{1}{s(s-4)}]+\mathcal{L}^{-1}[\frac{1}{s-4}] = \frac{1}{4}(e^{4t}-1)+e^{4t} 
```

(Readers can use integrating factor to have the same answer). 
:::


:::{admonition} Example 2
Solve $y^{''}+4y^{'}+3y=e^{t}; \; y(0)=0; \; y^{'}(0)=2$

First, we take Laplace Transform of each term
```{math}
\begin{align}
\mathcal{L}[y^{''}](s) &= s^2Y(s)-sy(0)-y^{'}(0) \\
4\mathcal{L}[y^{'}](s) &= 4sY(s)-4y(0) \\
3\mathcal{L}[y](s) &= 3Y(s) \\
\mathcal{L}[e^{t}](s) &= \frac{1}{s-1}
\end{align}
```

which yields 

```{math}
s^2Y(s)-2+4sY(s)+3Y(s) = \frac{1}{s-1}
```

or 

```{math}
(s+1)(s+4)Y(s) = \frac{2s-1}{s-1}
```

reorganize the equation, we can have 

```{math}
Y(s) = \frac{1}{8}\frac{1}{s-1}+\frac{3}{4}\frac{1}{s+1}-\frac{7}{8}\frac{1}{s+3}
```

From the table above, we know 
```{math}
y(t) = \frac{1}{8}e^{t}+\frac{3}{4}e^{-t}-\frac{7}{8}e^{-3t}
```

(hint: Readers can solve this problem by using undetermined coefficients)
:::

## Heaviside Function and Shifting Theorem 
Laplace Transform is a unique way of solving an initial value ODE. What it provides that previous methods do not, is its ability of handling discontinuity. i.e., piecewise continuous. 

First, let's define what is a piecewise continuity. It can be any of the following three conditions. 

:::{admonition} Piecewise Continuity 
(1) $f$ is continuous at all but perhaps finitely many points of $[a,b]$ \
(2) if $f$ is not continuous, then $f(t)$ has a finite limits as t approaches $t_0$ from the left and from the right. \
(3) $f(t)$ has a finite limit as t approaches a from the right and also t approaches b from the left. 
:::

Mathematically, this discontinuity can arise from filtering the data. In a real-world case, the localized SST forcing or convection can yield similar effects. One of the impulse forcing experiments, which has been widely used in atmospheric science, is called _Green's function_ approach. Dr. Yen-Ting Hwang in our department is an expert in this field. On the hand, when we study how the convection over the hurricane eyewall drive the transverse circulation, we use similar approach to handle the problem (more details will be covered in the Chapter of special function). 

To know the connection between impulse and Laplace Transform, we will walk through two shifting theorems first. We will start with the second shifting theorem of Laplace Transform since it's physically more intuitive. Inspecting {eq}`eq62`, we can already find the condition for discontinuity in Laplace Transform. The lower bound of the integral starts from 0, indicating any $f(t)$ with $s<0$ is discarded. We can consider function $f(t)$ is multiplied by a unit step, i.e., Heaviside function, which is defined as 

```{math}
:label: eq69
H(t) = \begin{cases}
0 & \text{if $t<0$} \\
1 & \text{if $t>=1$} 
\end{cases}
```

From the figure below, we can see that Heaviside function acts as an filter. 


```{figure} Heaviside1.png
---
name: FIG3
---
A schematic of Heaviside function with a jump condition at $t=0$
```

Now we can picture a new function, where we shift both $f(t)$ (the function which we would like to apply an Laplace Transform to) and Heaviside function slightly to the right...let's say by 3. Using a $\sin(t)$ as an example, the figure below shows the function before (upper panel) and after the shifting (lower panel). 


```{figure} Heaviside2.png
---
name: FIG4
---
The second shifting theorem 
```

We can find the entire function remains the same shape. If we take the Laplace Transform of the function above, 

```{math}
:label: eq70
\mathcal{L}[H(t-3)f(t-3)](s) = \int_{0}^{\infty} e^{-st} H(t-3)f(t-3) dt = \int_{3}^{\infty} e^{-st} H(t-3)f(t-3) dt  
```
Here we define a new t coordinate, $t_{\textrm{new}}$, where $t_{\textrm{new}}=t-3$, and use it to rewrite {eq}`eq70` 

```{math}
:label: eq71
\begin{align}
\int_{0}^{\infty} e^{-s(t_{\textrm{new}}+3)} H(t_{\textrm{new}})f(t_{\textrm{new}}) dt_{\textrm{new}} &= e^{-3s}\int_{0}^{\infty} e^{-st_{\textrm{new}}} H(t_{\textrm{new}})f(t_{\textrm{new}}) dt_{\textrm{new}} \\
&= e^{-3s}\mathcal{L}(f(t)) \\
& = e^{-3s}F(s) 
\end{align}
```

:::{admonition} Example 3
Determine the Laplace Transform of 

```{math}
g(t) = \begin{cases}
0 & \text{if $t<2$} \\
t^2+1 & \text{if $t>=2$} 
\end{cases}
```

We can first observe where the discontinuity happens. According the jump condition above, we know it happens at $t=2$. Therefore, we can define a new $t_{\textrm{new}}=t-2$ and rewrite the equation. 

```{math}
g(t) = \begin{cases}
0 & \text{if $t<2$} \\
(t-2)^2+4(t-2)+5 & \text{if $t>=2$} 
\end{cases}
```

or 

```{math}
g(t) = \begin{cases}
0 & \text{if $t_{\textrm{new}}<0$} \\
t_{\textrm{new}}^2+4t_{\textrm{new}}+5 & \text{if $t_{\textrm{new}}>=0$} 
\end{cases}
```

Then we can take Laplace Transform for the equation above, which gives us 

```{math}
e^{-2s}\mathcal{L}[t_{\textrm{new}}^2+4t_{\textrm{new}}+5] = e^{-2s} (\frac{1}{s^3}+\frac{4}{s^2}+\frac{5}{s})
```
:::

Now let's look at a real ODE problem. 

:::{admonition} Example 4
Consider 
```{math}
y^{''}+4y=f(t); \; y(t=0)=y^{'}(t=0)=0
```
where

```{math}
f(t) = \begin{cases}
0 & \text{if $t<3$} \\
t & \text{if $t>=3$} 
\end{cases}
```

Before we dive right in, readers can think about what methods other than Laplace Transform can be applied to this problem. By observing the equation, we can find it can solved by constant coefficient method plus the undetermined coefficients. It is evident that the general solution would be something like, 


```{math}
y(t) = c_1 \cos(2t)+c_2 \sin(2t) + At^2+Bt+C
```

Now let's define a new coordinate $t_{\textrm{new}}=t-3$. This enables us to rewrite the equation 

```{math}
y^{''}+4y=f(t_{\textrm{new}}); \; y(t_{\textrm{new}}=-3)=y^{'}(t_{\textrm{new}}=-3)=0
```

and 

```{math}
f(t) = \begin{cases}
0 & \text{if $t_{\textrm{new}}<0$} \\
t_{\textrm{new}}+3 & \text{if $t_{\textrm{new}}>=0$} 
\end{cases}
```

Using Laplace Transform, 

```{math}
\mathcal{L}[y^{''}+4y] = s^2F(s)-f^{'}(0)-sf(0)+4F(s) = (s^2+4)F(s) = \frac{1}{s^2}+\frac{3}{s} 
```

Rearrange the equation, 
```{math}
F(s) = \frac{1}{s^2}\frac{1}{s^2+4}+\frac{3}{s}\frac{1}{s^2+4} 
```

Assume we can decompose the equation above into different linear components 
```{math}
F(s) = \frac{A}{s}+\frac{B}{s^2}+\frac{Cs+D}{s^2+4}
```

(There is a small tip, when we apply such decomposition, we can set the numerator has an order small than the denominator.)

We will find $A=\frac{3}{4}$, $B=\frac{1}{4}$, $C=-\frac{3}{4}$, and $D=-\frac{1}{4}$ i.e.,

```{math}
F(s) = \frac{3}{4}\frac{1}{s}+\frac{1}{4}\frac{1}{s^2}-\frac{3}{4}\frac{s}{s^2+4}-\frac{1}{4}\frac{1}{s^2+4}
```

Take the inverse Transform will lead to the final form of solution 

```{math}
f(t_{\textrm{new}}) = \frac{3}{4}-\frac{3}{4}\cos(2t_{\textrm{new}})+\frac{1}{4}t_{\textrm{new}}-\frac{1}{8}\sin(2t_{\textrm{new}}) 
```

and 
```{math}
y(t) = \begin{cases}
0 & \text{if $t<3$} \\
\frac{1}{8}[2t-6\cos(2(t-3))-\sin(2(t-3))] & \text{if $t>=3$} 
\end{cases}
```

which is consistent with our first guess based on constant coefficient + undetermined coefficient methods.  
:::


In addition to the shift of function over the t space, we can also shift it over the s space (complex space), which is the first shifting theorem. Inspecting {eq}`eq62`, we can find to shift a function on s space, we need a function of $e^{at}$. 

If we multiply the original function $f(t)$ by $e^{at}$ and take the Laplace Transform, 

```{math}
:label: eq72
\mathcal{L}[e^{at}f](s) = \int_0^{\infty} e^{-st}e^{at}f(t)dt=F(s-a)
```  
we can see the function was shifted by a in s space. 

From equations {eq}`eq72` and {eq}`eq71`, one can find some similarity between the first and the second shifting theorems. 
```{math}
\begin{align}
\textrm{First shifting theorem} \Longrightarrow & \mathcal{L}[e^{at}f](s) = \int_0^{\infty} e^{-st}e^{at}f(t)dt=F(s-a) \\
\textrm{Second shifting theorem} \Longrightarrow & \mathcal{L}[H(t-a)f(t-a)](s) = \int_a^{\infty} e^{-st}f(t-a)dt=e^{-as}F(s) \\
\end{align}
```

If we would like to shift a function over physical space (t space), we need to multiply its Laplace transform (s space counterpart) by $e^{-as}$. Similarly, if we would like to shift a function over complex space (s space), we need to multiply the function at the physical space by $e^{at}$. 

:::{admonition} Example 5
Determine 
```{math}
\mathcal{L}^{-1}[\frac{4}{s^2+4s+20}](t)
```

we can first organize the function into a form of $e^{at}f(t)$. Let $s_{\textrm{new}}=s+2$

```{math}
\begin{align}
\mathcal{L}^{-1}[\frac{4}{(s+2)^2+16}](t) & = \mathcal{L}^{-1}[\mathcal{L}[\sin(4t)](s_{\textrm{new}})](t) \\
                                          & = \mathcal{L}^{-1}[\mathcal{L}[e^{-2t}\sin(4t)](s)](t) \\
                                          & = e^{-2t}\sin(4t)
\end{align}
```
:::


## Convolution 
Another widely used technique in atmospheric science, signal processing, and machine learning is called convolution. In atmospheric science, it has another famous nickname: running average or running mean. 

We can start from {eq}`eq69`. Here we take the difference of two Heaviside function $H(t-a)$ and $H(t-b)$ with $b>a$. This will create two jump conditions: 

```{math}
:label: eq73
g(t)=H(t-a)-H(t-b) = \begin{cases}
0 & \text{if $t<a$} \\
1 & \text{if $t>=a$ and $t<b$} \\
0 & \text{if $t>=b$ }
\end{cases}
```

{eq}`eq73` is visualized in the figure below. 

```{figure} Heaviside3.png
---
name: FIG5
---
A schematic of Heaviside function with two jump conditions at $t=3$ and $t=7$
```

Now we introduce a new operator, _convolution_
Mathematically, it is written as 

```{math}
:label: eq74
(f*g)(t)=\int_0^{\tau} f(t-\tau)g(\tau)d\tau
```

{eq}`eq74` indicates if we first divide $f(t)$ into infinite data points $f(t_i)$ and take the average for the first data point from $t_1-\frac{b-a}{2}$ to $t_1+\frac{b-a}{2}$. After the calculation, we repeat it all over again for $t_2$ and so and so forth until we loop over the entire $f(t)$. This can be visualized with the figure below. 

```{figure} Heaviside4.png
---
name: FIG6
---
A schematic of how convolution works when we apply a Heaviside function with a window size of 4 to a noisy data. (top): The Heaviside function used for calculation (middle): the raw time series, $f(t)$, and (bottom): $(g*f)(t)$
```

From the figure above, one can find that $(f*g)(t)$ has less wiggles than the original time series $f(x)$ indicating when we apply a mask with a constant value for convolution, it can smooth out the high frequency noise. We can also have different shape of $g(t)$ for convolution. (see later examples)

There are some very important characteristics for the Laplace Transform of a convolution, so-called _convolution theorem_. 

```{math}
:label: eq75
\mathcal{L}[(f*g)]=F(s)G(s)
```

and 

```{math}
:label: eq76
\mathcal{L}^{-1}[F(s)G(s)]=(f*g)(t)
```

{eq}`eq75` tells us the Laplace Transform of convolution of two functions equals the multiplication of each other's transforms. so does the inverse transform. 


:::{admonition} Example 6
Find 

```{math}
\mathcal{L}^{-1}[\frac{1}{s(s-4)^2}]
```

From the convolution theorem, we can assume 

```{math}
\begin{align}
F(s) & = \frac{1}{s} \\
G(s) & = \frac{1}{(s-4)^2} 
\end{align}
```

where 

```{math}
\begin{align}
\mathcal{L}^{-1}[F(s)] & = \mathcal{L}^{-1}[\frac{1}{s}]=1 \\
\mathcal{L}^{-1}[G(s)] & = \mathcal{L}^{-1}[\frac{1}{(s-4)^2}] = te^{4t} 
\end{align}
```

then 

```{math}
\begin{align}
\mathcal{L}^{-1}[F(s)G(s)] & = 1*te^{4t} = \int_{0}^{t} \tau e^{4\tau}d\tau \\
                           & = \frac{1}{4}te^{4t}-\frac{1}{16}e^{4t}+\frac{1}{16}
\end{align}
```
:::


:::{admonition} Example 7
Find the solution of 

```{math}
y^{''}-2y^{'}-8y=f(t); \; y(0)=1,y^{'}(0)=0 
```

Take the Laplace transform of both sides of the equation, we have 
```{math}
\begin{align}
\mathcal{L}[y^{''}](s) &= s^2Y(s)-sy(0)-y^{'}(0) \\
-2\mathcal{L}[y^{'}](s) &= -2sY(s)+2y(0) \\
8\mathcal{L}[y](s) &= 8Y(s) \\
\mathcal{L}[f(t)](s) &= F(s)
\end{align}
```

and solve for $Y(s)$

```{math}
Y(s) = \frac{s-2}{s^2-2s-8}+\frac{1}{s^2-2s-8}F(s)
```

fractionalize $Y(s)$, we have 

```{math}
Y(s) = \frac{1}{3}\frac{1}{s-4}+\frac{2}{3}\frac{1}{s+2}+\frac{1}{6}\frac{1}{s-4}F(s)-\frac{1}{6}\frac{1}{s+2}F(s)
```

Apply with the inverse transform, 

```{math}
y(t) = \frac{1}{3}e^{4t}+\frac{2}{3}e^{-2t}+\frac{1}{6}e^{4t}*f(t)-\frac{1}{6}e^{-2t}*f(t)
```

Readers can substitute $f(t)$ with any function of interest to get the final solution. 

:::


## Impulse, Dirac Delta Function and Green's Function approach  
Following {eq}`eq73`, if one makes the gap between $a$ and $b$ extremely small, let's say nearly 0. Then, the window in {eq}`eq73` will become an _impulse_. Here we define a function, _Dirac Delta Function_. 

```{math}
:label: eq77
\delta_{\epsilon}(t)=\frac{1}{\epsilon}[H(t)-H(t-\epsilon)]
```

when the duration goes to $0$, the height of shock goes to infinity. The Laplace transform of the Dirac Delta function also has a very special quantity. For example, if today we have a shifted pulse, 


```{math}
:label: eq78
\delta_{\epsilon}(t)=\frac{1}{\epsilon}[H(t-a)-H(t-a-\epsilon)] = \begin{cases}
0 & \text{if $t<a$} \\
\frac{1}{\epsilon} & \text{if $t>=a$ and $t<a+\epsilon$} \\
0 & \text{if $t>=a+\epsilon$ }
\end{cases}
```

and take the Laplace transform, 
```{math}
:label: eq79
\begin{align}
\mathcal{L}[\delta_{\epsilon}(t-a)](s)&=\frac{1}{\epsilon}[\frac{1}{s}e^{-as}-\frac{1}{s}e^{-(a+\epsilon)s}] \\
                                      &=\frac{e^{-as}(1-e^{-\epsilon s})}{\epsilon s} \\
                                      &=e^{-as}
\end{align}
```


if we don't have any shift, $\mathcal{L}[\delta_{\epsilon}(t)](s) = 1$. 


:::{admonition} Example 8
Find the solution for the following 

```{math}
y^{''}+2y^{'}+2y=\delta (t-3); \; y(0)=y^{'}(0)=0
```

Applying Laplace Transform to the differential equation, we have 

```{math}
\begin{align}
\mathcal{L}[y^{''}](s) &= s^2Y(s)-sy(0)-y^{'}(0) \\
2\mathcal{L}[y^{'}](s) &= 2sY(s)-2y(0) \\
2\mathcal{L}[y](s) &= 2Y(s) \\
\mathcal{L}[\delta (t-3)](s) &= e^{-3s}
\end{align}
```

Reorganize the equation, 


```{math}
Y(s)=\frac{1}{s^2+2s+2}e^{-3s}=\frac{1}{(s+1)^2+1}e^{-3s}
```

Then 
```{math}
y(t) = \mathcal{L}^{-1}[Y(s)]
```

if we focus on $\frac{1}{(s+1)^2+1}$ first and leverage the first shifting theorem 

```{math}
\mathcal{L}^{-1}[\frac{1}{(s+1)^2+1}] = e^{-t}\sin(t)
```

then substitute $e^{-3s}$ back to the equation and use the second shifting theorem, we can have the final form of solution 
```{math}
y(t) = H(t-3)e^{-(t-3)}\sin(t-3)
```

:::