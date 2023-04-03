(Heat)=
# Week 11:  Heat Equation  
Heat diffusion equation has important applications in engineering and Earth science. Especially when it comes to studying the Earth energy balance. In this chapter, we will consider a very simple 1D thermal diffusion problem. We will also demonstrate that we can always go to higher-dimension problem using _separation of variable_. 


## History and Formula 
The heat equation was first developed by Joseph Fourier in 1822 for the purpose of modeling how heat diffuses over a certain material. Let's start from a 1D heat equation, which can be written as 

```{math}
:label: eq141
\frac{\partial u}{\partial t} = k \frac{\partial ^2 u}{\partial x^2}
```

Considering the solutions are bounded in a stick with $0\leq x \leq L $