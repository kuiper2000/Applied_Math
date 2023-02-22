(SpecialF)=
# Week 9 and 10:  Special Functions 
## Special Functions, Series Solutions and Recurrence Relations


Following the discussion of Sturm-Liouville form in the last chapter, we will go a step further to check the cases of  $r(x)\neq 1$,  $q(x)\neq 0$, and $p(x)\neq 0$. Indicating that we have inhomogeneous grid spacing. Some cases such as cylindrical coordinate, spherical harmonics and parabolic cylinder functions, all of their solutions correspond to a set of Sturm-Liouville equation with $r(x)\neq 1$,  $q(x)\neq 0$, and $p(x)\neq 0$. The reason why these special functions are important is that they have been widely used in geophysical fluid dynamics, quantum physics and electromagnetism. To solve these special functions, series solutions play the most important role where we will leverage it for finding the recurrence relation of each special function.

## Bessel Function 

Bessel function is the solution of Laplace equation on cylindrical coordinate. (or one can consider it as some kinds of diffusion problem on cylindrical cylindrical). Considering a diffusion equation on a cylindrical coordinate,  


```{math}
:label: eq94
\nabla^2 \psi = \frac{1}{x}\frac{\partial}{\partial x} (x \frac{\partial \psi}{\partial x}) + \frac{1}{r^2} \frac{\partial^2 \psi}{\partial \phi^2}+\frac{\partial^2 \psi}{\partial z^2}
```
