# Syllabus 

This Applied Math II course based on the book "Advanced Engineering Mathematics"{cite}`o2017advanced` will help you learn and (hopefully) equip with the basic skill sets in applied math, which are widely used in a variety of scientific research, such as physics, engineering and atmospheric science. There are three main topics in this class including (1) ordinary differential equations (2) Sturm-Liouville Theorem and (3) partial differential equation. Students are expected to learn how to (1) categorize the problem (2) find the solutions (both analytically and numerically) for 1st- and 2nd-order ODEs and even (3) formulate the questions mathematically. We will also walk through some very important examples in atmospheric science, which requires the techniques in this class for solving the problems. 


Prerequisites for this class includes Calculus, and Applied Math 1. 


## Course Outlines
__Part I: Ordinary Differential Equation__
* {ref}`1stODE`
	* What is an ODE ? (Different ODE forms)
    * 1st-order ODE
    * Categories of 1st-order ODE
        * Separation of variable
        * Linear ODE 
        * Potential Equation (for implicit ODE solution) 
        * Bernoulli form 
        * Riccati form 
* {ref}`2ndODE`
    * 2nd-order ODE
    * Some basic theorem of 2nd-order ODE 
        * Theorem 1: The Existence of Solutions of an Initial Value Problem
        * Theorem 2: The Uniqueness of Initial Value Problems and Wronskian Test for Independence 
        * Theorem 3: General and Particular Solutions 
    * Categories of 2nd-order ODE 
        * Homogeneous ODE with constant coefficients (general solution) 
        * Nonhomogeneous ODE (particular solution)  
        * Euler Form  
        * Series Solution: A more general approach for all of the problems above 

* {ref}`Laplace`
    * What is a Laplace Transform and how it works? (s space and t or x space)
    * Initial Value Problem
    * Heaviside Function and Shifting Theorem 
        * First Shifting Theorem 
        * Second Shifting Theorem
    * Convolution (Filtering) 
    * Impulse, Dirac Delta Function and Green's Function approach  

__Part II: Sturm-Liouville Theorem__
* {ref}`Sturm`
    * Sturm-Liouville Form 
    * Eigenfunction Expansions 
    * Fourier Analysis
        * Boundary Conditions and Types of Solutions 
        * Power Spectrum, Windows and Gibbs Phenomenon 

* {ref}`SpecialF`
    * Special Functions, Series Solutions and Recurrence Relations
    * Bessel Function 
    * Legendre Polynimial 
    * Hyperbolic Cylinder Function 

__Part III: Partial Differential Equation__
* {ref}`Heat`
    * History and Formula 
    * Forced Solutions 
    * Solutions on a Real String
    * Heat Diffusion on a 2D Plane 

* {ref}`Wave`
    * Problem Setup and Initial Conditions
    * Wave Solution with Space and Time Structures
    * Wave Motion in an Unbounded Medium
    * D’Alembert’s Solutions, Characteristic Lines, and Dispersion Relationship

* {ref}`Laplace`
    * Problem Setup and Initial Conditions
    * Dirichlet Problem for a Rectangle
    * Forced Solution-Poisson Equation
    * Green's Approach in Higher Dimension (Numerical) 


## Book
```{bibliography} references.bib
:filter: docname in docnames
```
