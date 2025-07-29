# SBP Operators

how do SPB Operators work... 

# SBP Operators

## 1. Introduction & Overview

Summation-by-Parts (SBP) operators are a powerful class of numerical methods that provide a systematic way to construct stable and conservative discretizations of partial differential equations. In DispersiveShallowWater.jl, SBP operators form the backbone of all spatial discretizations, ensuring that important physical properties like mass and energy conservation are preserved at the discrete level.

The key insight behind SBP operators is that they **mimic integration by parts discretely**. Just as integration by parts is fundamental to many stability proofs and conservation laws in continuous analysis, SBP operators provide the discrete analog that allows us to transfer these important properties to the numerical scheme.

### Why SBP Operators Matter

For dispersive shallow water equations, maintaining physical properties like mass and energy conservation is crucial for:
- **Long-time stability**: Simulations remain stable over extended time periods
- **Physical accuracy**: The discrete solution respects fundamental physical laws
- **Robustness**: Methods are less sensitive to parameter choices and grid resolution

SBP operators achieve this by providing:
- **Provable stability**: Mathematical guarantees about the behavior of the numerical method
- **Exact conservation**: Discrete conservation laws that hold to machine precision
- **Flexibility**: A unified framework that encompasses finite differences, finite elements, and spectral methods

### Historical Context

SBP operators were originally developed for finite difference methods to provide the same stability guarantees that finite element methods naturally possess through their variational formulation. However, the framework has since been extended to encompass a much broader class of methods, making it a unifying principle in numerical analysis.

## 2. Mathematical Foundation

### The Core SBP Property

The fundamental property that defines a first-derivative SBP operator is:

```math
MD + D^T M = t_R t_R^T - t_L t_L^T
```

where:
- ``D`` is the derivative operator (matrix)
- ``M`` is the symmetric, positive definite mass matrix
- ``t_L = (1,0,...,0)^T`` and ``t_R = (0,...,0,1)^T`` extract boundary values

This property is the discrete analog of integration by parts:

```math
\begin{aligned}
\underbrace{ \boldsymbol{u}^T M D \boldsymbol{v} + \boldsymbol{u}^T D^T M \boldsymbol{v} }_{\displaystyle \approx \int_{x_{\min}}^{x_{\max}} u\, (\partial_x v) + \int_{x_{\min}}^{x_{\max}} (\partial_x u)\, v } 
&= 
\underbrace{ \boldsymbol{u}^T \boldsymbol{t}_R \boldsymbol{t}_R^T \boldsymbol{v} - \boldsymbol{u}^T \boldsymbol{t}_L \boldsymbol{t}_L^T \boldsymbol{v} }_{\displaystyle \approx u(x_{\max})\,v(x_{\max}) - u(x_{\min})\,v(x_{\min}) }.
\end{aligned}
```

where for periodic SBP operators their property naturally simplifies to
```math
MD + D^T M = 0.
```

### Understanding the Mass Matrix

The mass matrix ``M`` approximates the continuous inner product:

```math
\langle u, v \rangle_M = u^T M v \approx \int_{x_{min}}^{x_{max}} u(x) v(x) dx
```

For the approximation to be meaningful, we require:
```math
1^T M 1 = x_{max} - x_{min}
```

This ensures that the discrete inner product correctly integrates constants.

### Understanding the Derivative Operator

The derivative operator ``D`` approximates the spatial derivative of a function on a discrete grid. That is, for a smooth function ``u(x)``, we want the discrete expression ``D \boldsymbol{u}`` to approximate ``\partial_x u(x)`` at the grid points.

A first-derivative SBP operator is said to be **$p$-th order accurate** if it satisfies:

```math
D\, \boldsymbol{x}^k = k\, \boldsymbol{x}^{k-1}, \quad \text{for all } k = 0, 1, ..., p,
```

where ``\boldsymbol{x}^k = (x_1^k, \dots, x_N^k)^T`` and the operations are performed pointwise. This means the discrete operator correctly differentiates monomials up to degree $p$. In particular, **consistency** requires at least zeroth-order accuracy, i.e.,

```math
D\, \boldsymbol{1} = \boldsymbol{0}.
```

However, approximating derivatives near boundaries poses a challenge: interior points can use standard finite difference stencils, but near the edges (e.g., at `x_1` or `x_N`), one must use specially designed one-sided approximations that still preserve stability and accuracy. SBP operators handle this via their defining identity:

```math
MD + D^T M = t_R t_R^T - t_L t_L^T,
```

which ensures that the discrete scheme mimics **integration by parts**, and hence inherits stability properties from the continuous problem.

In **periodic domains**, the boundary terms vanish, and we recover:

```math
MD + D^T M = 0,
```

which corresponds to skew-symmetry with respect to the discrete inner product.

Thus, the SBP derivative operator plays a crucial role in ensuring both **accuracy** and **energy stability** of the numerical method. Its structure is deliberately crafted to balance interior approximations with boundary treatments in a way that mirrors the analytical tools used in PDE theory.




### Discrete Integration by Parts in Action

To see how this works in practice, consider computing ``u^T M D v``:

```math
u^T M D v + u^T D^T M v = u^T t_R t_R^T v - u^T t_L t_L^T v = u_N v_N - u_1 v_1
```

This shows that:
```math
\langle u, Dv \rangle_M + \langle Du, v \rangle_M = u(x_{max})v(x_{max}) - u(x_{min})v(x_{min})
```

ver1
```math
\underbrace{ \underline{u}^T M D \underline{v} + \underline{u}^T D^T M \underline{v} }_{\substack{\textstyle \text{ } \\ \approx \int_{x_{\min}}^{x_{\max}} u \, (\partial_x v) + \int_{x_{\min}}^{x_{\max}} (\partial_x u)\, v }} 
= 
\underbrace{ \underline{u}^T \underline{t}_R \underline{t}_R^T \underline{v} - \underline{u}^T \underline{t}_L \underline{t}_L^T \underline{v} }_{\substack{\textstyle \text{ } \\ \approx u(x_{\max}) v(x_{\max}) - u(x_{\min}) v(x_{\min}) }}

```

ver2



which is exactly the discrete analog of integration by parts!

### Example: Central Difference Operator

A simple example is the classical second-order central difference operator on a uniform grid:

```math
D = \frac{1}{\Delta x} \begin{pmatrix}
-1 & 1 & & & \\
-1/2 & 0 & 1/2 & & \\
& -1/2 & 0 & 1/2 & \\
& & \ddots & \ddots & \ddots \\
& & & -1 & 1
\end{pmatrix}
```

with mass matrix:
```math
M = \Delta x \cdot \text{diag}(1/2, 1, 1, \ldots, 1, 1/2)
```

You can verify that this satisfies the SBP property and provides second-order accuracy in the interior with first-order accuracy at the boundaries.

## 3. Types of SBP Operators

### Central vs. Upwind SBP Operators

The choice between central and upwind SBP operators is crucial and depends on the nature of your problem.

#### Central SBP Operators
Central operators are symmetric and provide the "natural" discretization:
- **Best for**: Advection-dominated problems, hyperbolic systems
- **Properties**: Conservative, minimal dissipation
- **Stencil**: Symmetric around each point
- **Drawback**: Can suffer from spurious oscillations in certain contexts

#### Upwind SBP Operators
Upwind operators come in pairs ``D_+`` and ``D_-`` with additional dissipation:
- **Key Property**: ``M(D_+ - D_-)`` is negative semidefinite (provides controlled dissipation)
- **Best for**: Parabolic/elliptic problems, diffusion-dominated systems
- **Major Advantage**: **Higher resolution** - better at resolving sharp features
- **Construction**: ``D = (D_+ + D_-)/2`` gives the central operator

#### Why Upwind Operators Have Higher Resolution

The resolution advantage comes from the stencil structure. Consider discretizing a second derivative:

**Central approach** (like BR1 method):
```
(D²u)ᵢ ≈ (uᵢ₊₂ - 2uᵢ + uᵢ₋₂)/(4Δx²)  [Wide stencil!]
```
This wide stencil has **grid oscillations** of the form ``(-1,+1,-1,+1,...)`` in its null space.

**Upwind approach** (like LDG method):
```
(D²u)ᵢ ≈ (uᵢ₊₁ - 2uᵢ + uᵢ₋₁)/Δx²     [Narrow stencil!]
```
The narrow stencil naturally damps high-frequency oscillations, leading to **better resolution**.

### Different SBP Implementations

#### Finite Difference (FD) SBP
- **Grid**: Uniform or non-uniform node distributions
- **Mass matrix**: Usually diagonal (especially on uniform grids)
- **Accuracy**: Limited by stencil width
- **Best for**: Simple geometries, high efficiency

#### Discontinuous Galerkin (DG) SBP
- **Grid**: Elements with high-order polynomials within each element
- **Mass matrix**: Block diagonal (one block per element)
- **Accuracy**: Very high order possible
- **Best for**: Complex geometries, hp-adaptivity
- **Special feature**: Natural upwind operators through numerical fluxes

#### Continuous Galerkin (CG) SBP
- **Grid**: Continuous finite element basis
- **Mass matrix**: Sparse but not diagonal
- **Accuracy**: High order with fewer degrees of freedom than DG
- **Best for**: Smooth solutions, memory efficiency

#### Fourier/Spectral SBP
- **Grid**: Fourier modes or other orthogonal polynomials
- **Mass matrix**: Often identity or diagonal in spectral space
- **Accuracy**: Exponential convergence for smooth periodic problems
- **Best for**: Problems with high regularity and periodic boundaries

### Choosing the Right Operator Type

| Problem Type | Recommended Operator | Reason |
|--------------|---------------------|--------|
| Wave propagation | Central FD/DG | Minimal dissipation preserves waves |
| Diffusion/Heat equation | Upwind DG | Better resolution of boundary layers |
| Shock-dominated flows | Upwind DG with limiters | Handles discontinuities |
| Periodic smooth solutions | Fourier | Spectral accuracy |
| Complex geometries | DG | Geometric flexibility |

### Resolution Comparison Example

Here's a simple comparison showing the resolution difference:

```julia
# Central operators - may have oscillations
solver_central = Solver(mesh, 4)  # Uses central differences

# Upwind operators - higher resolution
using SummationByPartsOperators: couple_discontinuously, PeriodicUpwindOperators
D_legendre = legendre_derivative_operator(-1.0, 1.0, 4)
uniform_mesh = UniformPeriodicMesh1D(xmin, xmax, div(N, 4))
minus = couple_discontinuously(D_legendre, uniform_mesh, Val(:minus))
plus = couple_discontinuously(D_legendre, uniform_mesh, Val(:plus))
central = couple_discontinuously(D_legendre, uniform_mesh)

D1 = PeriodicUpwindOperators(minus, central, plus)
D2 = sparse(plus) * sparse(minus)  # Creates narrow-stencil second derivative
solver_upwind = Solver(D1, D2)
```

The upwind solver will typically provide sharper resolution of features like wave fronts and boundary layers, at the cost of slightly more computational work per time step.



## [Customize solver](@id customize_solver)

In the semidiscretization created above, we used the default SBP operators, which are periodic finite difference operators. Using different SBP operators for the
semidiscretization can be done leveraging [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl/), which needs to be imported first:

```
using SummationByPartsOperators: legendre_derivative_operator, UniformPeriodicMesh1D, couple_discontinuously, PeriodicUpwindOperators
```

As an example, let us create a semidiscretization based on discontinuous Galerkin (DG) upwind operators. A semidiscretization implemented in DispersiveShallowWater.jl
needs one first-derivative and one second-derivative SBP operator. To build the first-derivative operator, we first create a `LegendreDerivativeOperator` with polynomial
degree 3 on a reference element `[-1.0, 1.0]` and a `UniformPeriodicMesh1D` for the coupling.

```
mesh = Mesh1D(coordinates_min, coordinates_max, N)
accuracy_order = 4
D_legendre = legendre_derivative_operator(-1.0, 1.0, accuracy_order)
uniform_mesh = UniformPeriodicMesh1D(mesh.xmin, mesh.xmax, div(mesh.N, accuracy_order))
```

Upwind DG operators in negative, central and positive operators can be obtained by `couple_discontinuously`

```
central = couple_discontinuously(D_legendre, uniform_mesh)
minus = couple_discontinuously(D_legendre, uniform_mesh, Val(:minus))
plus = couple_discontinuously(D_legendre, uniform_mesh, Val(:plus))
D1 = PeriodicUpwindOperators(minus, central, plus)
```

In order to still have an entropy-conserving semidiscretization the second-derivative SBP operator needs to be

```
using SparseArrays: sparse
D2 = sparse(plus) * sparse(minus)
```

The [`Solver`](@ref) object can now be created by passing the two SBP operators to the constructor, which, in turn, can be used to construct a `Semidiscretization`:

```
solver = Solver(D1, D2)
semi = Semidiscretization(mesh, equations, initial_condition, solver, boundary_conditions = boundary_conditions)
```

As before, we can run the simulation by

```
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 velocity, entropy),
                                     io = devnull)
relaxation_callback = RelaxationCallback(invariant = entropy)
callbacks = CallbackSet(relaxation_callback, analysis_callback)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
anim = @animate for step in 1:length(sol.u)
    plot(semi => sol, plot_initial = true, conversion = waterheight_total, step = step, xlim = (-50, 20), ylims = (-0.8, 0.1))
end
gif(anim, "shoaling_solution_dg.gif", fps = 25)
nothing # hide
```


shoaling solution DG_shoaling_solution_dg.gif


For more details see also the [documentation of SummationByPartsOperators.jl](https://ranocha.de/SummationByPartsOperators.jl/stable/)