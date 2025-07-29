# [Solvers](@id solvers)

This chapter covers different solvers and how to use them in DispersivShallowWater.jl.

To learn more about the analytical and mathematical background, go to the chapter about [Summation by Parts Operators](@ref sbp_operates).

## Basic Summation by Parts Operator

TALK ABOUT IT USING [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl/) AS A SHORT INTRO

As we have already seen in the [basic example](@ref basic_example), the easiest way to create a solver is to pass both a mesh and the wanted accuracy order to the [`Solver`](@ref) function, which creates a first,- second- and third periodic summation by parts operators of given accuracy order on the given mesh.

```julia
accuracy_order = 4
solver = Solver(mesh, accuracy_order)
```

This however only work when working with periodic boundary condition. If a solver containing periodic SBP operators is passed in to `Semidiscretization` with non-peridoic boundary conditions, it will throw an error.


## Creating your own solver

*how to create your own solver by doing solver = Solver(D1, D2) with different examples for
- upwind operators for the KdV eq


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