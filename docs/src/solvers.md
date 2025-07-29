# [Solvers](@id solvers)

This chapter covers different solvers and how to use them in DispersivShallowWater.jl.

To learn more about the analytical and mathematical background, go to the chapter about [Summation by Parts Operators](@ref sbp_operates).

## Basic Summation by Parts Operator

TALK ABOUT THAT THIS REPO IS USING [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl/) AS A SHORT INTRO

As we have already seen in the [basic example](@ref basic_example), the easiest way to create a solver is to pass both a mesh and the wanted accuracy order to the [`Solver`](@ref) function, which creates a first,- second- and third periodic summation by parts operators of given accuracy order on the given mesh.

```julia
accuracy_order = 4
solver = Solver(mesh, accuracy_order)
```

This however only work when working with periodic boundary condition. If a solver containing periodic SBP operators is passed in to `Semidiscretization` with non-peridoic boundary conditions, it will throw an error.


## [Customize solver](@id customize_solver)

explain here how to create your own solver by doing solver = Solver(D1, D2) with different examples for

- normal SPB operators when dealing with reflecting boundary conditions (as in DispersiveShallowWater.jl/examples/svaerd_kalisch_1d/svaerd_kalisch_1d_basic_reflecting.jl)
- upwind operators for the SGN eq, where they even have their own special semi discretization specialised on Upwind SBP Operators (as in DispersiveShallowWater.jl/examples/serre_green_naghdi_1d/serre_green_naghdi_soliton_upwind.jl for periodic BC and DispersiveShallowWater.jl/examples/serre_green_naghdi_1d/serre_green_naghdi_manufactured_reflecting_upwind.jl for reflecting BC)
- KdV eq. where also the thired direivate is just, so one can pass solver = Solver(D1, nothing, D3) (as in DispersiveShallowWater.jl/examples/kdv_1d/kdv_1d_fourier.jl)
- example cases  where the second deriv operator is a sparse matrix (here explain that D1 always needs to be a SBP Operator but D2 and D3 can be of supertype AbstractMatrix) but also highlit that SBPOperator*vector is about a order of magnitide faster then sparsematrix*vector (as can be seen in https://ranocha.de/SummationByPartsOperators.jl/stable/benchmarks/) (as in DispersiveShallowWater.jl/examples/bbm_bbm_1d/bbm_bbm_1d_upwind_relaxation.jl)
- find example from example folder and do then for legendre_derivative_operator and couple_discontinuously (as in DispersiveShallowWater.jl/examples/bbm_bbm_1d/bbm_bbm_1d_dg.jl)
- fourier SBP operators as in (DispersiveShallowWater.jl/examples/hyperbolic_serre_green_naghdi_1d/hyperbolic_serre_green_naghdi_soliton_relaxation.jl)

where each example should be its own subsection.
Of course dont just copy the whole example, only show what is relevant.

At the end say that this solver object can then be but into `Semidiscretization`.



