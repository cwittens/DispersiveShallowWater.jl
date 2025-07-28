# Miscellaneous

## Equation Support Overview

The following table provides an overview of the supported features for each equation type in DispersiveShallowWater.jl:

TODO: check for correctness!!

| Equation | Variables | Periodic boundary conditions | Reflecting boundary conditions | Flat Bathymetry | Mild-slope Bathymetry* | Variable Bathymetry | Relaxation | Source Terms |
|----------|:---------:|:-----------:|:-------------:|:---------------:|:---------------------:|:-------------------:|:----------:|:------------:|
| [`BBMEquation1D`](@ref) | 1 (`η`) | ✅ | ❌ | ✅ | ❌ | ❌ | ✅ | ✅ |
| [`BBMBBMEquations1D`](@ref) | 3 (`η`, `v`, `D`) | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ |
| [`KdVEquation1D`](@ref) | 1 (`η`) | ✅ | ❌ | ✅ | ❌ | ❌ | ✅ | ✅ |
| [`SvaerdKalischEquations1D`](@ref) | 3 (`η`, `v`, `D`) | ✅ | ✅**| ❌ | ❌ | ✅ | ✅ | ✅ |
| [`SerreGreenNaghdiEquations1D`](@ref) | 3 (`η`, `v`, `D`) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| [`HyperbolicSerreGreenNaghdiEquations1D`](@ref) | 5 (`η`, `v`, `D`, `w`, `H`) | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ |

*\* Variable bathymetry with mild-slope approximation (some higher order terms are neglected)*

*\*\*Reflecting boundary conditions for Svärd-Kalisch equations require `alpha = gamma = 0`*


| Equation | Variables | Periodic boundary conditions | Reflecting boundary conditions | Flat Bathymetry | Mild-slope Bathymetry* | Variable Bathymetry | Relaxation | Source Terms |
|----------|:---------:|:-----------:|:-------------:|:----:|:-----------:|:--------:|:----------:|:-------:|
| [`BBM`](@ref BBMEquation1D) | 1 (`η`) | ✅ | ❌ | ✅ | ❌ | ❌ | ✅ | ✅ |
| [`BBM-BBM`](@ref BBMBBMEquations1D) | 3 (`η`, `v`, `D`) | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ |
| [`KdV`](@ref KdVEquation1D) | 1 (`η`) | ✅ | ❌ | ✅ | ❌ | ❌ | ✅ | ✅ |
| [`Svärd-Kalisch`](@ref SvaerdKalischEquations1D) | 3 (`η`, `v`, `D`) | ✅ | ✅** | ❌ | ❌ | ✅ | ✅ | ✅ |
| [`Serre-Green-Naghdi`](@ref SerreGreenNaghdiEquations1D) | 3 (`η`, `v`, `D`) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| [`Hyperbolic SGN`](@ref HyperbolicSerreGreenNaghdiEquations1D) | 5 (`η`, `v`, `D`, `w`, `H`) | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ |


## AnalysisCallback (maybe put it somewhere else but not in overview!)
And also summary_callback = SummaryCallback()?

Additionally, we can analyze the numerical solution using an [`AnalysisCallback`](@ref).
The analysis includes computing the ``L^2`` error and ``L^\infty`` error of the different solution's variables compared to the initial condition (or, if available,
at the same time analytical solution). Additional errors can be passed by the keyword argument `extra_analysis_errors`. Additional integral quantities that should
be analyzed can be passed by keyword argument `extra_analysis_integrals`. In this example we pass the `conservation_error`, which computes the temporal change of
the total amount (i.e. integral) of the different variables over time. In addition, the integrals of the total water height ``\eta`` [`waterheight_total`](@ref),
the [`velocity`](@ref) and the [`entropy`](@ref) are computed and saved for each time step. The total water height and the total velocity are linear invariants of
the BBM-BBM equations, i.e. they do not change over time. The total entropy

```math
\mathcal E(t; \eta, v) = \frac{1}{2}\int_\Omega g\eta^2 + (\eta + D)v^2\textrm{d}x
```

is a nonlinear invariant and should be constant over time as well. During the simulation, the `AnalysisCallback` will print the results to the terminal.


```
tspan = (0.0, 25.0)
ode = semidiscretize(semi, tspan)
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 velocity, entropy),
                                     io = devnull)
callbacks = CallbackSet(analysis_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
```


The errors and
integrals recorded by the `AnalysisCallback` can be obtained as `NamedTuple`s by [`errors(analysis_callback)`](@ref) and [`integrals(analysis_callback)`](@ref).

Often, it is interesting to have a look at how the quantities that are recorded by the `AnalysisCallback` evolve in time. To this end, you can `plot` the `AnalysisCallback` by

```
plot(analysis_callback)
savefig("analysis_callback.png") # hide
nothing # hide
```

This creates the following figure:

analysis callback_analysis_callback.png

You can see that the linear invariants ``\int_\Omega\eta\textrm{d}x`` and ``\int_\Omega v\textrm{d}x`` are indeed conserved exactly. The entropy, however, starts
growing at around ``t = 17``  and rises up to approximately `5e-5`. This is because of the fact that, during the time integration, a nonlinear invariant is not
necessarily conserved, even if the semidiscretization conserves the quantity exactly. How to obtain a fully-discrete structure-preserving numerical scheme is explained
in the following section.


## conversion functions

You can also provide a `conversion` function that converts the solution. A conversion function should take the values
of the primitive variables `q` at one node, and the `equations` as input and should return an `SVector` of any length as output. For a user defined conversion function,
there should also exist a function `varnames(conversion, equations)` that returns a `Tuple` of the variable names used for labelling. The conversion function can, e.g.,
be [`prim2cons`](@ref) or [`waterheight_total`](@ref) if one only wants to plot the total water height. The resulting plot will have one subplot for each of the returned
variables of the conversion variable. By default, the conversion function is just [`prim2phys`](@ref), which computes the physical variables
from the primitive ones. For most equations this is the identity, but for hyperbolic approximations like [`HyperbolicSerreGreenNaghdiEquations1D`](@ref), it returns a reduced set of variables (excluding auxiliary variables like `w` and `H`).
