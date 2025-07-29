# Miscellaneous


## Callbacks 

### Summary Callback

HERE PUT A SHORT DISCRITOPN ABOUT THE SummaryCallback()


Eg that it will return something like
```
──────────────────────────────────────────────────────────────────────────────────────
           DispersiveSWE                     Time                    Allocations
                                    ───────────────────────   ────────────────────────
         Tot / % measured:              7.04ms /  54.6%            420KiB /   0.3%

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
rhs!                         3.56k   3.84ms  100.0%  1.08μs   1.25KiB  100.0%    0.36B
  hyperbolic                 3.56k   1.34ms   34.8%   377ns     0.00B    0.0%    0.00B
  third-order derivatives    3.56k   1.32ms   34.2%   370ns     0.00B    0.0%    0.00B
  ~rhs!~                     3.56k   1.11ms   28.8%   312ns   1.25KiB  100.0%    0.36B
  source terms               3.56k   81.2μs    2.1%  22.8ns     0.00B    0.0%    0.00B
──────────────────────────────────────────────────────────────────────────────────────
```

in the end. So show an overview of the time it took and the number of alliocations etc

use it by doing:
```
summary_callback = SummaryCallback()
sol = solve(ode, Tsit5(), callback = summary_callback)
```


### Analysis Callback

Additionally, we can analyze the numerical solution using an [`AnalysisCallback`](@ref).
The analysis includes computing the ``L^2`` error and ``L^\infty`` error of the different solution's variables compared to the initial condition (or, if available,
at the same time analytical solution). Additional errors can be passed by the keyword argument `extra_analysis_errors`. Additional integral quantities that should
be analyzed can be passed by keyword argument `extra_analysis_integrals`. In this example we pass the `conservation_error`, which computes the temporal change of
the total amount (i.e. integral) of the different variables over time. In addition, the integrals of the total water height ``\eta`` [`waterheight_total`](@ref),
the [`velocity`](@ref) and the [`entropy`](@ref) are computed and saved for each time step. The total water height and the total velocity are linear invariants of
the BBM-BBM equations, i.e. they do not change over time. The total entropy
HERE MAKE CLEAR THAT THIS VARIES FOR THE DIFFERENT MODELS
```math
\mathcal E(t; \eta, v) = \frac{1}{2}\int_\Omega g\eta^2 + (\eta + D)v^2\textrm{d}x
```

is a nonlinear invariant and should be constant over time as well. During the simulation, the `AnalysisCallback` will print the results to the terminal.

THIS WAS ORIGINIALLY PART OF THE OVERVIEW. NOW THAT IT IS STAND ALONE MAKE SURE ALL THE RELEVANT CODE FROM BEFORE TO DEFINE THE SETUP WILL ALSO BE HERE AND EXECUTED WHEN DOING DOCUMENTER.JL BY DOING ```@example callback ...
USE THE SGN EQUATIONS HERE.
```
tspan = (0.0, 25.0)
ode = semidiscretize(semi, tspan)
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 velocity, entropy),
                                     io = devnull)
summary_callback = SummaryCallback() # also include the summary callback here
callbacks = CallbackSet(analysis_callback, summary_callback)

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

![analysis callback](analysis_callback.png)

You can see that the linear invariants ``\int_\Omega\eta\textrm{d}x`` and ``\int_\Omega v\textrm{d}x`` are indeed conserved exactly. The entropy, however, starts
growing at around ``t = 17``  and rises up to approximately `5e-5`. This is because of the fact that, during the time integration, a nonlinear invariant is not
necessarily conserved, even if the semidiscretization conserves the quantity exactly. How to obtain a fully-discrete structure-preserving numerical scheme is explained
in the following section.


### Relaxation Callback

## Use entropy-conserving time integration

To obtain entropy-conserving time-stepping schemes DispersiveShallowWater.jl uses the relaxation method introduced in [^Ketcheson2019] and further developed in
[^RanochaSayyariDalcinParsaniKetcheson2020]. The relaxation method is implemented as a [`RelaxationCallback`](@ref), which takes a function representing the conserved
quantity as the keyword argument `invariant`. Therefore, we can run the same example as above, but using relaxation on the entropy by simply adding another callback
to the `CallbackSet`:

USE THE SETUP FROM BEFORE TO DO THE RELAXATION CALLBACK NOW

```@example callback
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 velocity, entropy),
                                     io = devnull)
relaxation_callback = RelaxationCallback(invariant = entropy)
# Always put relaxation_callback before analysis_callback to guarantee conservation of the invariant
callbacks = CallbackSet(relaxation_callback, analysis_callback)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
```

When you use both, an `AnalysisCallback` and a `RelaxationCallback`, note that the `relaxation_callback` needs to come first inside the `CallbackSet` as it needs to be
invoked prior to the `analysis_callback`, such that the `analysis_callback` analyzes the solution with the already updated values.

Plotting the `analysis_callback` again, we can see that now also the `entropy` is conserved up to machine precision.

```@example callback
plot(analysis_callback, ylims = (-5e-16, 5e-16))
savefig("analysis_callback_relaxation.png") # hide
nothing # hide
```

![analysis callback relaxation](analysis_callback_relaxation.png)






[^RanochaSayyariDalcinParsaniKetcheson2020]:
    Ranocha, Sayyari, Dalcin, Parsani, Ketcheson (2020):
    Relaxation Runge–Kutta Methods: Fully-Discrete Explicit Entropy-Stable Schemes for the Compressible Euler and Navier–Stokes Equations
    [DOI: 10.1137/19M1263480](https://doi.org/10.1137/19M1263480)


[^Ketcheson2019]:
    Ketcheson (2019):
    Relaxation Runge-Kutta Methods: Conservation and stability for Inner-Product Norms.
    [DOI: 10.1137/19M1263662](https://doi.org/10.1137/19M1263662)




## conversion functions

You can also provide a `conversion` function that converts the solution. A conversion function should take the values
of the primitive variables `q` at one node, and the `equations` as input and should return an `SVector` of any length as output. For a user defined conversion function,
there should also exist a function `varnames(conversion, equations)` that returns a `Tuple` of the variable names used for labelling. The conversion function can, e.g.,
be [`prim2cons`](@ref) or [`waterheight_total`](@ref) if one only wants to plot the total water height. The resulting plot will have one subplot for each of the returned
variables of the conversion variable. By default, the conversion function is just [`prim2phys`](@ref), which computes the physical variables
from the primitive ones. For most equations this is the identity, but for hyperbolic approximations like [`HyperbolicSerreGreenNaghdiEquations1D`](@ref), it returns a reduced set of variables (excluding auxiliary variables like `w` and `H`).
