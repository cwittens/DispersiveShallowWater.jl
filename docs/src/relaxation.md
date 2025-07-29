## Use entropy-conserving time integration

## [Relaxation](@id relaxation_explained)

delete this file???



To obtain entropy-conserving time-stepping schemes DispersiveShallowWater.jl uses the relaxation method introduced in [^Ketcheson2019] and further developed in
[^RanochaSayyariDalcinParsaniKetcheson2020]. The relaxation method is implemented as a [`RelaxationCallback`](@ref), which takes a function representing the conserved
quantity as the keyword argument `invariant`. Therefore, we can run the same example as above, but using relaxation on the entropy by simply adding another callback
to the `CallbackSet`:

```
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

```
plot(analysis_callback, ylims = (-5e-16, 5e-16))
savefig("analysis_callback_relaxation.png") # hide
nothing # hide
```

analysis callback relaxation_analysis_callback_relaxation.png



[^RanochaSayyariDalcinParsaniKetcheson2020]:
    Ranocha, Sayyari, Dalcin, Parsani, Ketcheson (2020):
    Relaxation Runge–Kutta Methods: Fully-Discrete Explicit Entropy-Stable Schemes for the Compressible Euler and Navier–Stokes Equations
    [DOI: 10.1137/19M1263480](https://doi.org/10.1137/19M1263480)


[^Ketcheson2019]:
    Ketcheson (2019):
    Relaxation Runge-Kutta Methods: Conservation and stability for Inner-Product Norms.
    [DOI: 10.1137/19M1263662](https://doi.org/10.1137/19M1263662)
