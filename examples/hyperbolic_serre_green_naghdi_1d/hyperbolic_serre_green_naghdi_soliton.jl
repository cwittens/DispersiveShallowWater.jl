using OrdinaryDiffEqLowStorageRK
using DispersiveShallowWater

###############################################################################
# Semidiscretization of the hyperbolic Serre-Green-Naghdi equations

equations = HyperbolicSerreGreenNaghdiEquations1D(bathymetry_type = bathymetry_mild_slope,
                                                  lambda = 12000.0,
                                                  gravity = 9.81)

initial_condition = initial_condition_soliton
boundary_conditions = boundary_condition_reflecting
# boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -50.0
coordinates_max = 50.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with periodic SBP operators of accuracy order 4
accuracy_order = 2
solver = Solver(mesh, accuracy_order)

accuracy_order = 2
D1 = derivative_operator(MattssonNordström2004();
                         derivative_order = 1, accuracy_order,
                         xmin = xmin(mesh), xmax = xmax(mesh), N = N)

solver = Solver(D1)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation
tspan = (0.0, (xmax(mesh) - xmin(mesh)) / sqrt(1.2 * gravity(equations))) # one period
tspan = (0.0, 20.00)
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total, energy_total_modified,
                                                                 entropy_modified))
callbacks = CallbackSet(analysis_callback, summary_callback)

saveat = range(tspan..., length = 100)
# optimized time integration methods like this one are much more efficient
# for stiff problems (λ big) than standard methods like Tsit5()
alg = RDPK3SpFSAL35()
hyper_sol = solve(ode, alg, abstol = 1e-8, reltol = 1e-8,
            save_everystep = false, callback = callbacks, saveat = saveat);
plot(semi => hyper_sol, legend =false)

anim = @animate for step in 1:length(hyper_sol.u)
    plot(semi => hyper_sol, plot_initial = true, legend =false, step = step,)
end
gif(anim, fps = 20)

print("\a")

integrals(analysis_callback)

integrals(analysis_callback).energy_total_modified |> extrema |> x -> last(x) - first(x)

N = 10
tr, tl = zeros(N), zeros(N)
tr[end] = 1
tl[1] = 1

RBC_M = sparse(tr * tr' - tl * tl')

tmp = get_tmp(semi.cache.BC_Ref_Matrix, ode.u0.x[1])
tmp[end,end]

