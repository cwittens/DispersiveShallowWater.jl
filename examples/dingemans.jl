using OrdinaryDiffEqTsit5
using DispersiveShallowWater
using Plots
###############################################################################

bbm = BBMBBMEquations1D(bathymetry_type = bathymetry_variable,
                        gravity = 9.81, eta0 = 0.0)

sk = SvaerdKalischEquations1D(gravity = 9.81, eta0 = 0.8, alpha = 0.0,
                              beta = 0.27946992481203003, gamma = 0.0521077694235589)

sgn = SerreGreenNaghdiEquations1D(bathymetry_type = bathymetry_variable,
                                  gravity = 9.81)

hysgn = HyperbolicSerreGreenNaghdiEquations1D(bathymetry_type = bathymetry_mild_slope,
                                              lambda = 500.0,
                                              gravity = 9.81)

initial_condition = initial_condition_dingemans
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -138.0
coordinates_max = 46.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with periodic SBP operators of accuracy order 4
accuracy_order = 4
solver = Solver(mesh, accuracy_order)

tspan = (0.0, 70.0)
saveat = range(tspan..., length = 500)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi_bbm = Semidiscretization(mesh, bbm, initial_condition, solver,
                              boundary_conditions = boundary_conditions)

semi_sk = Semidiscretization(mesh, sk, initial_condition, solver,
                             boundary_conditions = boundary_conditions)

semi_sgn = Semidiscretization(mesh, sgn, initial_condition, solver,
                              boundary_conditions = boundary_conditions)

semi_hysgn = Semidiscretization(mesh, hysgn, initial_condition, solver,
                                boundary_conditions = boundary_conditions)

ode_bbm = semidiscretize(semi_bbm, tspan)
ode_sk = semidiscretize(semi_sk, tspan)
ode_sgn = semidiscretize(semi_sgn, tspan)
ode_hysgn = semidiscretize(semi_hysgn, tspan)

sol_bbm = solve(ode_bbm, Tsit5(), abstol = 1e-7, reltol = 1e-7,
                save_everystep = false, saveat = saveat)
sol_sk = solve(ode_sk, Tsit5(), abstol = 1e-7, reltol = 1e-7,
               save_everystep = false, saveat = saveat)
sol_sgn = solve(ode_sgn, Tsit5(), abstol = 1e-7, reltol = 1e-7,
                save_everystep = false, saveat = saveat)
sol_hysgn = solve(ode_hysgn, Tsit5(), abstol = 1e-7, reltol = 1e-7,
                  save_everystep = false, saveat = saveat)

# BBM-BBM equations need to be translated vertically
# in order to do this, we define a new conversion function
shifted_waterheight(q, equations) = waterheight_total(q, equations) + 0.8
DispersiveShallowWater.varnames(shifted_waterheight, equations) = ("η",)

times = [14.0, 28.0, 42.0, 70.0]
PLOTS = [plot(), plot(), plot(), plot()]

for i in 1:4
    p = PLOTS[i]
    step_idx = argmin(abs.(saveat .- times[i]))

    plot!(p, semi_bbm => sol_bbm, label = "BBM", conversion = shifted_waterheight,
          suptitle = "Dingemans", title = "", plot_bathymetry = false, subplot = 1,
          step = step_idx)
    plot!(p, semi_sk => sol_sk, label = "Svärd-Kalisch", conversion = waterheight_total,
          plot_bathymetry = false, step = step_idx)
    plot!(p, semi_sgn => sol_sgn, label = "Serre-Green-Naghdi",
          conversion = waterheight_total, plot_bathymetry = false, step = step_idx)
    plot!(p, semi_hysgn => sol_hysgn, label = "Hyperbolic Serre-Green-Naghdi",
          conversion = waterheight_total, plot_bathymetry = true,
          suptitle = "Dingemans at t = $(times[i])", title = "", legend = false,
          step = step_idx, ylims = (-0.025239684626944536, 0.8665625055250958))
end

# extra plot which is just a legend
p_legend = plot([], [], label = "BBM", legend=:top, grid = false, showaxis = false, size = (100, 100))
plot!(p_legend, [], [],label = "Svärd-Kalisch")
plot!(p_legend, [], [],label = "Serre-Green-Naghdi")
plot!(p_legend, [], [],label = "Hyperbolic Serre-Green-Naghdi")
plot!(p_legend, [], [],label = "Bathymetry", color = :black,)
push!(PLOTS, p_legend)

plot(PLOTS..., size = (900, 900), layout = @layout([a b; c d; e]), 
     xlabel = "x", ylabel = "η")



# Define simulation parameters
times = [14.0, 28.0, 42.0, 70.0]
y_limits = (-0.03, 0.87)

# Model configurations: (semidiscretization, solution, label, conversion_function)
models = [
    (semi_bbm, sol_bbm, "BBM", shifted_waterheight),
    (semi_sk, sol_sk, "Svärd-Kalisch", waterheight_total),
    (semi_sgn, sol_sgn, "Serre-Green-Naghdi", waterheight_total),
    (semi_hysgn, sol_hysgn, "Hyperbolic Serre-Green-Naghdi", waterheight_total)
]

# Create snapshot plots for each time
snapshot_plots = []
for time_val in times
    step_idx = argmin(abs.(saveat .- time_val))
    p = plot(title = "t = $time_val", ylims = y_limits)
    
    for (i, (semi, sol, label, conversion)) in enumerate(models)

        plot!(p, semi => sol, 
              label = label,
              conversion = conversion,
              plot_bathymetry = true,
              step = step_idx,
              legend = false, suptitle = "Dingemans at t = $(time_val)", title = "")
    end
    
    push!(snapshot_plots, p)
end

# Create legend plot ()
legend_plot = plot(legend=:top, framestyle = :none, legendfontsize = 11)

for (_, _, label, _) in models
    plot!(legend_plot, [], [], label = label)
end
plot!(legend_plot, [], [], label = "Bathymetry", color = :black)

# Combine all plots
all_plots = [snapshot_plots..., legend_plot]
plot(all_plots..., 
     size = (900, 800), 
     layout = @layout([a b; c d; e]),
)