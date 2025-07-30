# [Equations](@id equations)

DispersiveShallowWater.jl provides six different dispersive shallow water equation systems for modeling water waves. Each equation system offers different levels of physical accuracy, computational complexity, and supports various boundary conditions and bathymetry types.

## [Supported Models and Features](@id eq_overview)

The following table provides an overview of all available equation systems and their supported features:

| Equation | Variables | [Periodic boundary conditions](@ref boundary_condition_periodic) | [Reflecting boundary conditions](@ref boundary_condition_reflecting) | [Flat Bathymetry](@ref bathymetry_flat) | [Mild-slope Bathymetry](@ref bathymetry_mild_slope) | [Variable Bathymetry](@ref bathymetry_variable) | Relaxation | Source Terms |
|----------|:---------:|:-----------:|:-------------:|:----:|:-----------:|:--------:|:----------:|:-------:|
| [`KdV`](@ref KdVEquation1D) | ``(\eta)`` | ✅ | ❌ | ✅ | ❌ | ❌ | ✅ | ✅ |
| [`BBM`](@ref BBMEquation1D) | ``(\eta)`` | ✅ | ❌ | ✅ | ❌ | ❌ | ✅ | ✅ |
| [`BBM-BBM`](@ref BBMBBMEquations1D) | ``(\eta, v, D)`` | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ |
| [`Svärd-Kalisch`](@ref SvaerdKalischEquations1D) | ``(\eta, v, D)`` | ✅ | ✅ᵃ | ❌ | ❌ | ✅ | ✅ | ✅ |
| [`Serre-Green-Naghdi`](@ref SerreGreenNaghdiEquations1D) | ``(\eta, v, D)`` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| [`Hyperbolic SGN`](@ref HyperbolicSerreGreenNaghdiEquations1D) |``(\eta, v, D, w, H)`` | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ |

*ᵃReflecting boundary conditions for Svärd-Kalisch equations require `alpha = gamma = 0`*

### Variable Descriptions

- ``\eta``: Total water height
- ``v``: Velocity in horizontal direction  
- ``D``: Still-water depth
- ``w``: Auxiliary variable in hyperbolic approximation (``\approx -h v_x``)
- ``H``: Auxiliary variable in hyperbolic approximation (``\approx h``)

## Equation Categories

The equations can be categorized based on their complexity and number of variables:

### Single Variable Equations

These equations describe wave propagation using a single primary variable (the total water height η) and are computationally efficient but with limited physical accuracy:

- **KdV (Korteweg-de Vries)**: The classic integrable equation for weakly nonlinear, long waves
- **BBM (Benjamin-Bona-Mahony)**: An alternative to KdV with better stability properties

### Multi-Variable Systems

These systems use multiple variables for more accurate wave modeling, capturing both water height and velocity dynamics:

- **BBM-BBM**: A coupled system extending BBM to include velocity evolution
- **Svärd-Kalisch**: A recent model with tunable parameters for optimal dispersion properties
- **Serre-Green-Naghdi**: A comprehensive system with excellent dispersion characteristics
- **Hyperbolic SGN**: A hyperbolic approximation of Serre-Green-Naghdi with additional auxiliary variables

## Detailed Documentation

Each equation system below includes its complete mathematical formulation, physical background, implementation details, conservation properties, and relevant literature references.

## Benjamin-Bona-Mahony (BBM)

```@docs
DispersiveShallowWater.BBMEquation1D
```

## Korteweg–De Vries (KdV)

```@docs
DispersiveShallowWater.KdVEquation1D
```

## BBM-BBM

```@docs
DispersiveShallowWater.BBMBBMEquations1D
```

## Svärd-Kalisch

```@docs
DispersiveShallowWater.SvaerdKalischEquations1D
DispersiveShallowWater.SvärdKalischEquations1D
```

## Serre-Green-Naghdi

```@docs
DispersiveShallowWater.SerreGreenNaghdiEquations1D
```

## Hyperbolic Serre-Green-Naghdi

```@docs
DispersiveShallowWater.HyperbolicSerreGreenNaghdiEquations1D
```
