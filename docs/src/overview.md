# Overview

## [Models in DispersiveShallowWater.jl](@id eq_overview)

DispersiveShallowWater.jl currently provides six different dispersive shallow water equation systems that can be used for modeling water waves. Each equation system offers different levels of physical accuracy, computational complexity, and supports various boundary conditions and bathymetry types. To learn more about the mathematical definition of each model, go to the chapter [Equations](@id all_eq). The following table provides an overview of the supported features for each equation type:

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

# [Equations](@id all_eq)

DispersiveShallowWater.jl provides several dispersive shallow water equation systems, each with different levels of physical accuracy and computational complexity.

## Overview

The package currently implements the following equation systems:

- [`KdVEquation1D`](@ref): Korteweg-de Vries equation
- [`BBMEquation1D`](@ref): Benjamin-Bona-Mahony equation  
- [`BBMBBMEquations1D`](@ref): BBM-BBM system
- [`SvaerdKalischEquations1D`](@ref): Svärd-Kalisch equations
- [`SerreGreenNaghdiEquations1D`](@ref): Serre-Green-Naghdi equations
- [`HyperbolicSerreGreenNaghdiEquations1D`](@ref): Hyperbolic approximation of Serre-Green-Naghdi equations

## Single Variable Equations

```@docs
DispersiveShallowWater.KdVEquation1D
DispersiveShallowWater.BBMEquation1D
```

## Multi-Variable Systems

```@docs
DispersiveShallowWater.BBMBBMEquations1D
DispersiveShallowWater.SvaerdKalischEquations1D
DispersiveShallowWater.SerreGreenNaghdiEquations1D
DispersiveShallowWater.HyperbolicSerreGreenNaghdiEquations1D
```
