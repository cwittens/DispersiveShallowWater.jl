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
