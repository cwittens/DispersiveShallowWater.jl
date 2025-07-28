# Miscellaneous

## Equation Support Overview

The following table provides an overview of the supported features for each equation type in DispersiveShallowWater.jl:

TODO: check for correctness!!

| Equation | Variables | Periodic boundary conditions | Reflecting boundary conditions | Flat Bathymetry | Mild-slope Bathymetry | Variable Bathymetry | Relaxation | Source Terms |
|----------|:---------:|:-----------:|:-------------:|:---------------:|:---------------------:|:-------------------:|:----------:|:------------:|
| [`BBMEquation1D`](@ref) | 1 (`η`) | ✅ | ❌ | ✅ | ❌ | ❌ | ✅ | ✅ |
| [`BBMBBMEquations1D`](@ref) | 3 (`η`, `v`, `D`) | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ |
| [`KdVEquation1D`](@ref) | 1 (`η`) | ✅ | ❌ | ✅ | ❌ | ❌ | ✅ | ✅ |
| [`SvaerdKalischEquations1D`](@ref) | 3 (`η`, `v`, `D`) | ✅ | ✅* | ❌ | ❌ | ✅ | ✅ | ✅ |
| [`SerreGreenNaghdiEquations1D`](@ref) | 3 (`η`, `v`, `D`) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| [`HyperbolicSerreGreenNaghdiEquations1D`](@ref) | 5 (`η`, `v`, `D`, `w`, `H`) | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ |

*\* Reflecting boundary conditions for Svärd-Kalisch equations require `alpha = gamma = 0`*

## Legend

- **Periodic BC**: Periodic boundary conditions
- **Reflecting BC**: Reflecting boundary conditions  
- **Variables**: Number and type of primitive variables
- **Flat Bathymetry**: Constant bathymetry support
- **Mild-slope Bathymetry**: Variable bathymetry with mild-slope approximation support
- **Variable Bathymetry**: General variable bathymetry support (no approximation)
- **Relaxation**: Support for relaxation methods to preserve invariants
- **Source Terms**: Support for source terms in the equations

## Bathymetry Types

The package supports three types of bathymetry:

- [`bathymetry_flat`](@ref): Flat bathymetry (typically `b = 0` everywhere)
- [`bathymetry_mild_slope`](@ref): Variable bathymetry with mild-slope approximation (some terms like `b_x²` are neglected)
- [`bathymetry_variable`](@ref): General variable bathymetry without approximation

## Conservation Properties

All equations in DispersiveShallowWater.jl are designed to preserve certain quantities:

### Linear Invariants
- **Mass conservation**: All equations conserve the total water mass
- **Momentum conservation**: Most multi-variable equations conserve momentum (for appropriate boundary conditions)

### Nonlinear Invariants  
- **Energy/Entropy**: All equations can preserve a form of energy or entropy using relaxation methods
- **Modified energy**: Some equations preserve modified energy that includes dispersive terms

## Boundary Condition Notes

### Periodic Boundary Conditions
- Supported by all equations
- Require periodic SBP operators from SummationByPartsOperators.jl
- Natural choice for problems with wave propagation in periodic domains

### Reflecting Boundary Conditions  
- Require non-periodic SBP operators
- Some equations have restrictions (e.g., Svárd-Kalisch requires specific parameter values)
- Useful for problems with solid walls or impermeable boundaries

## Variable Descriptions

- `η`: Total water height (waterheight_total)
- `v`: Velocity in horizontal direction  
- `D`: Still-water depth
- `w`: Auxiliary variable in hyperbolic approximation (≈ -h v_x)
- `H`: Auxiliary variable in hyperbolic approximation (≈ h)


