# Reactive Euler Solver

## Introduction

This project is a finite volume solver framework. It is tailored to do simulations for the euler equations with a kinetic source term and with respect to different multi-species ideal gas equation of states which depend on a specified reaction mechanism. 

The Reaction Euler Solver is designed to provide multiple exchangable grid backends which control how to parallelise or distribute the computation data across possibly many computation nodes for computer clusters.

## Project Overview

The solver is divided into 4 abstraction layers:

- Core: The core layer is lowest level of abstraction where fundamental vocabulary types are defined such as `span<T>` or `simd<T, Abi>`.
- Grid: The grid layer defines variable notation and the different grid backends can be found here.
- Solver: The solver layer contains source code for the numerical schemes which solve the differential equations for patch types defined by the grid layer.
- Driver: The driver layer provides bundled interfaces for each equation of state and dimension to solve a specified problem.
- App: The app layer finally defines executables to solve concrete initial value problems.

## Build Instructions

TODO

## Testing

TODO

## Contribution

TODO

## References

TODO