"""
    AbstractApproximation

Abstract supertype for phase-space approximation methods.

Approximation objects describe how the quantum problem is approximated.
They are intentionally separate from the physical model, initial state,
environment, solver configuration, and random-number generator.
"""
abstract type AbstractApproximation end
