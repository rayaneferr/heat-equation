"""Heat equation solvers (1D & 2D) using finite differences."""

from .heat import solve_heat_1d, solve_heat_2d, cfl_number_1d, cfl_number_2d

__all__ = [
    "solve_heat_1d",
    "solve_heat_2d",
    "cfl_number_1d",
    "cfl_number_2d",
]
