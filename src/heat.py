"""
Résolution numérique de l'équation de la chaleur par différences finies.

Équation de la chaleur (diffusion thermique) :

    1D :  du/dt = alpha * d2u/dx2
    2D :  du/dt = alpha * (d2u/dx2 + d2u/dy2)

où `u` est la température, `alpha` la diffusivité thermique.

On utilise le schéma explicite FTCS (Forward-Time Central-Space) :
la dérivée temporelle est approchée à l'avant, les dérivées spatiales
par un Laplacien centré. Ce schéma est simple mais *conditionnellement
stable* : il faut respecter la condition CFL (voir `cfl_number_*`).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


# --------------------------------------------------------------------------- #
# Conditions de stabilité (CFL)
# --------------------------------------------------------------------------- #
def cfl_number_1d(alpha: float, dt: float, dx: float) -> float:
    """Nombre de Fourier r = alpha*dt/dx^2 ; le schéma 1D est stable si r <= 0.5."""
    return alpha * dt / dx**2


def cfl_number_2d(alpha: float, dt: float, dx: float, dy: float) -> float:
    """Critère de stabilité 2D : alpha*dt*(1/dx^2 + 1/dy^2) <= 0.5."""
    return alpha * dt * (1.0 / dx**2 + 1.0 / dy**2)


def _auto_dt(margin: float, *terms: float) -> float:
    """dt maximal stable (= 0.5 / sum(terms)) réduit d'une marge de sécurité."""
    return margin * 0.5 / sum(terms)


# --------------------------------------------------------------------------- #
# Solveur 1D
# --------------------------------------------------------------------------- #
def solve_heat_1d(
    *,
    length: float = 1.0,
    nx: int = 201,
    alpha: float = 1.0,
    t_max: float = 0.1,
    dt: float | None = None,
    initial: Array | None = None,
    bc: tuple[float, float] = (0.0, 0.0),
    store_every: int = 1,
    safety: float = 0.9,
) -> tuple[Array, Array, Array]:
    """
    Résout l'équation de la chaleur 1D sur [0, length] avec conditions de
    Dirichlet (températures fixées aux bords).

    Paramètres principaux :
        length     : longueur du domaine
        nx         : nombre de points de la grille
        alpha      : diffusivité thermique
        t_max      : temps final de simulation
        dt         : pas de temps (si None, choisi automatiquement et stable)
        initial    : profil initial u(x, 0) (sinon créneau chaud au centre)
        bc         : températures (gauche, droite) imposées aux bords
        store_every: on ne mémorise qu'une frame sur `store_every` (animation)

    Retour : (x, times, frames) où frames a la forme (n_frames, nx).
    """
    x = np.linspace(0.0, length, nx)
    dx = x[1] - x[0]

    # Pas de temps : automatique et stable si non fourni
    if dt is None:
        dt = _auto_dt(safety, alpha / dx**2)

    r = cfl_number_1d(alpha, dt, dx)
    if r > 0.5:
        raise ValueError(
            f"Schéma instable : r = {r:.3f} > 0.5. Réduire dt ou augmenter dx."
        )

    n_steps = int(round(t_max / dt))

    # Profil initial : créneau chaud au centre par défaut
    if initial is None:
        u = np.zeros(nx)
        u[(x > 0.4 * length) & (x < 0.6 * length)] = 1.0
    else:
        u = np.asarray(initial, dtype=float).copy()

    # Conditions aux limites de Dirichlet
    u[0], u[-1] = bc

    frames = [u.copy()]
    times = [0.0]

    for step in range(1, n_steps + 1):
        # Laplacien discret centré : u[i-1] - 2 u[i] + u[i+1]
        laplacian = u[:-2] - 2.0 * u[1:-1] + u[2:]
        u[1:-1] += r * laplacian
        u[0], u[-1] = bc  # on réimpose les bords

        if step % store_every == 0:
            frames.append(u.copy())
            times.append(step * dt)

    return x, np.array(times), np.array(frames)


# --------------------------------------------------------------------------- #
# Solveur 2D
# --------------------------------------------------------------------------- #
def solve_heat_2d(
    *,
    width: float = 1.0,
    height: float = 1.0,
    nx: int = 101,
    ny: int = 101,
    alpha: float = 1.0,
    t_max: float = 0.05,
    dt: float | None = None,
    initial: Array | None = None,
    bc: float = 0.0,
    store_every: int = 5,
    safety: float = 0.9,
) -> tuple[Array, Array, Array, Array]:
    """
    Résout l'équation de la chaleur 2D sur [0, width] x [0, height] avec
    conditions de Dirichlet (température `bc` imposée sur tout le bord).

    Retour : (x, y, times, frames) où frames a la forme (n_frames, ny, nx).
    """
    x = np.linspace(0.0, width, nx)
    y = np.linspace(0.0, height, ny)
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    if dt is None:
        dt = _auto_dt(safety, alpha / dx**2, alpha / dy**2)

    cfl = cfl_number_2d(alpha, dt, dx, dy)
    if cfl > 0.5:
        raise ValueError(
            f"Schéma instable : CFL = {cfl:.3f} > 0.5. Réduire dt ou la résolution."
        )

    rx = alpha * dt / dx**2
    ry = alpha * dt / dy**2
    n_steps = int(round(t_max / dt))

    # Profil initial : source chaude gaussienne au centre par défaut
    if initial is None:
        xx, yy = np.meshgrid(x, y)
        u = np.exp(-(((xx - width / 2) ** 2 + (yy - height / 2) ** 2)) / (0.01))
    else:
        u = np.asarray(initial, dtype=float).copy()

    def apply_bc(field: Array) -> None:
        field[0, :] = field[-1, :] = bc
        field[:, 0] = field[:, -1] = bc

    apply_bc(u)
    frames = [u.copy()]
    times = [0.0]

    for step in range(1, n_steps + 1):
        # Laplacien 2D centré, dérivées x et y séparées (axes : 0 = y, 1 = x)
        lap_x = u[1:-1, :-2] - 2.0 * u[1:-1, 1:-1] + u[1:-1, 2:]
        lap_y = u[:-2, 1:-1] - 2.0 * u[1:-1, 1:-1] + u[2:, 1:-1]
        u[1:-1, 1:-1] += rx * lap_x + ry * lap_y
        apply_bc(u)

        if step % store_every == 0:
            frames.append(u.copy())
            times.append(step * dt)

    return x, y, np.array(times), np.array(frames)
