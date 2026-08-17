# Heat Equation — Finite Difference Simulation

![Python](https://img.shields.io/badge/Python-3.12+-3776AB?logo=python&logoColor=white)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)

Solving the heat (diffusion) equation in 1D and 2D with an explicit finite-difference
scheme (FTCS), and rendering the result as animations. Parameters are exposed on the
command line, so you can experiment without touching the code.

<p align="center">
  <img src="figures/heat_1d.gif" width="48%"/>
  <img src="figures/heat_2d.gif" width="48%"/>
</p>

## Installation

Clone the repo, then use either method.

**With [uv](https://github.com/astral-sh/uv) (recommended).** It reads `uv.lock` and
reproduces the exact environment:

```bash
git clone https://github.com/rayaneferr/heat-equation.git && cd heat-equation
uv run main.py            # creates the env, installs deps and runs, all at once
```

If you don't have uv: `curl -LsSf https://astral.sh/uv/install.sh | sh`
(or `brew install uv` on macOS).

**With pip.** The project uses a standard `pyproject.toml`:

```bash
git clone https://github.com/rayaneferr/heat-equation.git && cd heat-equation
python -m venv .venv && source .venv/bin/activate
pip install .
python main.py
```

pip resolves the latest compatible versions and ignores `uv.lock`, so the pinning is a
bit looser than with uv — fine for this project.

## Usage

```bash
uv run main.py            # 1D + 2D with default parameters
uv run main.py --help     # list every parameter
```

The animations are written to `figures/` (`heat_1d.gif`, `heat_2d.gif`).

Pass flags to change the run — no need to edit the code:

```bash
# 1D only, slower diffusion, finer grid
uv run main.py --dim 1d --alpha 0.5 --nx 201 --t-max 0.2

# 2D only, finer grid, faster GIF, different colormap
uv run main.py --dim 2d --nx 201 --ny 201 --fps 20 --cmap viridis

# A hot bar: left edge held at 1.0, right at 0.0 (heat flows left to right)
uv run main.py --dim 1d --bc-left 1.0 --bc-right 0.0 --t-max 0.3
```

### Parameters

| Flag | Default | Applies to | Meaning |
|---|---|---|---|
| `--dim {1d,2d,both}` | `both` | — | which simulation(s) to run |
| `--alpha` | `1.0` | 1D & 2D | thermal diffusivity (higher = faster diffusion) |
| `--nx` | `121` | 1D & 2D | grid points along x |
| `--fps` | `12` | 1D & 2D | frames per second of the output GIF |
| `--length` | `1.0` | 1D | length of the bar |
| `--t-max` | `0.08` | 1D | final simulation time (1D) |
| `--bc-left` | `0.0` | 1D | fixed temperature at the left edge |
| `--bc-right` | `0.0` | 1D | fixed temperature at the right edge |
| `--store-every-1d` | `50` | 1D | keep 1 frame every N steps |
| `--width` / `--height` | `1.0` | 2D | plate dimensions |
| `--ny` | `121` | 2D | grid points along y |
| `--t-max-2d` | `0.015` | 2D | final simulation time (2D) |
| `--bc` | `0.0` | 2D | fixed temperature on the whole border |
| `--store-every-2d` | `24` | 2D | keep 1 frame every N steps |
| `--cmap` | `inferno` | 2D | Matplotlib colormap |

The time step `dt` is picked automatically to satisfy the CFL stability condition, so a
finer grid (`--nx`/`--ny`) or a larger `--alpha` uses a smaller `dt`. If you pass an
unstable `dt` yourself through the Python API, the solver raises an error rather than
returning garbage.

## Calling the solvers from Python

```python
from src import solve_heat_1d, solve_heat_2d

# 1D — returns (x, times, frames), frames.shape == (n_frames, nx)
x, times, frames = solve_heat_1d(alpha=1.0, t_max=0.08, nx=201)

# 2D — returns (x, y, times, frames), frames.shape == (n_frames, ny, nx)
x, y, times, frames = solve_heat_2d(alpha=1.0, t_max=0.02, bc=0.0)
```

Both accept a custom initial profile through the `initial=` argument.

## The math

The heat equation describes how temperature `u` diffuses over time:

$$
\frac{\partial u}{\partial t} = \alpha \, \nabla^2 u
\qquad\Longrightarrow\qquad
\begin{cases}
\dfrac{\partial u}{\partial t} = \alpha \dfrac{\partial^2 u}{\partial x^2} & \text{(1D)}\\[2ex]
\dfrac{\partial u}{\partial t} = \alpha\!\left(\dfrac{\partial^2 u}{\partial x^2} + \dfrac{\partial^2 u}{\partial y^2}\right) & \text{(2D)}
\end{cases}
$$

With Forward-Time, Central-Space (FTCS) differences, the 1D update is:

$$
u_i^{\,n+1} = u_i^{\,n} + r\left(u_{i+1}^{\,n} - 2u_i^{\,n} + u_{i-1}^{\,n}\right),
\qquad r = \frac{\alpha \, \Delta t}{\Delta x^2}
$$

The scheme is only conditionally stable — it needs the CFL condition to hold:

| Dimension | Stability condition |
|---|---|
| 1D | $r = \dfrac{\alpha \Delta t}{\Delta x^2} \le \dfrac{1}{2}$ |
| 2D | $\alpha \Delta t \left(\dfrac{1}{\Delta x^2} + \dfrac{1}{\Delta y^2}\right) \le \dfrac{1}{2}$ |

Boundaries use Dirichlet conditions (fixed temperature at the edges).

## Project layout

```
heat-equation/
├── README.md
├── pyproject.toml   # dependencies (managed by uv)
├── src/heat.py      # solvers (1D & 2D) + CFL helpers
├── main.py          # CLI demo that generates the animations
└── figures/         # output GIFs
```

A couple of quick checks that the physics behaves: in 1D the central pulse flattens
monotonically (peak 1.0 → 0.18), and in 2D the total heat decreases as the cold
boundaries absorb energy.

## License

MIT — see [LICENSE](LICENSE). Made by [Rayane Ferrat](https://github.com/rayaneferr).
</content>
