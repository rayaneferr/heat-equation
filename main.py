"""
Équation de la chaleur — démonstration animée (1D & 2D).

Tous les paramètres physiques et numériques sont réglables en ligne de commande
(voir `uv run main.py --help`). Exemples :

    uv run main.py                                  # 1D + 2D, valeurs par défaut
    uv run main.py --dim 1d --alpha 0.5 --t-max 0.2 # juste la 1D, diffusion plus lente
    uv run main.py --dim 2d --nx 201 --fps 20       # 2D plus fine, animation plus rapide

Les animations sont écrites dans le dossier `figures/`.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from src import solve_heat_1d, solve_heat_2d

FIGURES = Path(__file__).parent / "figures"


def animate_1d(args: argparse.Namespace) -> None:
    """Diffusion 1D : un créneau chaud s'aplatit au fil du temps."""
    x, times, frames = solve_heat_1d(
        length=args.length,
        nx=args.nx,
        alpha=args.alpha,
        t_max=args.t_max,
        bc=(args.bc_left, args.bc_right),
        store_every=args.store_every_1d,
    )

    fig, ax = plt.subplots(figsize=(8, 4.5))
    (line,) = ax.plot(x, frames[0], color="crimson", lw=2)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("position x")
    ax.set_ylabel("température u(x, t)")
    ax.grid(alpha=0.3)
    title = ax.set_title("")

    def update(i: int):
        line.set_ydata(frames[i])
        title.set_text(f"Équation de la chaleur 1D — t = {times[i]:.4f}")
        return line, title

    anim = FuncAnimation(fig, update, frames=len(frames), interval=60, blit=False)
    out = FIGURES / "heat_1d.gif"
    anim.save(out, writer="pillow", fps=args.fps)
    plt.close(fig)
    print(f"✅ {out}  ({len(frames)} frames)")


def animate_2d(args: argparse.Namespace) -> None:
    """Diffusion 2D : une tache de chaleur gaussienne se répand sur une plaque."""
    x, y, times, frames = solve_heat_2d(
        width=args.width,
        height=args.height,
        nx=args.nx,
        ny=args.ny,
        alpha=args.alpha,
        t_max=args.t_max_2d,
        bc=args.bc,
        store_every=args.store_every_2d,
    )

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(
        frames[0],
        origin="lower",
        extent=(x[0], x[-1], y[0], y[-1]),
        cmap=args.cmap,
        vmin=0.0,
        vmax=1.0,
    )
    fig.colorbar(im, ax=ax, label="température")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    title = ax.set_title("")

    def update(i: int):
        im.set_data(frames[i])
        title.set_text(f"Équation de la chaleur 2D — t = {times[i]:.4f}")
        return im, title

    anim = FuncAnimation(fig, update, frames=len(frames), interval=80, blit=False)
    out = FIGURES / "heat_2d.gif"
    anim.save(out, writer="pillow", fps=args.fps)
    plt.close(fig)
    print(f"✅ {out}  ({len(frames)} frames)")


def build_parser() -> argparse.ArgumentParser:
    """Définit tous les paramètres réglables depuis le terminal."""
    p = argparse.ArgumentParser(
        description="Simulation de l'équation de la chaleur (1D & 2D).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Quelle(s) simulation(s) lancer
    p.add_argument(
        "--dim",
        choices=("1d", "2d", "both"),
        default="both",
        help="dimension à simuler",
    )

    # Paramètres physiques communs
    phys = p.add_argument_group("physique")
    phys.add_argument("--alpha", type=float, default=1.0, help="diffusivité thermique")

    # Paramètres 1D
    g1 = p.add_argument_group("1D")
    g1.add_argument("--length", type=float, default=1.0, help="longueur de la barre")
    g1.add_argument("--t-max", type=float, default=0.08, help="temps final (1D)")
    g1.add_argument("--bc-left", type=float, default=0.0, help="température bord gauche")
    g1.add_argument("--bc-right", type=float, default=0.0, help="température bord droit")
    g1.add_argument("--store-every-1d", type=int, default=50, help="frames : 1 pas sur N")

    # Paramètres 2D
    g2 = p.add_argument_group("2D")
    g2.add_argument("--width", type=float, default=1.0, help="largeur de la plaque")
    g2.add_argument("--height", type=float, default=1.0, help="hauteur de la plaque")
    g2.add_argument("--ny", type=int, default=121, help="points de grille en y")
    g2.add_argument("--t-max-2d", type=float, default=0.015, help="temps final (2D)")
    g2.add_argument("--bc", type=float, default=0.0, help="température imposée au bord")
    g2.add_argument("--store-every-2d", type=int, default=24, help="frames : 1 pas sur N")
    g2.add_argument("--cmap", default="inferno", help="palette de couleurs Matplotlib")

    # Grille (partagé) & rendu
    p.add_argument("--nx", type=int, default=121, help="points de grille en x")
    p.add_argument("--fps", type=int, default=12, help="images par seconde des GIF")

    return p


def main() -> None:
    args = build_parser().parse_args()
    FIGURES.mkdir(exist_ok=True)

    print(f"Simulation (dim={args.dim}, alpha={args.alpha})…")
    if args.dim in ("1d", "both"):
        animate_1d(args)
    if args.dim in ("2d", "both"):
        animate_2d(args)
    print("Terminé.")


if __name__ == "__main__":
    main()
