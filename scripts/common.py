"""Shared paths, constants, data loaders and plot style for the Mira DR4 mock tests."""
import os
from pathlib import Path

# AGAMA integrates orbits with OpenMP. Default to every hardware thread unless the caller (or a
# scheduler) has already chosen: on an 8-core/16-thread i7-10700KF, 16 threads integrate at 426
# orbits/s against 356 at 8, i.e. hyperthreading is worth ~20%. Must happen before agama/numpy
# are imported, and common is imported first by every script, so this is the place for it.
os.environ.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 1))

import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
MOCK_DIR = ROOT / "mock"
FIG_DIR = ROOT / "figures"
RESULTS_DIR = ROOT / "results"
EXTERNAL_DIR = ROOT / "external"

# Bulky per-point outputs (trajectories, library sums), kept out of the repo. Portable default:
# ROOT/"runs", which can be a symlink to scratch space; override with $MIRA_DR4_RUNS.
RUNS_DIR = Path(os.environ.get("MIRA_DR4_RUNS") or ROOT / "runs")

MIRA_DR3_FILE = DATA_DIR / "gaia_dr3_lpv_mira.csv"

# The catalogue x, y, z are in the astropy-default Galactocentric frame:
# Sun at x = -R0, y towards l = 90 deg, z towards the NGP.
R0_KPC = 8.122          # GRAVITY Collaboration 2018
Z_SUN_KPC = 0.0208      # Bennett & Bovy 2019
V_SUN_KMS = (12.9, 245.6, 7.78)   # Drimmel & Poggio 2018 (astropy default)

# Bar geometry used only for annotation: angle between the bar major axis and the
# Sun-GC line (near end at l > 0, i.e. y > 0), and long-bar half-length.
BAR_ANGLE_DEG = 27.0    # Wegg & Gerhard 2013
BAR_HALF_LENGTH_KPC = 5.0  # Wegg, Gerhard & Portail 2015

# Plot palette: single-hue blue ramp for densities, orange for overlays,
# neutral inks for text.
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_MUTED = "#8a8984"
ACCENT = "#eb6834"
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CMAP_SEQ = LinearSegmentedColormap.from_list("blue_seq", BLUE_RAMP)
CMAP_SEQ.set_bad("white")
CMAP_SEQ.set_under("white")
# same hue starting from white, for continuous maps (e.g. simulation surface densities)
CMAP_SEQ_W = LinearSegmentedColormap.from_list("blue_seq_w", ["#ffffff"] + BLUE_RAMP)
CMAP_SEQ_W.set_bad("white")
CMAP_SEQ_W.set_under("white")

MIRA_COLUMNS = ["source_id", "l", "b", "period", "age_both", "luminosity_distance",
                "luminosity_distance_err", "x", "y", "z"]


def load_dr3_mira(columns=MIRA_COLUMNS):
    """Load the DR3 Mira catalogue (selected columns only) and add Galactocentric R, phi."""
    df = pd.read_csv(MIRA_DR3_FILE, usecols=columns)
    df = df[np.isfinite(df[["x", "y", "z", "period"]]).all(axis=1)].reset_index(drop=True)
    df["R"] = np.hypot(df["x"], df["y"])
    df["phi"] = np.arctan2(df["y"], df["x"])
    return df


def rel(path):
    """Path relative to ROOT if it is inside the repo, else unchanged (for printing).
    The un-resolved path is tried first, so a symlink inside the repo (e.g. mock/Nbody pointing at
    scratch space) still prints as its repo-relative name."""
    p = Path(path)
    for q in (p, p.resolve()):
        try:
            return q.relative_to(ROOT)
        except ValueError:
            continue
    return p


def set_style():
    mpl.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.edgecolor": INK_2,
        "axes.labelcolor": INK,
        "axes.linewidth": 0.8,
        "xtick.color": INK_2,
        "ytick.color": INK_2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "text.color": INK,
        "legend.frameon": False,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "figure.facecolor": "white",
    })


def savefig(fig, name):
    """Save under figures/<name>.png. `name` may contain "/" to group figures in a sub-directory,
    e.g. savefig(fig, "test2_schemeB/t7_selection") -> figures/test2_schemeB/t7_selection.png."""
    path = FIG_DIR / f"{name}.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    print(f"saved {rel(path)}")
    return path
