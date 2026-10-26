"""Trial potentials for the steady-state test: Au18 with its bar rescaled in strength and length and
rotating with a chosen pattern speed. Follows make_potential() in python_script/oCen_bar/oCen_bar.ipynb.

The axisymmetric part (DM + BH + axisymmetrised stars and gas) is fixed. The non-axisymmetric part
is the full baryonic potential minus its axisymmetrised version:

    Phi(x, t) = Phi_axi(x) + Scale[A, S] Rotate[psi(t)] Phi_bary(x) + Scale[-A, S] Phi_bary,axi(x)
              = Phi_axi(x) + bar_frac * [Phi_bary(R(-psi) x / S) - Phi_bary,axi(x / S)]

using AGAMA's scale modifier, Phi -> (A/S) Phi(x/S), with A = bar_frac * S as in oCen_bar, so that
bar_frac sets the amplitude and S the length of the bar. The rotation modifier turns the bar by
psi(t) = phi0 + Omega_b t, which is anticlockwise for Omega_b > 0 (AGAMA convention). In the MW-like
frame Au18's bar rotates clockwise, so the true Omega_b is negative.

Units: Msun, kpc, km/s; Omega_b in km/s/kpc; time in kpc/(km/s) = TIME_UNIT_GYR Gyr.
Truth: bar_frac = 1, S = 1, Omega_b = OMEGA_B_TRUE, phi0 = 0 (bar at -27 deg at t = 0).

Files (mock/au18_potential/, written by build_files(), mirroring the MWPotentialHunter24_*.ini set):
  Au18_full.ini         DM + BH + stars + gas (ground truth)
  Au18_axi.ini          DM + BH + axisymmetrised stars + gas
  Au18_baryon_full.ini  stars + gas
  Au18_baryon_axi.ini   axisymmetrised stars + gas (CylSpline, same grid, symmetry='axisymmetric')
"""
import json

import numpy as np
import agama

from common import ROOT
from build_au18_potential import POT_DIR, CYL_PARAMS

agama.setUnits(mass=1, length=1, velocity=1)

TIME_UNIT_GYR = 0.977792221      # 1 kpc / (1 km/s) in Gyr
with open(ROOT / "results" / "au18_ground_truth.json") as _f:
    OMEGA_B_TRUE = -json.load(_f)["Omega_p_kms_kpc"]      # clockwise in the MW-like frame

FILES = {k: POT_DIR / f"Au18_{k}.ini" for k in ("full", "axi", "baryon_full", "baryon_axi")}


def _concat_ini(paths, out):
    """Write several (possibly multi-section) potential .ini files as one composite file."""
    sections = []
    for p in paths:
        for block in open(p).read().split("[Potential")[1:]:
            sections.append(block.split("]", 1)[1])
    with open(out, "w") as f:
        for i, body in enumerate(sections):
            f.write(f"[Potential{i if i else ''}]" + body.rstrip("\n") + "\n\n")


def build_files():
    """Create the four files above from the component files of build_au18_potential.py."""
    _concat_ini([POT_DIR / "stars.ini", POT_DIR / "gas.ini"], FILES["baryon_full"])
    bary = agama.Potential(file=str(FILES["baryon_full"]))
    grid = {k: v for k, v in CYL_PARAMS.items() if k not in ("type", "symmetry", "mmax")}
    agama.Potential(type="CylSpline", potential=bary, symmetry="axisymmetric", **grid).export(
        str(FILES["baryon_axi"]))
    _concat_ini([POT_DIR / "dm.ini", POT_DIR / "bh.ini", FILES["baryon_axi"]], FILES["axi"])
    _concat_ini([POT_DIR / "dm.ini", POT_DIR / "bh.ini", POT_DIR / "stars.ini", POT_DIR / "gas.ini"],
                FILES["full"])


def load_axi():
    return agama.Potential(file=str(FILES["axi"]))


def make_potential(bar_frac=1.0, S=1.0, Omega_b=OMEGA_B_TRUE, phi0=0.0, bartime=None, phi_b=None,
                   pot_axi=None):
    """Barred Au18 potential with rescaled bar strength (bar_frac), length (S) and pattern speed.

    Constant parameters (bartime=None): the bar angle is psi(t) = phi0 + Omega_b t [rad].
    Time-dependent parameters, as in oCen_bar: bartime is an array of times, phi_b the bar angle at
    those times, and bar_frac and S may be arrays of the same length (Omega_b and phi0 are then unused).
    """
    if bartime is None:
        scale_full, scale_axi = [bar_frac * S, S], [-bar_frac * S, S]
        rotation = [[0.0, phi0], [1.0, phi0 + Omega_b]]
    else:
        bartime = np.asarray(bartime, dtype=float)
        A = np.broadcast_to(bar_frac, bartime.shape) * np.broadcast_to(S, bartime.shape)
        S_t = np.broadcast_to(S, bartime.shape)
        scale_full = np.column_stack([bartime, A, S_t])
        scale_axi = np.column_stack([bartime, -A, S_t])
        rotation = np.column_stack([bartime, phi_b])
    pot_axi = pot_axi if pot_axi is not None else load_axi()
    # full, scaled, rotated baryonic potential
    pot_bary_scaled = agama.Potential(file=str(FILES["baryon_full"]), scale=scale_full, rotation=rotation)
    # axisymmetrised, scaled baryonic potential with equal but negative mass
    pot_bary_axi_neg = agama.Potential(file=str(FILES["baryon_axi"]), scale=scale_axi)
    return agama.Potential(pot_axi, pot_bary_scaled, pot_bary_axi_neg)


def bar_angle(t, Omega_b=OMEGA_B_TRUE, phi0=0.0):
    """Rotation of the bar relative to t = 0 [rad] for constant Omega_b."""
    return phi0 + Omega_b * np.asarray(t)


def to_bar_frame(posvel, t, Omega_b=OMEGA_B_TRUE, phi0=0.0):
    """Rotate inertial positions and velocities (N x 6) at time t into the frame co-rotating with
    the bar, in which the bar stays at its t = 0 orientation. Velocities stay inertial."""
    a = -bar_angle(t, Omega_b, phi0)
    c, s = np.cos(a), np.sin(a)
    out = np.array(posvel, dtype=float, copy=True)
    for i in (0, 3):
        x, y = posvel[..., i], posvel[..., i + 1]
        out[..., i], out[..., i + 1] = c * x - s * y, s * x + c * y
    return out


if __name__ == "__main__":
    build_files()
    for k, p in FILES.items():
        print(f"{p.name}: {agama.Potential(file=str(p))}")
