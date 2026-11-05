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

Files (written by build_files(galaxy), mirroring the MWPotentialHunter24_*.ini set):
  au18  (mock/au18_potential/):  Au18_full = DM + BH + stars + gas; Au18_axi = DM + BH + axisymmetrised
        stars + gas; Au18_baryon_full = stars + gas; Au18_baryon_axi = axisymmetrised stars + gas
  nbody (mock/nbody_potential/): Nbody_full = spherical DM + bulge + disc; Nbody_axi = spherical DM +
        bulge_axi + disc_axi; Nbody_baryon_full = bulge + disc; Nbody_baryon_axi = bulge_axi + disc_axi,
        where each component is axisymmetrised with its own expansion type and grid (bulge Multipole,
        disc CylSpline). The spherical DM halo is never rescaled.
All functions take galaxy="au18" (default) or "nbody".
All functions take any Auriga halo ("au18", "au5", ...) or any N-body snapshot ("nbody",
"nbody_t<T>"). Every Auriga halo has the same component set as Au18 (dm, stars, gas, bh) and is
built the same way; the N-body snapshots share the nbody recipe.
"""
import json
import re

import numpy as np
import agama

from common import ROOT
from build_au18_potential import POT_DIR, CYL_PARAMS
from galaxies import galaxy as _galaxy

agama.setUnits(mass=1, length=1, velocity=1)

TIME_UNIT_GYR = 0.977792221      # 1 kpc / (1 km/s) in Gyr
# `verbose` was added to agama.orbit after 1.0; probe for it so the scripts run on either build.
ORBIT_KW = {"verbose": False} if "verbose" in (agama.orbit.__doc__ or "") else {}
FAMILY = {
    "au18": dict(file_prefix="Au18", fixed=("dm", "bh"), bar=("stars", "gas")),
    "nbody": dict(file_prefix="Nbody", fixed=("dm_spherical",), bar=("bulge", "disc")),
}


def _is_auriga(galaxy):
    return bool(re.fullmatch(r"au\d+", str(galaxy)))


def _family(galaxy):
    """Family recipe for a galaxy name.

    Every nbody snapshot (nbody, nbody_t<T>) shares the nbody recipe, and every Auriga halo shares
    Au18's -- they come out of the same build script with the same components and the same grids.
    Before this, any au<N> other than au18 fell through to the nbody recipe and went looking for
    dm_spherical.ini / bulge / disc, which no Auriga halo has.
    """
    if galaxy in FAMILY:
        return FAMILY[galaxy]
    if _is_auriga(galaxy):
        return dict(FAMILY["au18"], file_prefix=f"Au{str(galaxy)[2:]}")
    return FAMILY["nbody"]


def files(galaxy="au18"):
    g, fp = _galaxy(galaxy), _family(galaxy)["file_prefix"]
    return {k: g["pot_dir"] / f"{fp}_{k}.ini" for k in ("full", "axi", "baryon_full", "baryon_axi")}


def omega_true(galaxy="au18"):
    """True Omega_b in AGAMA's sign convention (negative: clockwise in the MW-like frame)."""
    with open(_galaxy(galaxy)["truth"]) as f:
        return -json.load(f)["Omega_p_kms_kpc"]


def _expansion(galaxy, comp):
    """(type, grid parameters without symmetry) of a component's potential expansion."""
    if _is_auriga(galaxy):
        prm = dict(CYL_PARAMS)      # identical in build_au18_potential and build_auriga_potential
    else:
        from build_nbody_potential import PARAMS
        prm = dict(PARAMS[comp])
    prm.pop("symmetry", None)
    return prm.pop("type"), prm


OMEGA_B_TRUE = omega_true("au18")
FILES = files("au18")


def _concat_ini(paths, out):
    """Write several (possibly multi-section) potential .ini files as one composite file."""
    sections = []
    for p in paths:
        for block in open(p).read().split("[Potential")[1:]:
            sections.append(block.split("]", 1)[1])
    with open(out, "w") as f:
        for i, body in enumerate(sections):
            f.write(f"[Potential{i if i else ''}]" + body.rstrip("\n") + "\n\n")


def build_files(galaxy="au18"):
    """Create the four family files from the component files of the galaxy's build script."""
    if not _is_auriga(galaxy):
        F, pd, fam = files(galaxy), _galaxy(galaxy)["pot_dir"], _family(galaxy)
        axi_paths = []
        for c in fam["bar"]:
            typ, prm = _expansion(galaxy, c)
            prm.pop("mmax", None)
            agama.Potential(type=typ, potential=agama.Potential(str(pd / f"{c}.ini")), symmetry="axisymmetric",
                            **prm).export(str(pd / f"{c}_axi.ini"))
            axi_paths.append(pd / f"{c}_axi.ini")
        _concat_ini([pd / f"{c}.ini" for c in fam["bar"]], F["baryon_full"])
        _concat_ini(axi_paths, F["baryon_axi"])
        _concat_ini([pd / f"{c}.ini" for c in fam["fixed"]] + axi_paths, F["axi"])
        _concat_ini([pd / f"{c}.ini" for c in fam["fixed"] + fam["bar"]], F["full"])
        return
    # Auriga: one merged mmax-8 CylSpline for stars + gas, dm and bh kept as they are. Identical
    # to what Au18 has always done -- POT_DIR/FILES were simply the au18 ones spelled out.
    F, pd = files(galaxy), _galaxy(galaxy)["pot_dir"]
    _concat_ini([pd / "stars.ini", pd / "gas.ini"], F["baryon_full"])
    bary = agama.Potential(file=str(F["baryon_full"]))
    grid = {k: v for k, v in CYL_PARAMS.items() if k not in ("type", "symmetry", "mmax")}
    agama.Potential(type="CylSpline", potential=bary, symmetry="axisymmetric", **grid).export(
        str(F["baryon_axi"]))
    _concat_ini([pd / "dm.ini", pd / "bh.ini", F["baryon_axi"]], F["axi"])
    _concat_ini([pd / "dm.ini", pd / "bh.ini", pd / "stars.ini", pd / "gas.ini"], F["full"])


def load_axi(galaxy="au18"):
    return agama.Potential(file=str(files(galaxy)["axi"]))


def make_potential(bar_frac=1.0, S=1.0, Omega_b=None, phi0=0.0, bartime=None, phi_b=None,
                   pot_axi=None, galaxy="au18"):
    """Barred Au18 potential with rescaled bar strength (bar_frac), length (S) and pattern speed.

    Constant parameters (bartime=None): the bar angle is psi(t) = phi0 + Omega_b t [rad].
    Time-dependent parameters, as in oCen_bar: bartime is an array of times, phi_b the bar angle at
    those times, and bar_frac and S may be arrays of the same length (Omega_b and phi0 are then unused).
    Omega_b = None means the galaxy's true pattern speed.
    """
    F = files(galaxy)
    Omega_b = omega_true(galaxy) if Omega_b is None else Omega_b
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
    pot_axi = pot_axi if pot_axi is not None else load_axi(galaxy)
    # full, scaled, rotated baryonic potential
    pot_bary_scaled = agama.Potential(file=str(F["baryon_full"]), scale=scale_full, rotation=rotation)
    # axisymmetrised, scaled baryonic potential with equal but negative mass
    pot_bary_axi_neg = agama.Potential(file=str(F["baryon_axi"]), scale=scale_axi)
    return agama.Potential(pot_axi, pot_bary_scaled, pot_bary_axi_neg)


def make_static_fast(bar_frac=1.0, S=1.0, phi0=0.0, mmax=4, symmetry="bisymmetric", pot_dm_bh=None,
                     grid=None, galaxy="au18"):
    """Fast, non-rotating version of make_potential for integration in the bar frame with
    agama.orbit(..., Omega=Omega_b). au18: the barred baryonic potential (axisymmetric + scaled bar,
    turned by phi0) is re-expanded into a single CylSpline; nbody: each bar component (bulge, disc)
    is re-expanded with its own expansion type and grid. The fixed components (DM, BH) are kept as
    they are (pot_dm_bh, if given, replaces them)."""
    if not _is_auriga(galaxy):
        pd, fam = _galaxy(galaxy)["pot_dir"], _family(galaxy)
        parts = []
        for c in fam["bar"]:
            pc = agama.Potential(agama.Potential(str(pd / f"{c}_axi.ini")),
                                 agama.Potential(file=str(pd / f"{c}.ini"), scale=[bar_frac * S, S], rotation=phi0),
                                 agama.Potential(file=str(pd / f"{c}_axi.ini"), scale=[-bar_frac * S, S]))
            typ, prm = _expansion(galaxy, c)
            prm["mmax"] = mmax
            if typ == "Multipole":
                prm["lmax"] = max(prm.get("lmax", mmax), mmax)
            parts.append(agama.Potential(type=typ, potential=pc, symmetry=symmetry, **prm))
        if pot_dm_bh is None:
            pot_dm_bh = agama.Potential(*[agama.Potential(str(pd / f"{c}.ini")) for c in fam["fixed"]])
        return agama.Potential(pot_dm_bh, *parts)
    F, pd = files(galaxy), _galaxy(galaxy)["pot_dir"]
    bary = agama.Potential(agama.Potential(file=str(F["baryon_axi"])),
                           agama.Potential(file=str(F["baryon_full"]), scale=[bar_frac * S, S],
                                           rotation=phi0),
                           agama.Potential(file=str(F["baryon_axi"]), scale=[-bar_frac * S, S]))
    grid = grid or {k: v for k, v in CYL_PARAMS.items() if k not in ("type", "symmetry", "mmax")}
    one = agama.Potential(type="CylSpline", potential=bary, mmax=mmax, symmetry=symmetry, **grid)
    if pot_dm_bh is None:
        pot_dm_bh = agama.Potential(agama.Potential(str(pd / "dm.ini")),
                                    agama.Potential(str(pd / "bh.ini")))
    return agama.Potential(pot_dm_bh, one)


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
    import argparse
    from galaxies import add_argument as add_galaxy_argument
    ap = argparse.ArgumentParser(description=__doc__)
    add_galaxy_argument(ap)
    name = ap.parse_args().galaxy
    build_files(name)
    for k, p in files(name).items():
        print(f"{p.name}: {agama.Potential(file=str(p))}")
