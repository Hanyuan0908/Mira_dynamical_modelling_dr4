"""Readers and frame utilities for Auriga z=0 snapshots (Arepo HDF5 + Subfind catalogues).

Snapshot units: length [Mpc/h comoving], mass [1e10 Msun/h], velocity [km/s sqrt(a)].
Everything returned here is converted to kpc, Msun and km/s (physical).
"""
from pathlib import Path

import numpy as np
import h5py


def _chunks(directory, stem):
    files = Path(directory).glob(f"{stem}.*.hdf5")
    return sorted(files, key=lambda p: int(p.stem.split(".")[-1]))


def read_header(snapdir):
    with h5py.File(_chunks(snapdir, "snapshot_*")[0], "r") as f:
        return dict(f["Header"].attrs)


def unit_factors(header):
    h, a = header["HubbleParam"], header["Time"]
    return dict(length=1e3 * a / h, mass=1e10 / h, velocity=np.sqrt(a), box=header["BoxSize"] * 1e3 * a / h)


def main_halo(groupdir, header):
    """Centre (potential minimum), bulk velocity and virial quantities of the most massive FoF group."""
    u = unit_factors(header)
    gfields = ["GroupFirstSub", "Group_M_Crit200", "Group_R_Crit200"]
    sfields = ["SubhaloPos", "SubhaloVel", "SubhaloMassType", "SubhaloHalfmassRadType"]
    g = {k: [] for k in gfields}
    s = {k: [] for k in sfields}
    for fn in _chunks(groupdir, "fof_subhalo_tab_*"):
        with h5py.File(fn, "r") as f:
            if "Group" in f and "GroupFirstSub" in f["Group"]:
                for k in gfields:
                    g[k].append(f["Group"][k][:])
            if "Subhalo" in f and "SubhaloPos" in f["Subhalo"]:
                for k in sfields:
                    s[k].append(f["Subhalo"][k][:])
    g = {k: np.concatenate(v) for k, v in g.items()}
    s = {k: np.concatenate(v) for k, v in s.items()}
    i = int(g["GroupFirstSub"][0])
    return dict(
        pos=s["SubhaloPos"][i] * u["length"],
        vel=s["SubhaloVel"][i],
        mstar=s["SubhaloMassType"][i, 4] * u["mass"],
        rhalf_star=s["SubhaloHalfmassRadType"][i, 4] * u["length"],
        M200=g["Group_M_Crit200"][0] * u["mass"],
        R200=g["Group_R_Crit200"][0] * u["length"],
    )


def read_particles(snapdir, ptype, fields, centre, rmax, header):
    """Read `fields` of PartType`ptype` within `rmax` kpc of `centre`, one chunk at a time.

    Coordinates are returned relative to `centre` in kpc; Masses in Msun (from the mass
    table if the block is absent); Velocities in km/s. Other fields are returned as stored.
    """
    u = unit_factors(header)
    out = {k: [] for k in ["Coordinates", *fields]}
    for fn in _chunks(snapdir, "snapshot_*"):
        with h5py.File(fn, "r") as f:
            grp = f.get(f"PartType{ptype}")
            if grp is None:
                continue
            dx = grp["Coordinates"][:].astype(np.float64) * u["length"] - centre
            dx = (dx + 0.5 * u["box"]) % u["box"] - 0.5 * u["box"]
            sel = np.einsum("ij,ij->i", dx, dx) < rmax**2
            out["Coordinates"].append(dx[sel].astype(np.float32))
            del dx
            for k in fields:
                if k == "Masses" and k not in grp:
                    out[k].append(np.full(sel.sum(), header["MassTable"][ptype] * u["mass"], np.float32))
                    continue
                v = grp[k][:][sel]
                if k == "Masses":
                    v = v * u["mass"]
                elif k == "Velocities":
                    v = v * u["velocity"]
                out[k].append(v.astype(np.float32) if v.dtype.kind == "f" else v)
    return {k: np.concatenate(v) for k, v in out.items()}


def stellar_age_gyr(a_form, header):
    """Age at z=0 [Gyr] from formation scale factor, for the snapshot's flat LCDM cosmology."""
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=100 * header["HubbleParam"], Om0=header["Omega0"])
    a_grid = np.linspace(1e-3, 1.0, 4000)
    lookback = cosmo.lookback_time(1 / a_grid - 1).value
    return np.interp(a_form, a_grid, lookback)


def rotation_to(vector, target):
    """Proper rotation matrix R with R @ vector_hat = target_hat (Rodrigues)."""
    a = vector / np.linalg.norm(vector)
    b = target / np.linalg.norm(target)
    v, c = np.cross(a, b), np.dot(a, b)
    if np.isclose(c, -1.0):            # antiparallel: rotate by pi about any perpendicular axis
        p = np.cross(a, [1.0, 0, 0]) if abs(a[0]) < 0.9 else np.cross(a, [0, 1.0, 0])
        p /= np.linalg.norm(p)
        return 2 * np.outer(p, p) - np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx / (1 + c)


def rotation_z(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


def fourier_m2(x, y, mass, R_edges):
    """Mass-weighted m=2 amplitude A2/A0 and phase (rad, in [-pi/2, pi/2)) in annuli."""
    R, phi = np.hypot(x, y), np.arctan2(y, x)
    idx = np.digitize(R, R_edges) - 1
    ok = (idx >= 0) & (idx < len(R_edges) - 1)
    n = len(R_edges) - 1
    a0 = np.bincount(idx[ok], weights=mass[ok], minlength=n)
    c2 = np.bincount(idx[ok], weights=(mass * np.cos(2 * phi))[ok], minlength=n)
    s2 = np.bincount(idx[ok], weights=(mass * np.sin(2 * phi))[ok], minlength=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        A2 = np.hypot(c2, s2) / a0
    return A2, 0.5 * np.arctan2(s2, c2)
