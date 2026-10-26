"""Ground-truth potential of Au18 (z = 0) in the MW-like frame of mock/au18_stars_mwframe.hdf5.

Components (AGAMA units: Msun, kpc, km/s):
  dark matter : Multipole, spherical symmetry (all high-res DM within R_SEL)
  stars       : CylSpline, no symmetry (star particles of au18_stars_mwframe.hdf5, wind excluded)
  gas         : CylSpline, no symmetry (gas cells; CylSpline only uses those inside its grid)
  black hole  : Plummer sphere with the gravitational softening as scale radius

Hyperparameters (see PROGRESS.md for the reasoning):
  Multipole follows AGAMA's py/example_gizmo_snapshot.py, the AGAMA example for a cosmological
  zoom-in galaxy (FIRE): rmin = 0.1 kpc, rmax ~ 500 kpc, default gridSizeR = 25.
  CylSpline follows the grids in python_script/Nexus/get_potential*.py (gridSizeR = 30,
  Rmin = 0.15, zmin = 0.1, Rmax = 50, mmax = 8), with zmax raised from 20 to Rmax so the grid
  encloses the whole sphere r < 50 kpc. CylSpline ignores mass outside its grid, and the
  excluded outer gas is then nearly spherical; its force is checked below.

Output: mock/au18_potential/{dm,stars,gas,bh}.ini and summary.json
"""
import json
import resource
import time

import numpy as np
import h5py
import agama

import auriga as au
from common import MOCK_DIR
from prepare_au18_mwframe import OUT_FILE as STARS_FILE, HALO_DIR

agama.setUnits(mass=1, length=1, velocity=1)

POT_DIR = MOCK_DIR / "au18_potential"
R_SEL = 500.0                 # kpc, selection radius for DM and gas
BH_SOFTENING = 0.369          # kpc, Auriga level-4 softening at z = 0 (Grand et al. 2017)
DM_PARAMS = dict(type="Multipole", symmetry="spherical", gridSizeR=25, rmin=0.1, rmax=R_SEL)
CYL_PARAMS = dict(type="CylSpline", symmetry="none", mmax=8, gridSizeR=30, gridSizez=30,
                  Rmin=0.15, Rmax=50.0, zmin=0.1, zmax=50.0)


def load_frame():
    with h5py.File(STARS_FILE, "r") as f:
        return {k: f.attrs[k] for k in ("centre_box_kpc", "vcentre_kms", "rotation")}


def read_mwframe(ptype, fields, rmax, frame, header):
    """Particles of PartType`ptype` within rmax, positions (and velocities) in the MW-like frame."""
    d = au.read_particles(HALO_DIR / "snapdir_127", ptype, fields, frame["centre_box_kpc"], rmax, header)
    d["pos"] = d.pop("Coordinates").astype(np.float64) @ frame["rotation"].T
    if "Velocities" in d:
        d["vel"] = (d.pop("Velocities").astype(np.float64) - frame["vcentre_kms"]) @ frame["rotation"].T
    return d


def load_potential():
    """Components {dm, stars, gas, bh} and their sum, from the exported coefficient files."""
    comps = {k: agama.Potential(str(POT_DIR / f"{k}.ini")) for k in ("dm", "stars", "gas", "bh")}
    return comps, agama.Potential(*comps.values())


def vcirc(pot, R, nphi=72):
    """Circular speed from the azimuthally averaged radial force at z = 0 (NaN where v_c^2 < 0)."""
    phi = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    RR, PP = np.meshgrid(R, phi, indexing="ij")
    pts = np.column_stack([(RR * np.cos(PP)).ravel(), (RR * np.sin(PP)).ravel(), np.zeros(RR.size)])
    f = pot.force(pts)
    vc2 = -(f[:, 0] * pts[:, 0] + f[:, 1] * pts[:, 1]).reshape(RR.shape).mean(axis=1)
    with np.errstate(invalid="ignore"):
        return np.sqrt(vc2)


def timed(label, fn):
    t0 = time.time()
    out = fn()
    print(f"  {label}: {time.time() - t0:.1f} s")
    return out


def main():
    POT_DIR.mkdir(exist_ok=True)
    hdr = au.read_header(HALO_DIR / "snapdir_127")
    frame = load_frame()
    summary = dict(units="Msun, kpc, km/s", frame="MW-like frame of au18_stars_mwframe.hdf5",
                   dm_params=DM_PARAMS, cyl_params=CYL_PARAMS, R_sel=R_SEL)

    print("dark matter")
    dm = read_mwframe(1, ["Masses"], R_SEL, frame, hdr)
    pot_dm = timed("Multipole", lambda: agama.Potential(particles=(dm["pos"], dm["Masses"]), **DM_PARAMS))
    summary["dm"] = dict(N=len(dm["Masses"]), mass=float(dm["Masses"].sum()))
    del dm

    print("stars")
    with h5py.File(STARS_FILE, "r") as f:
        pos = np.column_stack([f[k][:].astype(np.float64) for k in "xyz"])
        m = f["mass"][:].astype(np.float64)
    inside = (np.hypot(pos[:, 0], pos[:, 1]) < CYL_PARAMS["Rmax"]) & (np.abs(pos[:, 2]) < CYL_PARAMS["zmax"])
    pot_st = timed("CylSpline", lambda: agama.Potential(particles=(pos[inside], m[inside]), **CYL_PARAMS))
    summary["stars"] = dict(N=int(inside.sum()), mass=float(m[inside].sum()),
                            mass_outside_grid=float(m[~inside].sum()))
    del pos, m

    print("gas")
    gas = read_mwframe(0, ["Masses"], R_SEL, frame, hdr)
    pos, m = gas["pos"], gas["Masses"].astype(np.float64)
    inside = (np.hypot(pos[:, 0], pos[:, 1]) < CYL_PARAMS["Rmax"]) & (np.abs(pos[:, 2]) < CYL_PARAMS["zmax"])
    pot_gas = timed("CylSpline", lambda: agama.Potential(particles=(pos[inside], m[inside]), **CYL_PARAMS))
    # force at the solar circle from the gas the CylSpline grid leaves out (spherical expansion)
    pot_out = agama.Potential(type="Multipole", particles=(pos[~inside], m[~inside]), lmax=6,
                              symmetry="none", rmin=1.0, rmax=R_SEL)
    summary["gas"] = dict(N=int(inside.sum()), mass=float(m[inside].sum()),
                          mass_outside_grid=float(m[~inside].sum()))
    del gas, pos, m

    print("black hole")
    bh = read_mwframe(5, ["Masses"], 5.0, frame, hdr)
    i = np.argmin(np.linalg.norm(bh["pos"], axis=1))
    pot_bh = agama.Potential(type="Plummer", mass=float(bh["Masses"][i]), scaleRadius=BH_SOFTENING)
    summary["bh"] = dict(mass=float(bh["Masses"][i]), offset_kpc=float(np.linalg.norm(bh["pos"][i])),
                         scaleRadius=BH_SOFTENING)

    comps = dict(dm=pot_dm, stars=pot_st, gas=pot_gas, bh=pot_bh)
    for name in ("dm", "stars", "gas"):
        comps[name].export(str(POT_DIR / f"{name}.ini"))
    # export() does not store the parameters of analytic potentials, so write the Plummer by hand
    with open(POT_DIR / "bh.ini", "w") as f:
        f.write(f"[Potential]\ntype=Plummer\nmass={summary['bh']['mass']!r}\nscaleRadius={BH_SOFTENING!r}\n")
    comps_loaded, total = load_potential()
    for name, p in comps.items():                  # the exported files must reproduce the potentials
        pt = np.array([[0.5, 0.3, 0.1], [8.0, -1.0, 0.2]])
        assert np.allclose(comps_loaded[name].potential(pt), p.potential(pt), rtol=1e-6), name

    # force check at R = 8 kpc, averaged over azimuth
    phi = np.linspace(0, 2 * np.pi, 72, endpoint=False)
    pts = np.column_stack([8 * np.cos(phi), 8 * np.sin(phi), 0 * phi])

    def vc2(p):
        f = p.force(pts)
        return float(np.mean(-(f[:, 0] * pts[:, 0] + f[:, 1] * pts[:, 1])))   # R * (-dPhi/dR)
    summary["vc2_R8"] = {k: vc2(p) for k, p in comps.items()} | {"total": vc2(total),
                                                                   "gas_outside_grid": vc2(pot_out)}
    with open(POT_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps({k: summary[k] for k in ("dm", "stars", "gas", "bh")}, indent=1))
    v = summary["vc2_R8"]
    print("v_c(R=8 kpc) [km/s]: total {:.1f}; dm {:.1f}, stars {:.1f}, gas {:.1f}, bh {:.1f}".format(
        *(np.sqrt(max(v[k], 0)) for k in ("total", "dm", "stars", "gas", "bh"))))
    print(f"gas outside the CylSpline grid: v_c^2 contribution at R = 8 kpc = {v['gas_outside_grid']:.1f} "
          f"(km/s)^2, i.e. {100 * v['gas_outside_grid'] / v['total']:.2f}% of the total")
    print(f"peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
