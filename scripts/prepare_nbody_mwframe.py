"""Place an isolated N-body snapshot (mock/Nbody/t_t0_<T>) in the same MW-like frame as Au18.

Same recipe as prepare_au18_mwframe.py: centre on the stars (shrinking spheres), remove the bulk
velocity of the stars within R_DISC, spin the disc to L_z < 0 (clockwise seen from +z), and rotate
the bar (m = 2 phase of the stars within the bar) to -BAR_ANGLE_DEG. Masses are converted to Msun.

Usage: python prepare_nbody_mwframe.py [--time T]      (T is model time, 1 unit = 0.9778 Gyr)

Outputs (<tag> is empty for the fiducial T = 7 and "_t<T>" otherwise):
  mock/nbody<tag>_stars_mwframe.hdf5  disc + bulge particles (comp 0 = disc, 1 = bulge), with the
                                      simulation's own potential and acceleration (pot_sim, acc_sim)
  mock/nbody<tag>_dm_mwframe.hdf5     dark-matter particles, same fields
"""
import argparse
import resource

import numpy as np
import h5py

import auriga as au
import nbody
from common import BAR_ANGLE_DEG, savefig, rel
from galaxies import galaxy, nbody_name, add_time_argument
from prepare_au18_mwframe import shrinking_sphere, measure_bar, plot_bar_diagnostic, R_DISC
from common import set_style


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_time_argument(ap)
    args = ap.parse_args()
    set_style()
    g = galaxy(nbody_name(args.time))
    snap = nbody.snapshot(args.time)
    print(f"reading {rel(snap)}")
    d = nbody.read_nemo(snap)
    comp = nbody.components(d["Mass"])
    pos = d["Position"].astype(np.float64)
    vel = d["Velocity"].astype(np.float64)
    acc = d["Acceleration"].astype(np.float64)
    m = d["Mass"].astype(np.float64) * nbody.MASS_UNIT_MSUN
    st = comp < 2
    print(f"t = {d['Time']:.3f} ({d['Time'] * nbody.TIME_UNIT_GYR:.2f} Gyr); disc {np.sum(comp == 0):,}, "
          f"bulge {np.sum(comp == 1):,}, DM {np.sum(comp == 2):,}; M_disc = {m[comp == 0].sum():.3e}, "
          f"M_bulge = {m[comp == 1].sum():.3e}, M_DM = {m[comp == 2].sum():.3e} Msun")

    c = shrinking_sphere(pos[st], m[st], radii=(20, 10, 5, 3, 2, 1, 0.5))
    pos -= c
    s = st & (np.linalg.norm(pos, axis=1) < R_DISC)
    vc = np.average(vel[s], axis=0, weights=m[s])
    vel -= vc
    L = np.sum(m[s, None] * np.cross(pos[s], vel[s]), axis=0)
    R1 = au.rotation_to(L, np.array([0.0, 0.0, -1.0]))
    x1 = pos[st] @ R1.T
    bar0 = measure_bar(x1[:, 0], x1[:, 1], x1[:, 2], m[st])
    R_tot = au.rotation_z(np.radians(-BAR_ANGLE_DEG) - bar0["phase_bar"]) @ R1
    pf, vf, af = pos @ R_tot.T, vel @ R_tot.T, acc @ R_tot.T
    bar = measure_bar(pf[st, 0], pf[st, 1], pf[st, 2], m[st])
    print(f"centre offset from mean = {c.round(3)} kpc, bulk velocity = {vc.round(2)} km/s")
    print(f"bar: A2_max = {bar['A2_max']:.3f} at R = {bar['R_peak']:.2f} kpc, R_bar = {bar['R_bar']:.2f} kpc, "
          f"phase {np.degrees(bar0['phase_bar']):.1f} -> {np.degrees(bar['phase_bar']):.1f} deg")
    Rc = np.hypot(pf[:, 0], pf[:, 1])
    ring = st & (np.abs(Rc - 8) < 0.5) & (np.abs(pf[:, 2]) < 0.5)
    vphi = -(pf[ring, 0] * vf[ring, 1] - pf[ring, 1] * vf[ring, 0]) / Rc[ring]
    print(f"mean v_phi (along rotation) of stars at R = 8 kpc: {vphi.mean():.1f} km/s")

    attrs = dict(source=str(rel(snap)), snapshot_time=float(args.time), time_model=d["Time"], time_Gyr=d["Time"] * nbody.TIME_UNIT_GYR,
                 centre_model_kpc=c, vcentre_kms=vc, rotation=R_tot,
                 note_rotation="x_new = rotation @ (x_snap - centre); velocities after subtracting vcentre",
                 mass_unit_Msun=nbody.MASS_UNIT_MSUN, bar_angle_deg=BAR_ANGLE_DEG,
                 bar_half_length_kpc=bar["R_bar"], bar_A2_max=bar["A2_max"],
                 note_pot="pot_sim [(km/s)^2] and acc_sim [(km/s)^2/kpc] are the simulation's own (softened) values")
    for fname, sel in ((g["stars"], st), (g["dm"], comp == 2)):
        with h5py.File(fname, "w") as f:
            for i, k in enumerate("xyz"):
                f[k] = pf[sel, i].astype(np.float32)
                f["v" + k] = vf[sel, i].astype(np.float32)
            f["mass"] = m[sel].astype(np.float32)
            f["comp"] = comp[sel]
            f["pot_sim"] = d["Potential"][sel]
            f["acc_sim"] = af[sel].astype(np.float32)
            f.attrs.update(attrs)
        print(f"wrote {rel(fname)}: {sel.sum():,} particles")
    plot_bar_diagnostic(bar0, bar, name=f"{g['prefix']}_bar_fourier", label=g["label"])
    print(f"peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
