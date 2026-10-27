"""Place the isolated N-body galaxy (mock/Nbody/t_t0_7) in the same MW-like frame as Au18.

Same recipe as prepare_au18_mwframe.py: centre on the stars (shrinking spheres), remove the bulk
velocity of the stars within R_DISC, spin the disc to L_z < 0 (clockwise seen from +z), and rotate
the bar (m = 2 phase of the stars within the bar) to -BAR_ANGLE_DEG. Masses are converted to Msun.

Outputs:
  mock/nbody_stars_mwframe.hdf5  disc + bulge particles (comp 0 = disc, 1 = bulge), with the simulation's
                                 own potential and acceleration at each particle (pot_sim, acc_sim)
  mock/nbody_dm_mwframe.hdf5     dark-matter particles, same fields
"""
import resource

import numpy as np
import h5py

import auriga as au
import nbody
from common import BAR_ANGLE_DEG, savefig
from galaxies import galaxy
from prepare_au18_mwframe import shrinking_sphere, measure_bar, plot_bar_diagnostic, R_DISC
from common import set_style


def main():
    set_style()
    g = galaxy("nbody")
    d = nbody.read_nemo()
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

    attrs = dict(source=str(nbody.SNAPSHOT), time_model=d["Time"], time_Gyr=d["Time"] * nbody.TIME_UNIT_GYR,
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
        print(f"wrote {fname.name}: {sel.sum():,} particles")
    plot_bar_diagnostic(bar0, bar, name="nbody_bar_fourier", label=g["label"])
    print(f"peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
