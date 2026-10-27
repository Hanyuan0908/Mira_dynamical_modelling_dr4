"""Check: where do the m = 2 kinematic patterns point when the trial bar angle is wrong?

The potential's bar is always along x (bar frame); the stars are rotated by the trial angle phi_b. For trial angles
17, 27 (truth) and 37 deg, the m = 2 coefficients of <v_R>, <v_phi>, sigma_R, sigma_phi (GLS fit per ring with
1, cos/sin 2phi, cos/sin 4phi; R 1.5-5 kpc, |z| < 1) are compared with those at the true angle, giving the rotation of
each pattern. Expected: library patterns stay on x (rotation ~0); data patterns rotate by (true - trial) angle.
Uses the pilot libraries (pilot_04 = 17 deg, pilot_00 = truth, pilot_05 = 37 deg). Prints a table.
"""
import numpy as np
import h5py

import steady_state_likelihood as ssl
from test1_nbody_moments import boot_errors, N_MIN, SEED
from test1_nbody_qmc import MOCK
from test1_fourier_compression import lib, NR, NP, NZ, PHI_C
from test1_fourier_symmetry_terms import design, BOTH

FIELDS = {0: "<v_R>", 1: "<v_phi>", 3: "sigma_R", 4: "sigma_phi"}
LIBS = {17.0: "pilot_04.npz", 27.0: "pilot_00.npz", 37.0: "pilot_05.npz"}
RC = 0.5 * (ssl.CELL_R[1:] + ssl.CELL_R[:-1])
SEL = (RC >= 1.5)[:, None] & (RC <= 5.0)[:, None] & (np.arange(NZ) < 2)[None, :]


def m2_coeffs(val, err, good):
    """Complex m = 2 coefficient (c - i s) per (ring, layer)."""
    Z = np.full((NR, NZ), np.nan + 0j)
    for r in range(NR):
        for z in range(2):
            av = good[r, :, z]
            if av.sum() < 6:
                continue
            F = design(PHI_C[av], BOTH)
            w = 1 / err[r, av, z]**2
            c = np.linalg.solve(F.T @ (w[:, None] * F), F.T @ (w * val[r, av, z]))
            Z[r, z] = c[1] - 1j * c[2]
    return Z


def rotation_deg(Z, Z0):
    """Angle delta by which pattern Z is rotated relative to Z0 (Z = Z0 exp(-2 i delta)), |Z0|^2-weighted."""
    return -0.5 * np.degrees(np.angle(np.nansum((Z * np.conj(Z0))[SEL])))


def fields(obs, lib_m1, lib_s, f):
    N, m1, s, em, es = obs
    good = ((N >= N_MIN) & np.all(em > 0, 1)).reshape(NR, NP, NZ)
    E = np.column_stack([em, es])[:, f].reshape(NR, NP, NZ)
    d = np.column_stack([m1, s])[:, f].reshape(NR, NP, NZ)
    l = np.column_stack([lib_m1, lib_s])[:, f].reshape(NR, NP, NZ)
    return m2_coeffs(d, E, good), m2_coeffs(l, E, good)


def main():
    with h5py.File(MOCK, "r") as fh:
        stars = np.column_stack([fh[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    obs = {p: boot_errors(ssl.rotate_z(stars, np.radians(p)), np.random.default_rng(SEED)) for p in LIBS}
    print("m = 2 pattern rotation relative to the true angle (27 deg), degrees; expected data rotation = trial − 27")
    print(f"{'field':10s} {'trial':>6s} {'data':>8s} {'library':>8s}")
    for f, name in FIELDS.items():
        Zd0, Zl0 = fields(obs[27.0], *lib(LIBS[27.0]), f)
        for p in LIBS:
            Zd, Zl = fields(obs[p], *lib(LIBS[p]), f)
            print(f"{name:10s} {p:6.0f} {rotation_deg(Zd, Zd0):+8.1f} {rotation_deg(Zl, Zl0):+8.1f}")


if __name__ == "__main__":
    main()
