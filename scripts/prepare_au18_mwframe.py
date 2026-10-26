"""Extract Au18 (level 4, z=0) star particles and place them in a Milky-Way-like frame.

Frame: astropy-like Galactocentric axes, i.e. the disc rotates clockwise seen from +z
(L_z < 0, so a 'Sun' at x = -R0 moves towards +y), and the bar major axis is rotated to
BAR_ANGLE_DEG from the -x (Sun-GC) direction, with its near end at y > 0.

Output: mock/au18_stars_mwframe.hdf5 (stars within RMAX of the centre) and a bar-measurement
diagnostic figure.
"""
import resource

import numpy as np
import h5py
import matplotlib.pyplot as plt

import auriga as au
from common import (MOCK_DIR, set_style, savefig, BAR_ANGLE_DEG, INK_2, INK_MUTED, ACCENT,
                    BLUE_RAMP)

HALO_DIR = MOCK_DIR / "halo_18"
OUT_FILE = MOCK_DIR / "au18_stars_mwframe.hdf5"
RMAX = 50.0          # kpc, extraction radius
R_DISC = 10.0        # kpc, sphere used for the bulk velocity and angular momentum
Z_BAR = 2.0          # kpc, |z| slab used for the m=2 analysis
DR_BAR = 0.25        # kpc, annulus width for the m=2 profile
R_PROFILE = 12.0     # kpc


def shrinking_sphere(x, m, radii=(10, 5, 3, 2, 1)):
    c = np.zeros(3)
    for r in radii:
        s = np.linalg.norm(x - c, axis=1) < r
        c = np.average(x[s], axis=0, weights=m[s])
    return c


def measure_bar(x, y, z, m):
    """A2 profile, bar half-length (A2 falls below half its peak) and bar phase within it."""
    R_edges = np.arange(0.0, R_PROFILE + 1e-9, DR_BAR)
    slab = np.abs(z) < Z_BAR
    A2, phase = au.fourier_m2(x[slab], y[slab], m[slab], R_edges)
    i_peak = int(np.nanargmax(A2))
    below = np.where(A2[i_peak:] < 0.5 * A2[i_peak])[0]
    R_bar = R_edges[i_peak + below[0]] if len(below) else np.nan
    inside = slab & (np.hypot(x, y) < R_bar)
    _, ph = au.fourier_m2(x[inside], y[inside], m[inside], np.array([0.0, R_bar]))
    return dict(R_edges=R_edges, A2=A2, phase=phase, A2_max=A2[i_peak],
                R_peak=0.5 * (R_edges[i_peak] + R_edges[i_peak + 1]), R_bar=R_bar, phase_bar=ph[0])


def plot_bar_diagnostic(bar_before, bar_after):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5.2, 5.0), sharex=True,
                                   gridspec_kw=dict(hspace=0.08))
    Rm = 0.5 * (bar_after["R_edges"][1:] + bar_after["R_edges"][:-1])
    ax1.plot(Rm, bar_after["A2"], color=BLUE_RAMP[4], lw=2)
    ax1.axhline(0.5 * bar_after["A2_max"], color=INK_MUTED, lw=0.8, ls=":")
    ax1.set_ylabel(r"$A_2/A_0$ ($|z|<2$ kpc)")
    ax2.plot(Rm, np.degrees(bar_before["phase"]), color=INK_MUTED, lw=1.2, label="as aligned")
    ax2.plot(Rm, np.degrees(bar_after["phase"]), color=BLUE_RAMP[4], lw=2,
             label="after rotation to MW angle")
    ax2.axhline(-BAR_ANGLE_DEG, color=ACCENT, lw=1, ls="--")
    ax2.text(0.2, -BAR_ANGLE_DEG + 3, f"target $-{BAR_ANGLE_DEG:.0f}°$",
             color=INK_2, ha="left", va="bottom", fontsize=8)
    ax2.set_ylabel(r"$m=2$ phase [deg]")
    ax2.set_xlabel(r"$R$ [kpc]")
    ax2.set_ylim(-90, 90)
    ax2.legend(loc="lower right", fontsize=8)
    for a in (ax1, ax2):
        a.axvline(bar_after["R_bar"], color=INK_MUTED, lw=0.8, ls="-.")
    ax1.text(bar_after["R_bar"] + 0.15, 0.95 * bar_after["A2_max"],
             rf"$R_{{\rm bar}}$ = {bar_after['R_bar']:.2f} kpc", color=INK_2, fontsize=8, va="top")
    ax1.set_title("Au18 bar: $m=2$ Fourier profile of the stellar mass")
    ax1.set_xlim(0, R_PROFILE)
    savefig(fig, "au18_bar_fourier")


def main():
    set_style()
    hdr = au.read_header(HALO_DIR / "snapdir_127")
    halo = au.main_halo(HALO_DIR / "groups_127", hdr)
    print(f"z = {hdr['Redshift']:.3g}, M200 = {halo['M200']:.3e} Msun, R200 = {halo['R200']:.1f} kpc")

    st = au.read_particles(HALO_DIR / "snapdir_127", 4,
                           ["Velocities", "Masses", "GFM_StellarFormationTime", "GFM_Metallicity",
                            "ParticleIDs"], halo["pos"], RMAX, hdr)
    keep = st["GFM_StellarFormationTime"] > 0          # drop wind particles
    st = {k: v[keep] for k, v in st.items()}
    x = st["Coordinates"].astype(np.float64)
    v = st["Velocities"].astype(np.float64)
    m = st["Masses"].astype(np.float64)

    c = shrinking_sphere(x, m)
    x -= c
    s = np.linalg.norm(x, axis=1) < R_DISC
    vc = np.average(v[s], axis=0, weights=m[s])
    v -= vc

    # disc spin along -z (clockwise seen from +z, as the MW in the astropy frame)
    L = np.sum(m[s, None] * np.cross(x[s], v[s]), axis=0)
    R1 = au.rotation_to(L, np.array([0.0, 0.0, -1.0]))
    x1 = x @ R1.T
    bar0 = measure_bar(x1[:, 0], x1[:, 1], x1[:, 2], m)

    # rotate about z so that the bar major axis sits at -BAR_ANGLE_DEG (near end at y > 0)
    dphi = np.radians(-BAR_ANGLE_DEG) - bar0["phase_bar"]
    R_tot = au.rotation_z(dphi) @ R1
    xf, vf = x @ R_tot.T, v @ R_tot.T
    bar = measure_bar(xf[:, 0], xf[:, 1], xf[:, 2], m)
    print(f"bar: A2_max = {bar['A2_max']:.3f} at R = {bar['R_peak']:.2f} kpc, "
          f"R_bar = {bar['R_bar']:.2f} kpc, phase {np.degrees(bar0['phase_bar']):.1f} -> "
          f"{np.degrees(bar['phase_bar']):.1f} deg")

    Rcyl = np.hypot(xf[:, 0], xf[:, 1])
    ring = (np.abs(Rcyl - 8.0) < 0.5) & (np.abs(xf[:, 2]) < 0.5)
    vphi = (xf[ring, 0] * vf[ring, 1] - xf[ring, 1] * vf[ring, 0]) / Rcyl[ring]
    print(f"mean v_phi (x vy - y vx)/R at R = 8 kpc: {vphi.mean():.1f} km/s (negative = MW sense)")

    age = au.stellar_age_gyr(st["GFM_StellarFormationTime"], hdr)
    with h5py.File(OUT_FILE, "w") as f:
        for i, k in enumerate("xyz"):
            f[k] = xf[:, i].astype(np.float32)
            f["v" + k] = vf[:, i].astype(np.float32)
        f["mass"] = m.astype(np.float32)
        f["age"] = age.astype(np.float32)
        f["Z"] = st["GFM_Metallicity"]
        f["ParticleIDs"] = st["ParticleIDs"]
        f.attrs.update(dict(
            description="Au18 level-4 z=0 star particles (no wind) within RMAX of the centre, "
                        "in a MW-like Galactocentric frame; kpc, km/s, Msun, Gyr.",
            centre_box_kpc=halo["pos"] + c, vcentre_kms=vc, rotation=R_tot,
            note_rotation="x_new = rotation @ (x_box - centre_box_kpc); same for velocities after "
                          "subtracting vcentre_kms",
            RMAX=RMAX, bar_angle_deg=BAR_ANGLE_DEG, bar_half_length_kpc=bar["R_bar"],
            bar_A2_max=bar["A2_max"], M200=halo["M200"], R200=halo["R200"],
            HubbleParam=hdr["HubbleParam"], Redshift=hdr["Redshift"]))
    print(f"wrote {OUT_FILE.name}: {len(m):,} stars, {OUT_FILE.stat().st_size / 1e6:.0f} MB")

    plot_bar_diagnostic(bar0, bar)
    print(f"peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
