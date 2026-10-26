"""Benchmarks of the Au18 ground-truth potential (scripts/build_au18_potential.py).

1. Stellar density: density of the stellar CylSpline (Poisson-consistent with its potential) vs
   the star particles. Face-on and edge-on projections, Sigma(R), rho(z), m = 2 amplitude and phase.
2. Rotation curve: v_c from the total potential, and the stellar <v_phi> it implies through the
   axisymmetric Jeans equation (asymmetric drift), compared with the measured stellar <v_phi>.
   Gas <v_phi> is shown as an independent near-circular tracer.
3. Extra: AGAMA potential vs the potential Arepo stores for each star particle (up to a constant).
"""
import json
import resource

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import auriga as au
from common import (ROOT, set_style, savefig, CMAP_SEQ_W, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP)
from prepare_au18_mwframe import OUT_FILE, HALO_DIR
from build_au18_potential import load_potential, vcirc, read_mwframe, CYL_PARAMS
from plot_au18_velocity_maps import CMAP_DIV

RESULTS_FILE = ROOT / "results" / "au18_potential_checks.json"
MAP_LIM, MAP_PIX = 12.0, 0.25          # kpc, face-on map
EDGE_ZLIM, EDGE_ZPIX = 4.0, 0.1        # kpc, edge-on map
Z_SLAB = 0.5                           # kpc, midplane slab for kinematics
DR_KIN = 0.5                           # kpc, radial bins for kinematics


def los_grid(zmax=CYL_PARAMS["zmax"], n=60):
    """Symmetric line-of-sight grid (log-spaced away from the plane) and trapezoid weights."""
    h = np.concatenate([[0], np.logspace(-2.5, np.log10(zmax), n)])
    g = np.concatenate([-h[:0:-1], h])
    w = np.zeros_like(g)
    w[1:] += 0.5 * np.diff(g)
    w[:-1] += 0.5 * np.diff(g)
    return g, w


def project(pot, a, b, axis):
    """Sigma(a, b) = integral of rho along `axis` ('z' for face-on x-y, 'y' for edge-on x-z)."""
    g, w = los_grid()
    A, B, G = (np.repeat(a[:, None], len(g), 1), np.repeat(b[:, None], len(g), 1),
               np.repeat(g[None, :], len(a), 0))
    xyz = np.column_stack([A.ravel(), B.ravel(), G.ravel()] if axis == "z" else
                          [A.ravel(), G.ravel(), B.ravel()])
    return (pot.density(xyz).reshape(A.shape) * w).sum(axis=1)


def density_benchmark(pot_st, x, y, z, m):
    out = {}
    e = np.arange(-MAP_LIM, MAP_LIM + 1e-9, MAP_PIX)
    c = 0.5 * (e[1:] + e[:-1])
    Hp, _, _ = np.histogram2d(x, y, bins=[e, e], weights=m)
    out["face_p"] = Hp.T / MAP_PIX**2
    X, Y = np.meshgrid(c, c)
    out["face_m"] = project(pot_st, X.ravel(), Y.ravel(), "z").reshape(X.shape)
    out["face_N"] = np.histogram2d(x, y, bins=[e, e])[0].T
    ez = np.arange(-EDGE_ZLIM, EDGE_ZLIM + 1e-9, EDGE_ZPIX)
    cz = 0.5 * (ez[1:] + ez[:-1])
    He, _, _ = np.histogram2d(x, z, bins=[e, ez], weights=m)
    out["edge_p"] = He.T / (MAP_PIX * EDGE_ZPIX)
    XE, ZE = np.meshgrid(c, cz)
    out["edge_m"] = project(pot_st, XE.ravel(), ZE.ravel(), "y").reshape(XE.shape)
    out.update(e=e, c=c, ez=ez, cz=cz)

    # Sigma(R) and m = 2 Fourier terms on a polar grid
    R_e = np.arange(0.0, 20.0 + 1e-9, 0.25)
    R_c = 0.5 * (R_e[1:] + R_e[:-1])
    R = np.hypot(x, y)
    out["SigR_p"] = np.histogram(R, R_e, weights=m)[0] / (np.pi * np.diff(R_e**2))
    phi = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    RR, PP = np.meshgrid(R_c, phi, indexing="ij")
    S = project(pot_st, (RR * np.cos(PP)).ravel(), (RR * np.sin(PP)).ravel(), "z").reshape(RR.shape)
    out["SigR_m"] = S.mean(axis=1)
    c2 = (S * np.cos(2 * PP)).mean(axis=1)
    s2 = (S * np.sin(2 * PP)).mean(axis=1)
    out["A2_m"] = np.hypot(c2, s2) / S.mean(axis=1)       # same convention as au.fourier_m2
    out["ph2_m"] = 0.5 * np.arctan2(s2, c2)
    out["A2_p"], out["ph2_p"] = au.fourier_m2(x, y, m, R_e)
    out["R_c"] = R_c

    # rho(z) at a few radii, annulus +/- 0.25 kpc, folded in |z|
    out["rho_z"] = {}
    z_e = np.arange(0, 3.0 + 1e-9, 0.1)
    z_c = 0.5 * (z_e[1:] + z_e[:-1])
    for R0 in (2.0, 5.0, 8.0):
        s = np.abs(R - R0) < 0.25
        vol = np.pi * ((R0 + 0.25)**2 - (R0 - 0.25)**2) * 2 * np.diff(z_e)
        rp = np.histogram(np.abs(z[s]), z_e, weights=m[s])[0] / vol
        ph = np.linspace(0, 2 * np.pi, 32, endpoint=False)
        P, Z = np.meshgrid(ph, z_c, indexing="ij")
        rm = pot_st.density(np.column_stack([(R0 * np.cos(P)).ravel(), (R0 * np.sin(P)).ravel(),
                                             Z.ravel()])).reshape(P.shape).mean(axis=0)
        out["rho_z"][R0] = (z_c, rp, rm)
    return out


def jeans_vphi(R_c, x, y, z, vx, vy, vz, m, vc):
    """Stellar <v_phi> predicted by the axisymmetric Jeans equation in the midplane:
    v_c^2 - <v_phi>^2 = sigma_phi^2 - sigma_R^2 - sigma_R^2 dln(nu sigma_R^2)/dlnR - R d<v_R v_z>/dz."""
    R = np.hypot(x, y)
    vR = (x * vx + y * vy) / R
    vphi = -(x * vy - y * vx) / R                    # positive along (clockwise) rotation
    idx = np.digitize(R, np.append(R_c - DR_KIN / 2, R_c[-1] + DR_KIN / 2)) - 1
    n = len(R_c)
    out = {k: np.full(n, np.nan) for k in ("vphi", "sR", "sphi", "nu", "dvRvz_dz", "vphi_pred", "N")}
    for i in range(n):
        s = idx == i
        if s.sum() < 200:
            continue
        mid = s & (np.abs(z) < Z_SLAB)
        w = m[mid]
        mu_R, mu_phi = np.average(vR[mid], weights=w), np.average(vphi[mid], weights=w)
        out["vphi"][i] = mu_phi
        out["sR"][i] = np.sqrt(np.average((vR[mid] - mu_R)**2, weights=w))
        out["sphi"][i] = np.sqrt(np.average((vphi[mid] - mu_phi)**2, weights=w))
        out["nu"][i] = w.sum() / (np.pi * ((R_c[i] + DR_KIN / 2)**2 - (R_c[i] - DR_KIN / 2)**2) * 2 * Z_SLAB)
        out["N"][i] = mid.sum()
        thick = s & (np.abs(z) < 2 * Z_SLAB)       # slope of <v_R v_z> in z through the midplane
        wt, zt = m[thick], z[thick]
        out["dvRvz_dz"][i] = np.sum(wt * zt * vR[thick] * vz[thick]) / np.sum(wt * zt**2)
    lnR = np.log(R_c)
    dln = np.gradient(np.log(out["nu"] * out["sR"]**2), lnR)
    drift = out["sphi"]**2 - out["sR"]**2 - out["sR"]**2 * dln - R_c * out["dvRvz_dz"]
    out["drift"] = drift
    with np.errstate(invalid="ignore"):
        out["vphi_pred"] = np.sqrt(vc**2 - drift)
    return out


def main():
    set_style()
    comps, total = load_potential()
    with h5py.File(OUT_FILE, "r") as f:
        x, y, z, vx, vy, vz, m = (f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz", "mass"))
        R_bar = float(f.attrs["bar_half_length_kpc"])
    results = {}

    # ---------- 1. stellar density ----------
    d = density_benchmark(comps["stars"], x, y, z, m)
    good = d["face_N"] >= 100
    ratio = d["face_m"] / d["face_p"]
    inner = good & (np.hypot(*np.meshgrid(d["c"], d["c"])) < 10)
    results["face_on_model_over_particles_median"] = float(np.median(ratio[inner]))
    results["face_on_abs_log_ratio_68pct"] = float(np.percentile(np.abs(np.log10(ratio[inner])), 68))
    sel = (d["R_c"] > 0.5) & (d["R_c"] < 15)
    results["SigR_ratio_range_0.5_15kpc"] = [float(np.min(d["SigR_m"][sel] / d["SigR_p"][sel])),
                                             float(np.max(d["SigR_m"][sel] / d["SigR_p"][sel]))]
    rb = (d["R_c"] > 1) & (d["R_c"] < 4)
    results["A2_peak_particles_model"] = [float(np.nanmax(d["A2_p"])), float(np.nanmax(d["A2_m"]))]
    results["bar_phase_1_4kpc_particles_model_deg"] = [float(np.degrees(np.median(d["ph2_p"][rb]))),
                                                       float(np.degrees(np.median(d["ph2_m"][rb])))]

    fig = plt.figure(figsize=(15.5, 9.6))
    gs = fig.add_gridspec(2, 3, hspace=0.32, wspace=0.42)
    vmax = d["face_p"].max()
    norm = LogNorm(vmin=vmax / 1e3, vmax=vmax)
    ax = fig.add_subplot(gs[0, 0])
    im = ax.pcolormesh(d["e"], d["e"], np.ma.masked_equal(d["face_p"], 0), cmap=CMAP_SEQ_W, norm=norm, rasterized=True)
    lev = 10 ** np.arange(np.log10(vmax) - 2.75, np.log10(vmax), 0.25)
    ax.contour(d["c"], d["c"], d["face_p"], levels=lev, colors=INK_2, linewidths=0.5)
    ax.contour(d["c"], d["c"], d["face_m"], levels=lev, colors=ACCENT, linewidths=0.9, linestyles="--")
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$ [kpc]")
    ax.set_ylabel(r"$y$ [kpc]")
    ax.set_title("Face-on Σ: particles (grey) vs potential (orange)", color=INK_2, fontsize=9.5)
    cb = fig.colorbar(im, cax=ax.inset_axes([1.03, 0, 0.045, 1]))
    cb.set_label(r"$\Sigma_\star$ [M$_\odot$ kpc$^{-2}$]", fontsize=9)

    ax = fig.add_subplot(gs[0, 1])
    im = ax.pcolormesh(d["e"], d["e"], np.ma.masked_where(~good, ratio - 1), cmap=CMAP_DIV, vmin=-0.2,
                       vmax=0.2, rasterized=True)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$ [kpc]")
    ax.set_title(r"Face-on $\Sigma_{\rm pot}/\Sigma_{\rm part} - 1$ (pixels with N ≥ 100)", color=INK_2, fontsize=9.5)
    cb = fig.colorbar(im, cax=ax.inset_axes([1.03, 0, 0.045, 1]), extend="both")
    cb.set_label("fractional residual", fontsize=9)

    ax = fig.add_subplot(gs[0, 2])
    vmax_e = d["edge_p"].max()
    lev_e = 10 ** np.arange(np.log10(vmax_e) - 2.75, np.log10(vmax_e), 0.25)
    ax.pcolormesh(d["e"], d["ez"], np.ma.masked_equal(d["edge_p"], 0), cmap=CMAP_SEQ_W,
                  norm=LogNorm(vmin=vmax_e / 1e3, vmax=vmax_e), rasterized=True)
    ax.contour(d["c"], d["cz"], d["edge_p"], levels=lev_e, colors=INK_2, linewidths=0.5)
    ax.contour(d["c"], d["cz"], d["edge_m"], levels=lev_e, colors=ACCENT, linewidths=0.9, linestyles="--")
    ax.set_aspect(1.5)
    ax.set_xlabel(r"$x$ [kpc]")
    ax.set_ylabel(r"$z$ [kpc] ($\times$1.5)")
    ax.set_title("Edge-on Σ (x–z): particles vs potential", color=INK_2, fontsize=9.5)

    ax = fig.add_subplot(gs[1, 0])
    ax.semilogy(d["R_c"], d["SigR_p"], color=INK_2, lw=2, label="particles")
    ax.semilogy(d["R_c"], d["SigR_m"], color=ACCENT, lw=1.3, ls="--", label="potential")
    ax.set_xlabel(r"$R$ [kpc]")
    ax.set_ylabel(r"$\Sigma_\star(R)$ [M$_\odot$ kpc$^{-2}$]")
    ax.set_xlim(0, 20)
    ax.legend(fontsize=8)
    axr = ax.inset_axes([0.45, 0.55, 0.5, 0.38])
    axr.plot(d["R_c"], d["SigR_m"] / d["SigR_p"] - 1, color=ACCENT, lw=1)
    axr.axhline(0, color=INK_MUTED, lw=0.6)
    axr.set_ylim(-0.15, 0.15)
    axr.set_xlim(0, 20)
    axr.tick_params(labelsize=7)
    axr.set_ylabel("ratio − 1", fontsize=7)

    ax = fig.add_subplot(gs[1, 1])
    for k, (R0, (zc, rp, rm)) in enumerate(d["rho_z"].items()):
        col = BLUE_RAMP[2 + 2 * k]
        ax.semilogy(zc, rp, color=col, lw=2, label=f"R = {R0:g} kpc")
        ax.semilogy(zc, rm, color=ACCENT, lw=1.1, ls="--")
    ax.set_xlabel(r"$|z|$ [kpc]")
    ax.set_ylabel(r"$\rho_\star$ [M$_\odot$ kpc$^{-3}$] (azimuthal mean)")
    ax.legend(fontsize=8, title="particles; dashed = potential", title_fontsize=7.5)

    ax = fig.add_subplot(gs[1, 2])
    ax.plot(d["R_c"], d["A2_p"], color=INK_2, lw=2, label=r"$A_2$ particles")
    ax.plot(d["R_c"], d["A2_m"], color=ACCENT, lw=1.3, ls="--", label=r"$A_2$ potential")
    ax.set_xlim(0, 12)
    ax.set_xlabel(r"$R$ [kpc]")
    ax.set_ylabel(r"$A_2$ (projected)")
    ax2 = ax.twinx()
    ax2.plot(d["R_c"], np.degrees(d["ph2_p"]), color=BLUE_RAMP[4], lw=1.5, label=r"$\psi_2$ particles")
    ax2.plot(d["R_c"], np.degrees(d["ph2_m"]), color=ACCENT, lw=1, ls=":", label=r"$\psi_2$ potential")
    ax2.set_ylim(-90, 90)
    ax2.set_ylabel(r"$\psi_2$ [deg]", color=BLUE_RAMP[4])
    ax.legend(loc="upper right", fontsize=8)
    ax2.legend(loc="lower right", fontsize=8)
    fig.suptitle("Au18 stellar CylSpline (mmax = 8, 30×30 grid): density implied by the potential vs "
                 "star particles", y=0.96, fontsize=11)
    savefig(fig, "au18_potential_density_check")

    # ---------- 2. rotation curve and asymmetric drift ----------
    R_c = np.arange(0.5, 20.0, DR_KIN) + DR_KIN / 2
    vc = {k: vcirc(p, R_c) for k, p in comps.items()}
    vc["total"] = vcirc(total, R_c)
    jv = jeans_vphi(R_c, x, y, z, vx, vy, vz, m, vc["total"])

    hdr = au.read_header(HALO_DIR / "snapdir_127")
    with h5py.File(OUT_FILE, "r") as f:
        frame = {k: f.attrs[k] for k in ("centre_box_kpc", "vcentre_kms", "rotation")}
    gas = read_mwframe(0, ["Velocities", "Masses"], 25.0, frame, hdr)
    gp, gv, gm = gas["pos"], gas["vel"], gas["Masses"].astype(np.float64)
    gR = np.hypot(gp[:, 0], gp[:, 1])
    gvphi = -(gp[:, 0] * gv[:, 1] - gp[:, 1] * gv[:, 0]) / gR
    gidx = np.digitize(gR, np.append(R_c - DR_KIN / 2, R_c[-1] + DR_KIN / 2)) - 1
    gas_vphi = np.full(len(R_c), np.nan)
    for i in range(len(R_c)):
        s = (gidx == i) & (np.abs(gp[:, 2]) < Z_SLAB)
        if s.sum() >= 30:
            gas_vphi[i] = np.average(gvphi[s], weights=gm[s])

    outside = R_c > R_bar
    frac = (jv["vphi"] - jv["vphi_pred"]) / vc["total"]
    gas_frac = (gas_vphi - vc["total"]) / vc["total"]
    results["vc_R8"] = float(np.interp(8.0, R_c, vc["total"]))
    results["jeans_resid_outside_bar_median_rms"] = [float(np.nanmedian(frac[outside & (R_c < 15)])),
                                                    float(np.sqrt(np.nanmean(frac[outside & (R_c < 15)]**2)))]
    results["gas_resid_outside_bar_median"] = float(np.nanmedian(gas_frac[outside & (R_c < 15)]))
    results["asym_drift_R8_kms"] = float(np.interp(8.0, R_c, vc["total"] - jv["vphi"]))

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.4, 7.6), sharex=True,
                                 gridspec_kw=dict(height_ratios=[2.6, 1], hspace=0.06))
    a1.plot(R_c, vc["total"], color=INK, lw=2.2, label=r"$v_c$ total (potential)")
    for k, ls in (("dm", "--"), ("stars", "-."), ("gas", ":"), ("bh", (0, (1, 3)))):
        a1.plot(R_c, vc[k], color=INK_MUTED, lw=1.2, ls=ls, label=f"$v_c$ {k}")
    a1.plot(R_c, jv["vphi_pred"], color=ACCENT, lw=2, label=r"stellar $\langle v_\phi\rangle$ predicted (Jeans)")
    a1.plot(R_c, jv["vphi"], "o", color=BLUE_RAMP[4], ms=4.5, label=r"stellar $\langle v_\phi\rangle$ measured")
    a1.plot(R_c, gas_vphi, "s", color=BLUE_RAMP[1], mec=BLUE_RAMP[5], ms=4, label=r"gas $\langle v_\phi\rangle$")
    a1.axvspan(0, R_bar, color="#f0efec", lw=0, zorder=0)
    a1.text(0.3, 280, "bar region\n(axisymmetric Jeans\nnot valid)", fontsize=8, color=INK_2)
    a1.set_ylabel(r"velocity [km s$^{-1}$]")
    a1.set_ylim(0, 420)
    a1.legend(fontsize=7.5, loc="upper right", ncol=2, bbox_to_anchor=(1.0, 1.0))
    a1.set_title(rf"Au18 rotation curve and asymmetric drift ($|z| < {Z_SLAB}$ kpc)", color=INK_2)
    a2.axhline(0, color=INK_MUTED, lw=0.8)
    a2.plot(R_c, frac, "o-", color=ACCENT, ms=3.5, lw=1, label="stars: (measured − Jeans) / $v_c$")
    a2.plot(R_c, gas_frac, "s", color=BLUE_RAMP[2], mec=BLUE_RAMP[5], ms=3.5, label="gas: ($v_\\phi$ − $v_c$) / $v_c$")
    a2.axvspan(0, R_bar, color="#f0efec", lw=0, zorder=0)
    a2.set_ylim(-0.15, 0.15)
    a2.set_xlim(0, 20)
    a2.set_xlabel(r"$R$ [kpc]")
    a2.set_ylabel("fractional residual")
    a2.legend(fontsize=7.5, loc="lower right")
    savefig(fig, "au18_rotation_curve_check")

    # ---------- 3. extra: AGAMA vs Arepo potential at the star particles ----------
    st = au.read_particles(HALO_DIR / "snapdir_127", 4, ["Potential", "GFM_StellarFormationTime"],
                           frame["centre_box_kpc"], 30.0, hdr)
    keep = st["GFM_StellarFormationTime"] > 0
    pos = st["Coordinates"][keep].astype(np.float64) @ frame["rotation"].T
    phi_sim = st["Potential"][keep].astype(np.float64) / hdr["Time"]        # (km/s)^2
    phi_ag = total.potential(pos)
    r = np.linalg.norm(pos, axis=1)
    s = (r > 1) & (r < 20)
    slope, offset = np.polyfit(phi_sim[s], phi_ag[s], 1)
    resid = phi_ag - phi_sim - np.median((phi_ag - phi_sim)[s])
    results["phi_fit_slope"] = float(slope)
    results["phi_resid_rms_1_20kpc"] = float(np.std(resid[s]))
    results["phi_range_1_20kpc"] = float(np.ptp(np.percentile(phi_ag[s], [1, 99])))

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw=dict(wspace=0.35))
    r_e = np.logspace(-1, np.log10(30), 40)
    r_c = np.sqrt(r_e[1:] * r_e[:-1])
    ri = np.digitize(r, r_e) - 1
    q = np.array([np.percentile(resid[ri == i], [16, 50, 84]) if np.sum(ri == i) > 50 else [np.nan] * 3
                  for i in range(len(r_c))])
    a1.fill_between(r_c, q[:, 0], q[:, 2], color=BLUE_RAMP[0], lw=0, label="16–84%")
    a1.plot(r_c, q[:, 1], color=BLUE_RAMP[4], lw=2, label="median")
    a1.axhline(0, color=INK_MUTED, lw=0.8)
    a1.set_xscale("log")
    a1.set_xlabel(r"$r$ [kpc]")
    a1.set_ylabel(r"$\Phi_{\rm AGAMA} - \Phi_{\rm Arepo}$ − const [(km s$^{-1}$)$^2$]")
    a1.set_title(f"Residual at star particles (fit slope {slope:.4f})", color=INK_2)
    a1.legend(fontsize=8)
    e = np.arange(-MAP_LIM, MAP_LIM + 1e-9, 0.5)
    mid = np.abs(pos[:, 2]) < 1
    W = np.histogram2d(pos[mid, 0], pos[mid, 1], bins=[e, e])[0]
    Q = np.histogram2d(pos[mid, 0], pos[mid, 1], bins=[e, e], weights=resid[mid])[0]
    with np.errstate(invalid="ignore"):
        mp = np.where(W >= 20, Q / W, np.nan).T
    lim = np.nanpercentile(np.abs(mp), 98)
    im = a2.pcolormesh(e, e, mp, cmap=CMAP_DIV, vmin=-lim, vmax=lim, rasterized=True)
    a2.set_aspect("equal")
    a2.set_xlabel(r"$x$ [kpc]")
    a2.set_ylabel(r"$y$ [kpc]")
    a2.set_title(r"Face-on mean residual, $|z| < 1$ kpc", color=INK_2)
    cb = fig.colorbar(im, cax=a2.inset_axes([1.03, 0, 0.045, 1]), extend="both")
    cb.set_label(r"[(km s$^{-1}$)$^2$]", fontsize=9)
    savefig(fig, "au18_potential_vs_arepo")

    with open(RESULTS_FILE, "w") as f:
        json.dump(results, f, indent=2)
    print(json.dumps(results, indent=1))
    print(f"peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
