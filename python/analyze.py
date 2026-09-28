"""Step 7: statistics and uncertainty quantification.

    python python/7_analyze/analyze.py

Reads results/results.csv, tension.csv, cij.csv (step 6) and writes to results/:
    summary.csv      per case and property: mean, SD, CV, bootstrap 95 % CI of
                     mean and SD, Shapiro-Wilk and Anderson-Darling normality
    fits.csv         normal / lognormal / Weibull fits (AIC) for E and sy
    trend.csv        linear trend of each property with alpha1 (80 boxes)
    tests.csv        each ordered case vs disordered: Welch t, Cohen d,
                     Brown-Forsythe (spread); plus ANOVA / BF over all 4 cases
    size_check.csv   large (82 A) vs main (41 A) boxes
    acar_check.csv   disordered case vs Billah & Acar (2024)
    noise_check.csv  tension E vs E from the 0 K elastic constants, box by box
                     (separates configurational scatter from MD noise)
    uq_poly.csv      uncertainty propagated to a polycrystal (E along the
                     loading direction): Monte Carlo and a simplified
                     analytical estimate (no normalization constraint), with the split between crystal and texture

Polycrystal model
    The single-crystal property is the directional Young's modulus E(theta),
    theta = angle between the loading axis and the c-axis, from each box's
    elastic constants. The texture is K = 50 orientation nodes with weights
    (volume fractions) A_i; the polycrystal value is E_poly = sum A_i E(theta_i),
    the same linear homogenization as Billah & Acar. Two textures are used:
    random, and basal (c-axes ~ perpendicular to the load, 20 deg spread).
    Texture uncertainty: each A_i has a 5 % relative SD (normal), then the
    weights are renormalized to sum to 1. Crystal uncertainty: the elastic
    constants are drawn from a multivariate normal fitted to the 20 boxes.
"""

import csv
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))

import settings as S

RNG = np.random.default_rng(20260924)
N_BOOT = 10000
N_MC = 100000
K_NODES = 50
TEX_SD = 0.05
BASAL_SPREAD_DEG = 20.0

PROPS = ["E_GPa", "sy_GPa", "nu_xy", "nu_xz", "C11h", "C12", "C13h", "C33", "C44h", "C66",
         "B_H", "G_H", "Ex_cij", "a_A", "c_A", "pe_atom_eV"]
CASES = list(S.ALPHA_TARGET)

ACAR = {"E_GPa": (123.344, 2.447), "sy_GPa": (5.314, 0.0472)}   # 160 random boxes, 1e10 /s


# ---------------------------------------------------------------- io

def read_csv(name):
    with open(S.RESULTS_DIR / name, newline="") as f:
        return list(csv.DictReader(f))


def write_csv(name, rows):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(S.RESULTS_DIR / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, (float, np.floating)) else v)
                        for k, v in r.items()})


def column(rows, key, case=None):
    return np.array([float(r[key]) for r in rows
                     if (case is None or r["case"] == case) and r.get(key, "") != ""])


# ---------------------------------------------------------------- statistics

def boot_ci(x, fn, n=N_BOOT):
    idx = RNG.integers(0, len(x), size=(n, len(x)))
    vals = fn(x[idx], axis=1) if fn is not np.std else np.std(x[idx], axis=1, ddof=1)
    return np.percentile(vals, [2.5, 97.5])


def cohen_d(a, b):
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / sp


def aic_fits(x):
    out = {}
    for name, dist, kw in [("normal", stats.norm, {}),
                           ("lognormal", stats.lognorm, {"floc": 0}),
                           ("weibull", stats.weibull_min, {"floc": 0})]:
        p = dist.fit(x, **kw)
        k = len(p) - len(kw)
        out[name] = 2 * k - 2 * np.sum(dist.logpdf(x, *p))
    return out


# ---------------------------------------------------------------- polycrystal model

VOIGT_PAIRS = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


def compliance_tensor(Sv):
    """6x6 Voigt compliance -> 3x3x3x3 tensor."""
    T = np.zeros((3, 3, 3, 3))
    for I, (i, j) in enumerate(VOIGT_PAIRS):
        for J, (k, l) in enumerate(VOIGT_PAIRS):
            f = (1 if I < 3 else 0.5) * (1 if J < 3 else 0.5)
            for a, b in {(i, j), (j, i)}:
                for c, d in {(k, l), (l, k)}:
                    T[a, b, c, d] = Sv[I, J] * f
    return T


def directional_E(Cv, thetas):
    """Young's modulus along n = (sin t, 0, cos t) in the crystal frame (z = c-axis)."""
    St = compliance_tensor(np.linalg.inv(Cv))
    n = np.stack([np.sin(thetas), np.zeros_like(thetas), np.cos(thetas)], axis=1)
    inv = np.einsum("pi,pj,pk,pl,ijkl->p", n, n, n, n, St)
    return 1.0 / inv


def cij_matrix(r):
    C = np.zeros((6, 6))
    for key, (i, j) in {"C11": (0, 0), "C22": (1, 1), "C33": (2, 2), "C12": (0, 1),
                        "C13": (0, 2), "C23": (1, 2), "C44": (3, 3), "C55": (4, 4),
                        "C66": (5, 5)}.items():
        C[i, j] = C[j, i] = float(r[key])
    return C


def texture_weights(kind, thetas):
    if kind == "random":
        w = np.ones_like(thetas)            # nodes are equal steps in cos(theta)
    else:                                   # basal: c-axes near 90 deg from the load
        w = np.exp(-0.5 * ((np.degrees(thetas) - 90.0) / BASAL_SPREAD_DEG) ** 2)
    return w / w.sum()


def propagate(C_boxes, thetas, mu_A):
    """Monte Carlo (total, crystal only, texture only) and the analytical formula."""
    keys = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2), (3, 3), (4, 4), (5, 5)]
    X = np.array([[C[i, j] for i, j in keys] for C in C_boxes])
    mu, cov = X.mean(0), np.cov(X, rowvar=False)

    def build(x):
        C = np.zeros((6, 6))
        for v, (i, j) in zip(x, keys):
            C[i, j] = C[j, i] = v
        return C

    n_cr = 4000                                             # crystal draws
    E_cr = np.array([directional_E(build(x), thetas)
                     for x in RNG.multivariate_normal(mu, cov, size=n_cr)])
    E_mean = directional_E(build(mu), thetas)

    A = mu_A * (1 + TEX_SD * RNG.standard_normal((N_MC, len(mu_A))))
    A = np.clip(A, 0, None)
    A /= A.sum(1, keepdims=True)

    pick = RNG.integers(0, n_cr, size=N_MC)
    total = np.einsum("ij,ij->i", E_cr[pick], A)
    crystal_only = E_cr @ mu_A
    texture_only = A @ E_mean

    # simplified first-order analytical estimate: nodes treated as independent and
    # no variance constraint for the normalization of the weights (Billah & Acar
    # add that constraint, eq. 9); shows how much the constraint matters
    E_box = np.array([directional_E(C, thetas) for C in C_boxes])
    muP, sdP = E_box.mean(0), E_box.std(0, ddof=1)
    sdA = TEX_SD * mu_A
    mean_an = np.sum(muP * mu_A)
    var_an = np.sum((sdP * mu_A) ** 2 + (muP * sdA) ** 2 + (sdP * sdA) ** 2)

    return {"mc_mean": total.mean(), "mc_sd": total.std(),
            "mc_ci_lo": np.percentile(total, 2.5), "mc_ci_hi": np.percentile(total, 97.5),
            "sd_crystal_only": crystal_only.std(), "sd_texture_only": texture_only.std(),
            "an_mean": mean_an, "an_sd": np.sqrt(var_an),
            "an_sd_crystal_term": np.sqrt(np.sum((sdP * mu_A) ** 2)),
            "an_sd_texture_term": np.sqrt(np.sum((muP * sdA) ** 2))}


# ---------------------------------------------------------------- main

def main():
    res = read_csv("results.csv")
    ten = read_csv("tension.csv")
    cij = read_csv("cij.csv")

    # 1. summary per case
    summ = []
    for p in PROPS:
        for c in CASES:
            x = column(res, p, c)
            if len(x) < 3:
                continue
            sw = stats.shapiro(x)
            ad = stats.anderson(x, "norm")
            ci_m, ci_s = boot_ci(x, np.mean), boot_ci(x, np.std)
            summ.append({"property": p, "case": c, "alpha": S.ALPHA_TARGET[c], "n": len(x),
                         "mean": x.mean(), "sd": x.std(ddof=1), "cv_pct": 100 * x.std(ddof=1) / abs(x.mean()),
                         "mean_ci_lo": ci_m[0], "mean_ci_hi": ci_m[1],
                         "sd_ci_lo": ci_s[0], "sd_ci_hi": ci_s[1],
                         "shapiro_p": sw.pvalue, "ad_stat": ad.statistic,
                         "ad_crit_5pct": ad.critical_values[2],
                         "normal_ok": bool(sw.pvalue > 0.05 and ad.statistic < ad.critical_values[2])})
    write_csv("summary.csv", summ)

    # 2. distribution fits for the two tension properties
    fits = []
    for p in ["E_GPa", "sy_GPa"]:
        for c in CASES:
            a = aic_fits(column(res, p, c))
            best = min(a, key=a.get)
            fits.append({"property": p, "case": c, **{f"aic_{k}": v for k, v in a.items()}, "best": best})
    write_csv("fits.csv", fits)

    # 3. trend with alpha1 over all 80 boxes
    trend = []
    alpha = column(res, "alpha1")
    for p in PROPS:
        y = column(res, p)
        if len(y) != len(alpha):
            continue
        lr = stats.linregress(alpha, y)
        base = column(res, p, "disordered").mean()
        trend.append({"property": p, "slope_per_unit_alpha": lr.slope, "slope_se": lr.stderr,
                      "change_per_-0.1_alpha": -0.1 * lr.slope,
                      "change_per_-0.1_alpha_pct": -10 * lr.slope / base,
                      "r2": lr.rvalue ** 2, "p": lr.pvalue})
    write_csv("trend.csv", trend)

    # 4. tests against the disordered case and across all four
    tests = []
    for p in PROPS:
        ref = column(res, p, "disordered")
        groups = [column(res, p, c) for c in CASES]
        anova = stats.f_oneway(*groups)
        bf_all = stats.levene(*groups, center="median")
        for c in CASES[1:]:
            x = column(res, p, c)
            t = stats.ttest_ind(x, ref, equal_var=False)
            bf = stats.levene(x, ref, center="median")
            tests.append({"property": p, "case": c, "mean_diff": x.mean() - ref.mean(),
                          "mean_diff_pct": 100 * (x.mean() - ref.mean()) / ref.mean(),
                          "welch_p": t.pvalue, "cohen_d": cohen_d(x, ref),
                          "sd_ratio": x.std(ddof=1) / ref.std(ddof=1), "bf_p": bf.pvalue,
                          "anova_p_all4": anova.pvalue, "bf_p_all4": bf_all.pvalue})
    write_csv("tests.csv", tests)

    # 5. size check
    size = []
    big = {("E_GPa", "B5b_tension80"): [], ("sy_GPa", "B5b_tension80"): []}
    for c in CASES:
        for p, rows, batch in [("E_GPa", ten, "B5b_tension80"), ("sy_GPa", ten, "B5b_tension80"),
                               ("C11h", cij, "B5c_cij80"), ("C33", cij, "B5c_cij80"),
                               ("C44h", cij, "B5c_cij80"), ("C66", cij, "B5c_cij80")]:
            xl = np.array([float(r[p]) for r in rows if r["case"] == c and r["batch"] == batch])
            xm = column(res, p, c)
            if len(xl) == 0:
                continue
            size.append({"case": c, "property": p, "large_mean": xl.mean(), "n_large": len(xl),
                         "main_mean": xm.mean(), "main_sd": xm.std(ddof=1),
                         "diff_pct": 100 * (xl.mean() - xm.mean()) / xm.mean(),
                         "diff_in_main_sd": (xl.mean() - xm.mean()) / xm.std(ddof=1)})
    write_csv("size_check.csv", size)

    # 6. comparison with Billah & Acar (disordered case)
    ac = []
    for p, (m, s) in ACAR.items():
        x = column(res, p, "disordered")
        ac.append({"property": p, "acar_mean": m, "acar_sd": s, "ours_mean": x.mean(),
                   "ours_sd": x.std(ddof=1), "diff_pct": 100 * (x.mean() - m) / m,
                   "sd_ratio": x.std(ddof=1) / s})
    write_csv("acar_check.csv", ac)

    # 6b. configurational scatter vs MD noise in E
    nc = []
    for c in CASES:
        E, Ex = column(res, "E_GPa", c), column(res, "Ex_cij", c)
        r, pv = stats.pearsonr(E, Ex)
        nc.append({"case": c, "sd_E_tension": E.std(ddof=1), "sd_Ex_cij": Ex.std(ddof=1),
                   "corr_r": r, "corr_p": pv,
                   "noise_sd_estimate": np.sqrt(max(E.var(ddof=1) - Ex.var(ddof=1), 0))})
    write_csv("noise_check.csv", nc)

    # 7. propagation to a polycrystal
    thetas = np.arccos((np.arange(K_NODES) + 0.5) / K_NODES)
    uq = []
    for c in CASES:
        C_boxes = [cij_matrix(r) for r in cij if r["case"] == c and r["batch"] == "B2_cij"]
        for tex in ("random", "basal"):
            mu_A = texture_weights(tex, thetas)
            out = propagate(C_boxes, thetas, mu_A)
            v_c, v_t, v = out["sd_crystal_only"] ** 2, out["sd_texture_only"] ** 2, out["mc_sd"] ** 2
            uq.append({"case": c, "alpha": S.ALPHA_TARGET[c], "texture": tex, **out,
                       "share_crystal_pct": 100 * v_c / v, "share_texture_pct": 100 * v_t / v,
                       "share_interaction_pct": 100 * (v - v_c - v_t) / v})
    write_csv("uq_poly.csv", uq)

    # ------------------------------------------------------------ printout
    def row(p, c):
        return next(r for r in summ if r["property"] == p and r["case"] == c)

    print("single crystal (20 boxes per case): mean, SD [95% CI of SD]")
    for p in ["E_GPa", "sy_GPa", "C44h", "C66", "G_H"]:
        print(f"  {p}")
        for c in CASES:
            r = row(p, c)
            print(f"    {c:<12} {r['mean']:8.3f}  SD {r['sd']:.3f} [{r['sd_ci_lo']:.3f}, {r['sd_ci_hi']:.3f}]"
                  f"  CV {r['cv_pct']:.2f}%  normal: {r['normal_ok']}")
    print("\ntrend with alpha1 (change for 0.1 more ordering)")
    for t in trend:
        if t["property"] in ("E_GPa", "sy_GPa", "C44h", "C66", "G_H", "B_H"):
            print(f"  {t['property']:<8} {t['change_per_-0.1_alpha']:+.3f} ({t['change_per_-0.1_alpha_pct']:+.2f} %)  "
                  f"R2 {t['r2']:.2f}  p {t['p']:.1e}")
    print("\nordered_a12 vs disordered")
    for t in tests:
        if t["case"] == "ordered_a12" and t["property"] in ("E_GPa", "sy_GPa", "G_H"):
            print(f"  {t['property']:<8} diff {t['mean_diff_pct']:+.2f} %  Welch p {t['welch_p']:.1e}  "
                  f"d {t['cohen_d']:.1f}  SD ratio {t['sd_ratio']:.2f}  BF p {t['bf_p']:.2f}")
    print("\npolycrystal E (MPa-level numbers in GPa)")
    for u in uq:
        print(f"  {u['case']:<12} {u['texture']:<7} mean {u['mc_mean']:.2f}  SD {u['mc_sd']:.3f} "
              f"(analytical {u['an_sd']:.3f})  crystal {u['share_crystal_pct']:.0f} % / texture {u['share_texture_pct']:.0f} %")


if __name__ == "__main__":
    main()
