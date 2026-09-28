"""Inverse problem: how ordered is the alloy, given a measured property?

    python python/7b_inverse/infer_alpha.py

Uses results/results.csv (no new simulations). The model itself is in
bayes_model.py; this script runs it and writes the tables and the figure.

Outputs
    results/inverse_variance.csv   Brown-Forsythe test and scatter model per property
    results/inverse_identify.csv   90 % interval width vs measurement noise
    results/inverse_loo.csv        leave-one-box-out check (coverage, error)
    results/inverse_examples.csv   posteriors for a few example measurements
    results/figures/fig8_inverse.png
"""

import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bayes_model import (S, PROPS, GRID, RANGE, ALPHA_MIN, ALPHA_MAX, DEG,
                         load, Forward, posterior, summary)

NOISE = [0.0, 0.005, 0.01, 0.02, 0.05]      # measurement error, fraction of value

# single properties and combinations tried as "measurements"
SETS = {
    "G_H": ["G_H"],
    "C66": ["C66"],
    "C44h": ["C44h"],
    "sy_GPa": ["sy_GPa"],
    "E_H": ["E_H"],
    "Ex_cij": ["Ex_cij"],
    "E_GPa": ["E_GPa"],
    "B_H": ["B_H"],
    "G_H+C44h+C66": ["G_H", "C44h", "C66"],
    "G_H+sy_GPa": ["G_H", "sy_GPa"],
}
LOO_SETS = ["G_H", "C66", "C44h", "sy_GPa", "E_GPa", "B_H", "G_H+C44h+C66"]

# worked examples: (property, measured value, measurement error)
EXAMPLES = [("G_H", 44.5, 0.01), ("G_H", 43.3, 0.01), ("G_H", 42.1, 0.01),
            ("G_H", 43.3, 0.02), ("G_H", 43.3, 0.05), ("E_GPa", 114.0, 0.01)]

FIG_DIR = S.RESULTS_DIR / "figures"
BLUE = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
CAT1, CAT2 = "#2a78d6", "#eb6834"

NAMES = list(PROPS)
COL = {p: NAMES.index(p) for p in NAMES}


def write(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"wrote {path.relative_to(S.ROOT)}")


def scatter_table(fwd):
    """0. which properties get an alpha1-dependent scatter"""
    rows = []
    for p in NAMES:
        j = COL[p]
        s0, s12 = fwd.sd(np.array([0.0, -0.12]))[:, j]
        model = "alpha1-dependent" if fwd.varies[j] else "constant"
        rows.append([p, f"{fwd.bf_p[j]:.4f}", model, f"{s0:.4f}", f"{s12:.4f}"])
    write(S.RESULTS_DIR / "inverse_variance.csv",
          ["property", "brown_forsythe_p", "scatter_model", "sd_at_alpha0", "sd_at_alpha-0.12"],
          rows)


def identifiability(fwd, a_true=-0.06):
    """1. how narrow is the 90 % interval for each property and noise level"""
    rows, width = [], {}
    for name, props in SETS.items():
        idx = [COL[p] for p in props]
        y = fwd.mean(a_true)[0, idx]         # noise-free "measurement" at a_true
        for nz in NOISE:
            m, sd, lo, hi = summary(posterior(fwd, idx, y, nz))
            width[(name, nz)] = hi - lo
            rows.append([name, nz, f"{m:.4f}", f"{sd:.4f}", f"{lo:.4f}", f"{hi:.4f}",
                         f"{hi - lo:.4f}", f"{(hi - lo) / RANGE:.3f}"])
    write(S.RESULTS_DIR / "inverse_identify.csv",
          ["properties", "meas_noise_frac", "post_mean", "post_sd", "lo90", "hi90",
           "width90", "width_over_range"], rows)
    return width


def leave_one_out(alpha, Y, case):
    """2. hide one box, infer its alpha1 from the other 79, repeat for all boxes"""
    variants = [(1, True, "linear"),
                (2, False, "quadratic, constant scatter"),
                (2, True, "quadratic")]
    rows, loo = [], {}
    for deg, varying, model in variants:
        for name in LOO_SETS:
            idx = [COL[p] for p in SETS[name]]
            est = []
            for k in range(len(alpha)):
                keep = np.arange(len(alpha)) != k
                f_k = Forward(alpha[keep], Y[keep], case[keep], deg, varying)
                m, sd, lo, hi = summary(posterior(f_k, idx, Y[k, idx], 0.0))
                est.append((alpha[k], m, lo, hi))
            est = np.array(est)
            inside = np.mean((est[:, 2] <= est[:, 0]) & (est[:, 0] <= est[:, 3]))
            err = est[:, 1] - est[:, 0]
            if deg == DEG and varying:
                loo[name] = est
            rows.append([model, name, len(alpha), f"{inside:.3f}",
                         f"{np.sqrt(np.mean(err ** 2)):.4f}", f"{np.mean(err):.4f}"])
    write(S.RESULTS_DIR / "inverse_loo.csv",
          ["model", "properties", "n_boxes", "coverage90", "rmse_alpha", "bias_alpha"], rows)
    return loo


def examples(fwd):
    """3. "measured G = ... GPa, how ordered is the alloy?" """
    rows, post = [], []
    for p, y, nz in EXAMPLES:
        pp = posterior(fwd, [COL[p]], y, nz)
        m, sd, lo, hi = summary(pp)
        post.append((y, pp))
        rows.append([p, y, nz, f"{m:.4f}", f"{sd:.4f}", f"{lo:.4f}", f"{hi:.4f}"])
    write(S.RESULTS_DIR / "inverse_examples.csv",
          ["property", "measured", "meas_noise_frac", "post_mean", "post_sd", "lo90", "hi90"],
          rows)
    return post


def make_figure(alpha, Y, case, fwd, post, width, loo):
    fig, ax = plt.subplots(2, 2, figsize=(10, 8))

    # (a) forward model for G with 90 % band
    a = ax[0, 0]
    j = COL["G_H"]
    for c, cc in zip(S.ALPHA_TARGET, BLUE):
        s = case == c
        a.plot(alpha[s], Y[s, j], "o", ms=4, color=cc, label=c.replace("ordered_", ""))
    mu = fwd.mean(GRID)[:, j]
    band = 1.645 * fwd.sd(GRID)[:, j] * np.sqrt(1 + fwd.fit_var(GRID))
    a.plot(GRID, mu, color="k", lw=1)
    a.fill_between(GRID, mu - band, mu + band, color="0.8", zorder=0)
    a.set_xlabel(r"$\alpha_1$")
    a.set_ylabel("G (Hill) [GPa]")
    a.set_title("(a) forward model, 90 % band")
    a.legend(frameon=False, fontsize=8)

    # (b) posteriors for three measured G values
    a = ax[0, 1]
    for (y, pp), cc in zip(post[:3], BLUE[1:]):
        a.plot(GRID, pp, color=cc, label=f"G = {y} GPa")
    a.set_xlabel(r"$\alpha_1$")
    a.set_ylabel("posterior density")
    a.set_title("(b) measured G (1 % noise)")
    a.legend(frameon=False, fontsize=8)

    # (c) interval width against measurement noise
    a = ax[1, 0]
    show = [("G_H", CAT1, "-"), ("C66", CAT1, "--"), ("sy_GPa", CAT2, "-"),
            ("E_GPa", CAT2, "--"), ("B_H", "0.5", ":")]
    x = np.array(NOISE) * 100
    for s, cc, ls in show:
        a.plot(x, [width[(s, nz)] for nz in NOISE], ls, marker="o", ms=4, color=cc,
               label=PROPS[s])
    a.axhline(0.9 * RANGE, color="k", lw=0.8, ls=":")
    a.text(0.1, 0.9 * RANGE + 0.002, "prior only", fontsize=8)
    a.set_xlabel("measurement noise [% of value]")
    a.set_ylabel(r"90 % interval width of $\alpha_1$")
    a.set_title(r"(c) which property pins down $\alpha_1$")
    a.legend(frameon=False, fontsize=8)

    # (d) leave-one-box-out check for G
    a = ax[1, 1]
    e = loo["G_H"]
    a.errorbar(e[:, 0], e[:, 1], yerr=[e[:, 1] - e[:, 2], e[:, 3] - e[:, 1]],
               fmt="o", ms=3, color=CAT1, ecolor="0.7", elinewidth=0.8)
    a.plot([ALPHA_MIN, ALPHA_MAX], [ALPHA_MIN, ALPHA_MAX], "k", lw=0.8)
    a.set_xlabel(r"true $\alpha_1$ of held-out box")
    a.set_ylabel(r"inferred $\alpha_1$ (mean, 90 %)")
    a.set_title("(d) leave-one-box-out, from G")

    fig.tight_layout()
    out = FIG_DIR / "fig8_inverse.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out.relative_to(S.ROOT)}")


def main():
    alpha, Y, case = load()
    fwd = Forward(alpha, Y, case)
    scatter_table(fwd)
    width = identifiability(fwd)
    loo = leave_one_out(alpha, Y, case)
    post = examples(fwd)
    make_figure(alpha, Y, case, fwd, post, width, loo)


if __name__ == "__main__":
    main()
