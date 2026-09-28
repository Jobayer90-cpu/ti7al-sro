"""Bayesian model for alpha1 (how ordered the alloy is) from measured properties.

Used by infer_alpha.py. The idea in four lines:

    forward model   property = f(alpha1)   quadratic fit to the 80 MD boxes
    scatter         spread of the boxes around f; constant, or linear in alpha1
                    when the Brown-Forsythe test says the spread changes
    prior           alpha1 equally likely anywhere in [ALPHA_MIN, ALPHA_MAX]
    posterior       prior x likelihood, evaluated on a fine grid of alpha1
"""

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.stats import levene

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import settings as S

# prior range: ALPHA_MIN = no Al-Al neighbours at 11.78 at% Al (-0.1336);
# ALPHA_MAX = +0.01 so that random boxes (alpha1 within +/-0.005) sit inside it
ALPHA_MIN = -S.X_AL / (1 - S.X_AL)
ALPHA_MAX = 0.01
RANGE = ALPHA_MAX - ALPHA_MIN
GRID = np.linspace(ALPHA_MIN, ALPHA_MAX, 1441)

LEVEL = 0.90        # credible level of the reported interval
DEG = 2             # forward model: 1 = straight line, 2 = quadratic
BF_ALPHA = 0.05     # Brown-Forsythe p below this -> scatter depends on alpha1

# properties that can be used as "measurements" (column name -> label)
PROPS = {
    "G_H":    "G (Hill)",
    "C66":    "C66",
    "C44h":   "C44",
    "sy_GPa": "yield (peak)",
    "E_H":    "E (Hill)",
    "Ex_cij": "E_x from Cij",
    "E_GPa":  "E from tension",
    "B_H":    "B (Hill)",
}


def load():
    """alpha1, property table (boxes x properties) and case name of every box."""
    with open(S.RESULTS_DIR / "results.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    alpha = np.array([float(r["alpha1"]) for r in rows])
    Y = np.array([[float(r[p]) for p in PROPS] for r in rows])
    case = np.array([r["case"] for r in rows])
    return alpha, Y, case


class Forward:
    """Forward model fitted to the boxes: mean curve and scatter of each property."""

    def __init__(self, alpha, Y, case, deg=DEG, varying=True):
        self.deg = deg

        # 1. least-squares fit of every property: y = c0 + c1*a (+ c2*a^2)
        X = self.design(alpha)
        self.XtXi = np.linalg.inv(X.T @ X)
        self.coef = self.XtXi @ X.T @ Y
        res = Y - X @ self.coef

        # 2. pooled scatter and correlation of the residuals (all 80 boxes)
        self.cov = res.T @ res / (len(alpha) - deg - 1)
        sd_pool = np.sqrt(np.diag(self.cov))
        self.corr = self.cov / np.outer(sd_pool, sd_pool)

        # 3. does the scatter change between ordering levels? (Brown-Forsythe)
        levels = list(dict.fromkeys(case))
        a_lev = np.array([alpha[case == c].mean() for c in levels])
        n = Y.shape[1]
        self.bf_p = np.empty(n)
        self.varies = np.zeros(n, bool)
        self.sd_coef = np.zeros((n, 2))          # scatter = s0 + s1 * alpha1
        for j in range(n):
            groups = [res[case == c, j] for c in levels]
            self.bf_p[j] = levene(*groups, center="median").pvalue
            if varying and self.bf_p[j] < BF_ALPHA:
                # straight line through the scatter of each level, scaled so
                # that its average size equals the pooled value
                sd_lev = np.array([np.sqrt(np.mean(g ** 2)) for g in groups])
                s1, s0 = np.polyfit(a_lev, sd_lev, 1)
                scale = sd_pool[j] / np.sqrt(np.mean((s0 + s1 * a_lev) ** 2))
                self.sd_coef[j] = (s0 * scale, s1 * scale)
                self.varies[j] = True
            else:
                self.sd_coef[j] = (sd_pool[j], 0.0)
        self.sd_floor = 0.25 * sd_pool           # keeps the line from reaching zero

    def design(self, a):
        """columns 1, a, a^2 ... for the polynomial fit"""
        return np.column_stack([a ** k for k in range(self.deg + 1)])

    def mean(self, a):
        """fitted property values at alpha1 = a (points x properties)"""
        return self.design(np.atleast_1d(a)) @ self.coef

    def sd(self, a):
        """box-to-box scatter at alpha1 = a (points x properties)"""
        a = np.atleast_1d(a)
        s = self.sd_coef[:, 0] + np.outer(a, self.sd_coef[:, 1])
        return np.maximum(s, self.sd_floor)

    def fit_var(self, a):
        """extra variance from the uncertainty of the fitted curve itself"""
        X = self.design(a)
        return np.einsum("ij,jk,ik->i", X, self.XtXi, X)


def posterior(fwd, idx, y, noise):
    """Posterior of alpha1 on GRID.

    fwd    fitted Forward model
    idx    column numbers of the measured properties
    y      measured values (same order as idx)
    noise  measurement error as a fraction of the value (0.01 = 1 %)
    """
    idx = list(idx)
    y = np.atleast_1d(y)
    mu = fwd.mean(GRID)[:, idx]              # predicted values at every alpha1
    R = fwd.corr[np.ix_(idx, idx)]           # correlation between properties
    D = fwd.sd(GRID)[:, idx]                 # scatter at every alpha1
    h = fwd.fit_var(GRID)
    Cm = np.diag((noise * y) ** 2)           # measurement error

    # Gaussian log-likelihood at every alpha1; the flat prior adds nothing
    logp = np.empty(len(GRID))
    for i in range(len(GRID)):
        V = R * np.outer(D[i], D[i]) * (1 + h[i]) + Cm
        r = y - mu[i]
        logdet = np.linalg.slogdet(V)[1]
        logp[i] = -0.5 * (r @ np.linalg.solve(V, r) + logdet)

    p = np.exp(logp - logp.max())
    return p / np.trapezoid(p, GRID)         # normalise to area 1


def summary(p):
    """Posterior mean, SD and 90 % highest-density interval (HDI).

    The HDI is used instead of equal tails because alpha1 = 0 sits near the
    edge of the prior: for a random box the posterior is cut off there, and
    an equal-tailed interval would miss 0 itself.
    """
    order = np.argsort(p)[::-1]
    mass = np.cumsum(p[order]) / p.sum()
    inside = order[: np.searchsorted(mass, LEVEL) + 1]
    lo, hi = GRID[inside].min(), GRID[inside].max()
    mean = np.trapezoid(GRID * p, GRID)
    sd = np.sqrt(np.trapezoid((GRID - mean) ** 2 * p, GRID))
    return mean, sd, lo, hi
