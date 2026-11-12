"""Minimal Cox PH (Breslow ties) + log-rank, numpy/scipy only."""
import numpy as np
from scipy import stats

def coxph(X, t, e, max_iter=50, tol=1e-9):
    X = np.asarray(X, float); t = np.asarray(t, float); e = np.asarray(e, float)
    n, p = X.shape
    mu, sd = X.mean(0), X.std(0); sd[sd == 0] = 1
    Z = (X - mu) / sd
    order = np.argsort(-t); Z, t, e = Z[order], t[order], e[order]  # descending time
    beta = np.zeros(p)
    for _ in range(max_iter):
        eta = Z @ beta; w = np.exp(eta)
        S0 = np.cumsum(w)                       # risk set sums (desc time => cumsum)
        S1 = np.cumsum(w[:, None] * Z, 0)
        S2 = np.cumsum(w[:, None, None] * Z[:, :, None] * Z[:, None, :], 0)
        # handle ties: use risk set at last index of each tie group
        ut, inv = np.unique(-t, return_inverse=True)
        last = np.zeros(len(ut), int)
        for i in range(n): last[inv[i]] = i
        S0g, S1g, S2g = S0[last][inv], S1[last][inv], S2[last][inv]
        Ebar = S1g / S0g[:, None]
        grad = (e[:, None] * (Z - Ebar)).sum(0)
        H = (e[:, None, None] * (S2g / S0g[:, None, None] - Ebar[:, :, None] * Ebar[:, None, :])).sum(0)
        step = np.linalg.solve(H, grad); beta = beta + step
        if np.abs(step).max() < tol: break
    cov = np.linalg.inv(H)
    b = beta / sd; se = np.sqrt(np.diag(cov)) / sd
    z = b / se; pval = 2 * stats.norm.sf(np.abs(z))
    return dict(coef=b, se=se, hr=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=pval)

def logrank(t, e, g):
    t = np.asarray(t, float); e = np.asarray(e, float); g = np.asarray(g)
    times = np.unique(t[e == 1]); O = E = V = 0.0
    for tt in times:
        at = t >= tt; n = at.sum(); n1 = (at & (g == 1)).sum()
        d = ((t == tt) & (e == 1)).sum(); d1 = ((t == tt) & (e == 1) & (g == 1)).sum()
        O += d1; E += d * n1 / n
        if n > 1: V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    chi = (O - E) ** 2 / V
    return stats.chi2.sf(chi, 1)
