from .hipparcos import hip_JD
from .gaia import gaia_JD
import numpy as np


#__all__="bin_by_time",

def bin_by_time(iad, gap=0.05):
    """
    Group observations whose consecutive timestamps differ by less than `gap` days
    (0.01 d ≈ 14 min), i.e. one bin per transit, and return the inverse-variance
    weighted mean of each bin. No rejection of points.
    Extra arrays in *cols (e.g. cpsi, spsi, parf) are averaged with the same weights.
    """

    t=gaia_JD(iad)
    y=iad[-2]
    sig=iad[-1]
    t, y, sig = (np.asarray(a, float) for a in (t, y, sig))
    cols = iad[0:-2]
    order = np.argsort(t)
    t, y, sig = t[order], y[order], sig[order]
    cols = [c[order] for c in cols]

    edges = np.r_[0, np.flatnonzero(np.diff(t) > gap) + 1, len(t)]
    out = []
    for i0, i1 in zip(edges[:-1], edges[1:]):
        sl = slice(i0, i1)
        w = 1.0 / sig[sl]**2
        W = w.sum()
        mean = lambda a: np.sum(w * a[sl]) / W
        yb = mean(y)
        n = i1 - i0
        scale = np.sqrt(np.sum(w * (y[sl] - yb)**2) / (n - 1)) if n > 1 else 1.0
        sig_b = W**-0.5 * max(1.0, scale)   # inflated when the scatter exceeds the formal errors
        out.append((mean(t), yb, sig_b, n, *[mean(c) for c in cols]))

    t_b, y_b, sig_b, n_obs, *cols_b = map(np.array, zip(*out))

    
    return np.vstack([cols_b, y_b,sig_b])