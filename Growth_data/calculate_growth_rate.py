import numpy as np
from sklearn.metrics import r2_score

def GR(
    t: np.ndarray,
    y: np.ndarray,
    window_size: int,
    min_r2: float = 0.9,
) -> tuple[np.ndarray, np.ndarray]:
    """
    For each rolling window of `window_size` consecutive points, fit y ≈ a·t + b
    and record slope `a` only if the linear fit has R² >= min_r2.

    Parameters
    ----------
    t
        Time values, same length as y.
    y
        Observed values. For exponential growth rate on OD, pass np.log(OD).
    window_size
        Number of points in each rolling window; must be >= 2.
    min_r2
        Minimum R² required to keep the slope.

    Returns
    -------
    t_center
        Time at the center of each valid window.
    slopes
        Least-squares slope in each window. Poor fits are discarded.
    """

    t = np.asarray(t)
    y = np.asarray(y)

    w = int(window_size)
    n = t.size

    if t.size != y.size:
        raise ValueError("t and y must have the same length")

    if w < 2:
        raise ValueError("window_size must be at least 2")

    if n < w:
        raise ValueError(f"Need at least {w} points; got {n}")

    n_out = n - w + 1

    slopes = []#np.full(n_out, np.nan)
    t_center = []#np.empty(n_out)
    #r2_values = np.full(n_out, np.nan)

    for k in range(n_out):
        sl = slice(k, k + w)
        tw, yw = t[sl], y[sl]

        slope, intercept = np.polyfit(tw, yw, 1)
        y_pred = slope * tw + intercept

        r2 = r2_score(yw, y_pred)

        if r2 >= min_r2:
            t_center.append(np.mean(tw))
            slopes.append(slope)
    
    t_center = np.array(t_center)
    slopes = np.array(slopes)

    return t_center, slopes