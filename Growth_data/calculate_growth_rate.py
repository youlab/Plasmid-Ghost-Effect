# define the function to calculate the growth rate by taking the slope of the log-transformed OD600 data
# with a rolling window
import numpy as np

def GR(
    t: np.ndarray,
    y: np.ndarray,
    window_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    For each rolling window of `window_size` consecutive points, fit y ≈ a·t + b
    and record slope `a`. Assign that slope to the time at the middle index of the window.

    Parameters
    ----------
    t
        Time (or any x-axis), same length as y.
    y
        Observed values. For exponential growth rate on OD, pass ``np.log(OD)`` (with a floor).
    window_size
        Number of points in each window; must be ≥ 2. Odd sizes place the center on an
        actual sample; even sizes use the mean value of the window as the center.

    Returns
    -------
    t_center
        Time at the center of each window, length ``len(t) - window_size + 1``.
    slope
        Least-squares slope in each window (growth rate in y-units per t-unit).
    """

    w = int(window_size)
    if w < 2:
        raise ValueError("window_size must be at least 2")
    n = t.size
    if n < w:
        raise ValueError(f"Need at least {w} points; got {n}")

    n_out = n - w + 1
    slopes = np.empty(n_out)
    t_center = np.empty(n_out)

    for k in range(n_out):
        sl = slice(k, k + w)
        tw, yw = t[sl], y[sl]
        slopes[k] = np.polyfit(tw, yw, 1)[0]
        t_center[k] = (tw[0]+tw[1])/2

    return t_center, slopes
