"""Statistical utilities for MD time-series analysis."""
import numpy as np
from scipy import stats


def effective_sample_size(data, max_lag=None):
    """Estimate effective sample size using Geyer's initial positive sequence estimator."""
    data = np.asarray(data, dtype=float)
    n = len(data)
    if n < 2:
        return 1.0

    # Normalized autocorrelation
    c0 = np.var(data, ddof=1)
    if c0 == 0:
        return float(n)

    max_lag = max_lag or min(n // 3, 1000)
    rho_sum = 0.0
    for t in range(1, max_lag + 1):
        c_t = np.mean((data[:-t] - np.mean(data)) * (data[t:] - np.mean(data)))
        rho_t = c_t / c0
        if rho_t < 0:
            break
        rho_sum += rho_t

    tau = 1.0 + 2.0 * rho_sum
    return max(1.0, n / tau)


def correlated_ttest(a, b):
    """Welch-style t-test for correlated time series (e.g. RMSD)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    diff = a - b
    n_eff = effective_sample_size(diff)
    t_stat = np.mean(diff) / (np.std(diff, ddof=1) / np.sqrt(n_eff))
    p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=n_eff - 1))
    return t_stat, p_value, np.mean(a), np.mean(b)
