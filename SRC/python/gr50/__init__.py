import pandas as pd
import numpy as np
import scipy.optimize, scipy.stats


__all__ = ['gr_metrics', 'logistic']


def logistic(x, params):
    einf, mid, slope = params
    emin = 1.0
    return ( (emin-einf) / (1 + ((x/mid)**slope) ) ) + einf

def _logistic_inv(y, params):
    einf, mid, slope = params
    emin = 1.0
    if y >= min(emin, einf) and y <= max(emin, einf):
        return mid * ( (y-emin) / (einf-y) ) ** (1/slope)
    else:
        return np.inf

def _flat(x, params):
    y, = params
    return y

def _rss(params, fn, xdata, ydata):
    rss = 0.0
    for x, y in zip(xdata, ydata):
        rss += (y - fn(x, params)) ** 2
    return rss

def _tss(ydata):
    tss = 0.0
    y_mean = ydata.mean()
    for y in ydata:
        tss += (y - y_mean) ** 2
    return tss

def _rsquare(params, fn, xdata, ydata):
    ss_res = _rss(params, fn, xdata, ydata)
    ss_tot = _tss(ydata)
    return 1 - ss_res / ss_tot

def _fit(fn, xdata, ydata, bounds):
    res = scipy.optimize.minimize(_rss, x0=np.mean(bounds, 1),
                                  args=(fn, xdata, ydata), bounds=bounds)
    return res

def _calculate_pval(logistic_result, flat_result, n):
    rss2 = logistic_result.fun
    rss1 = flat_result.fun
    df1 = len(logistic_result.x) - len(flat_result.x)
    df2 = n - len(logistic_result.x) + 1
    f = ( (rss1-rss2)/df1 ) / (rss2/df2)
    pval = 1 - scipy.stats.f.cdf(f, df1, df2)
    return pval

def _mklist(values):
    """Convert tuple to list, and anything else to a list with just that thing.

    This is a helper to fix an inconsistency with the group keys in a
    pandas.groupby object. When grouping by multiple columns, the keys are
    tuples of values. When grouping by a single column, even if specified as a
    single-element list, the keys are single values. We always want a list,
    which is what this function accomplishes."""
    if isinstance(values, tuple):
        return list(values)
    else:
        return [values]

def _metrics(df, alpha):
    conc_min = df.concentration.min()
    conc_max = df.concentration.max()
    bounds = np.array([[-1, 1], [conc_min/10, conc_max*10], [0.1, 5]])
    logistic_result = _fit(logistic, df.concentration, df.gr, bounds)
    flat_result = _fit(_flat, df.concentration, df.gr, bounds[[0]])
    pval = _calculate_pval(logistic_result, flat_result, len(df.concentration))
    if pval > alpha or not logistic_result.success:
        # Return values for the metrics such that the logistic function will
        # produce the same curve as the flat fit.
        if flat_result.x[0] > 0.5:
            gr50 = np.inf
        else:
            gr50 = -np.inf
        inf = flat_result.x[0]
        ec50 = 0.0
        # Must be non-zero or the logistic function will error.
        slope = 0.01
        r2 = _rsquare(flat_result.x, _flat, df.concentration, df.gr)
    else:
        gr50 = _logistic_inv(0.5, logistic_result.x)
        inf = logistic_result.x[0]
        ec50 = logistic_result.x[1]
        slope = logistic_result.x[2]
        r2 = _rsquare(logistic_result.x, logistic, df.concentration, df.gr)
    # Take the minimum across the highest 2 doses to minimize the effect of
    # outliers (robust minimum).
    max_ = min(df.gr[-2:])
    auc = np.trapz(1 - df.gr, df.concentration)
    return [gr50, max_, inf, slope, ec50, r2, auc]

def gr_metrics(data, alpha=0.05):
    """Compute Growth Response metrics.

    The input dataframe must contain a column named 'concentration' with the
    dose values of the perturbing agent and a column named 'gr' with the
    corresponding growth response (GR) values. Columns named 'cell_count',
    'cell_count__ctrl' and 'cell_count__time0' are ignored.

    Multiple dose-response experiments may be combined into a single dataframe
    by adding extra 'key' columns to distinguish them from each other. Each
    unique combination of values in the key columns will produce one
    corresponding row in the returned dataframe. The columns are the keys
    (if present) and metric names.

    The computed metrics are:

    * gr50: Dose at which GR reaches 0.5.
    * gr_max: Maximum observed GR effect (minimum value).
    * gr_inf: Extrapolated GR value at infinite dose.
    * slope: Hill slope of fitted logistic curve.
    * ec50: Dose at which GR is halfway between 1 and gr_inf.
    * r2: R squared of fitted logistic curve.
    * gr_auc: Area under the curve of observed data points. (Mathematically
      this is calculated as 1-AUC so that increasing gr_auc values correspond
      to increasing agent effect)

    The input data for each experiment are fitted with a logistic curve. An
    F-test is then performed with the null hypothesis being that there is no
    dose response, i.e. the data can be fitted well with a straight horizontal
    line. If the null hypothesis is not rejected, the returned metrics are
    chosen such that the logistic curve determined by gr_inf, ec50 and slope is
    equivalent to the horizontal line fit, and gr50 is infinite (potentially
    positive *or* negative).

    Parameters
    ----------
    data : pandas.DataFrame
        Input data on which to compute the metrics.
    alpha : Optional[float]
        Significance level for the F-test.

    Returns
    -------
    pandas.DataFrame
        The computed metrics.

    Example
    -------
    >>> import pandas
    >>> data = pandas.DataFrame(
    ...            [['A', 0.001, 0.965], ['A', 0.01, 0.953], ['A', 0.1, 0.533],
    ...             ['A', 1.0, 0.0976], ['A', 10.0, 0.0188],
    ...             ['B', 0.001, 0.985], ['B', 0.01, 0.916], ['B', 0.1, 0.978],
    ...             ['B', 1.0, 1.04], ['B', 10.0, 0.936]],
    ...            columns=['drug', 'concentration', 'gr'])
    >>> print gr_metrics(data)
      drug      gr50  gr_max    gr_inf     slope      ec50            r2    gr_auc
    0    A  0.114026  0.0188  0.018108  1.145268  0.110412  9.985790e-01  9.115929
    1    B       inf  0.9360  0.971000  0.010000  0.000000 -1.176836e-14  0.105115
    """
    non_keys = set(('concentration', 'cell_count', 'cell_count__ctrl',
                    'cell_count__time0', 'gr'))
    metric_columns = ['gr50', 'gr_max', 'gr_inf', 'slope', 'ec50', 'r2',
                      'gr_auc']
    keys = list(set(data.columns) - non_keys)
    gb = data.groupby(keys)
    data = [_mklist(k) + _metrics(v, alpha) for k, v in gb]
    df = pd.DataFrame(data, columns=keys + metric_columns)
    return df
