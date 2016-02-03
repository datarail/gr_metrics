from __future__ import division
import math

# Computing the full set of metrics requires several "big" packages, but we
# still want the basic GR computation (which only uses the stdlib) to be
# available regardless.
try:
    import pandas as pd
    import numpy as np
    import scipy.optimize, scipy.stats
    _packages_available = True
except ImportError:
    _packages_available = False


__all__ = ['compute_gr', 'compute_gr_single', 'gr_metrics', 'logistic']


def _normalize_log2(n, n_0_0):
    normalized = n / n_0_0
    return math.log(normalized, 2)

def compute_gr_single(record):
    """Compute Growth Response value for a single sample.

    The input is a namedtuple or pandas Series with at least the following
    numeric fields:

    * cell_count: Number of cells detected in this sample.
    * cell_count__time0: Number of cells in the time=0 control for this sample.
    * cell_count__ctrl: Number of cells in the no-perturbation control.

    Parameters
    ----------
    record : Union[namedtuple, pandas.Series]
        Input data on which to compute the GR value.

    Returns
    -------
    float
        The computed GR value.

    Example
    -------
    >>> from collections import namedtuple
    >>> Record = namedtuple('Record',
    ...              ['cell_count', 'cell_count__ctrl', 'cell_count__time0'])
    >>> rec = Record(cell_count=1710, cell_count__ctrl=1766.0,
    ...              cell_count__time0=492.8)
    >>> print compute_gr_single(rec)
    0.965305500206
    """
    cc_t0 = float(record.cell_count__time0)
    log2nn = _normalize_log2(float(record.cell_count), cc_t0)
    log2nn_ctrl = _normalize_log2(float(record.cell_count__ctrl), cc_t0)
    gr = 2 ** (log2nn / log2nn_ctrl) - 1
    return gr

def compute_gr(data):
    """Compute Growth Response value for an entire dataset.

    The input dataframe must contain at least the following numeric fields:

    * cell_count: Number of cells detected per sample.
    * cell_count__time0: Number of cells in the time=0 control for each sample.
    * cell_count__ctrl: Number of cells in the no-perturbation control.

    The input must not already contain a column named 'gr'.

    A new dataframe will be returned with the GR values stored in a new 'gr'
    column.

    Parameters
    ----------
    data : pandas.DataFrame
        Input data on which to compute the metrics.

    Returns
    -------
    pandas.DataFrame
        Shallow copy of input data with a 'gr' column appended.

    Example
    -------

    """
    if 'gr' in data:
        raise ValueError("Data already contains a 'gr' column; aborting")
    result = data.copy(deep=False)
    result['gr'] = data.apply(compute_gr_single, axis=1)
    return result

def logistic(x, params):
    """Evaluate logistic curve equation using log-transformed x_0.

    Parameters
    ----------
    x
        X-value at which to evaluate the logistic function.
    params : list
        * einf: maximum Y-value (effect)
        * log10_mid: log10-transformed X-value of midpoint
        * slope: Steepness of curve

    Returns
    -------
    float
        Y-value of logistic function.
    """
    einf, log10_mid, slope = params
    emin = 1.0
    mid = 10 ** log10_mid
    return ( (emin-einf) / (1 + ((x/mid)**slope) ) ) + einf

def _logistic_inv(y, params):
    einf, log10_mid, slope = params
    emin = 1.0
    mid = 10 ** log10_mid
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

def _fit(fn, xdata, ydata, prior, bounds):
    res = scipy.optimize.minimize(_rss, args=(fn, xdata, ydata),
                                  x0=prior, bounds=bounds)
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
    conc_min = df.concentration.min() / 100
    conc_max = df.concentration.max() * 100
    bounds = np.array([[-1, 1], np.log10([conc_min, conc_max]), [0.1, 5]])
    prior = np.array([0.1, np.log10(np.median(df.concentration)), 2])
    logistic_result = _fit(logistic, df.concentration, df.gr, prior, bounds)
    flat_result = _fit(_flat, df.concentration, df.gr, prior[[0]], bounds[[0]])
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
        ec50 = 10 ** logistic_result.x[1]
        slope = logistic_result.x[2]
        r2 = _rsquare(logistic_result.x, logistic, df.concentration, df.gr)
    # Take the minimum across the highest 2 doses to minimize the effect of
    # outliers (robust minimum).
    max_ = min(df.gr[-2:])
    log_conc = np.log10(df.concentration)
    # Normalize AUC by concentration range (width of curve).
    auc_width = log_conc.max() - log_conc.min()
    auc = np.trapz(1 - df.gr, log_conc) / auc_width
    return [gr50, max_, auc, ec50, inf, slope, r2, pval]

def gr_metrics(data, alpha=0.05):
    """Compute Growth Response metrics for an entire dataset.

    The input dataframe must contain a column named 'concentration' with the
    dose values of the perturbing agent and a column named 'gr' with the
    corresponding growth response (GR) values. Columns named 'cell_count',
    'cell_count__ctrl' and 'cell_count__time0', which are used by the compute_gr
    function, will be ignored if they are still present in your dataframe.

    Multiple dose-response experiments may be combined into a single dataframe
    by adding extra 'key' columns to distinguish them from each other. Each
    unique combination of values in the key columns will produce one
    corresponding row in the returned dataframe. The columns are the keys
    (if present) and metric names.

    The computed metrics are:

    * GR50: Dose at which GR reaches 0.5.
    * GRmax: Maximum observed GR effect (minimum value).
    * GR_AUC: Area under the curve of observed data points. Mathematically this
      is calculated as 1-AUC so that increasing gr_auc values correspond to
      increasing agent effect. Also note the x-axis (concentration) values are
      log10-transformed, and the entire area is normalized by the width
      (difference between maximum and minimum concentration).
    * EC50: Dose at which GR is halfway between 1 and gr_inf.
    * GRinf: Extrapolated GR value at infinite dose.
    * Hill: Hill slope of fitted logistic curve.
    * r2: R squared of fitted logistic curve.
    * pval: P-value from the F-test (see below).

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
      drug      GR50   GRmax    GR_AUC      EC50     GRinf     Hill            r2  \\
    0    A  0.114027  0.0188  9.115929  0.110413  0.018109  1.14526  9.985790e-01   
    1    B       inf  0.9360  0.105115  0.000000  0.971000  0.01000 -1.176836e-14   
    <BLANKLINE>
           pval  
    0  0.000054  
    1  1.000000  
    """
    if not _packages_available:
        raise RuntimeError("Please install numpy, scipy and pandas in order "
                           "to use this function")
    non_keys = set(('concentration', 'cell_count', 'cell_count__ctrl',
                    'cell_count__time0', 'gr'))
    metric_columns = ['GR50', 'GRmax', 'GR_AUC', 'EC50', 'GRinf', 'Hill', 'r2',
                      'pval']
    keys = list(set(data.columns) - non_keys)
    gb = data.groupby(keys)
    data = [_mklist(k) + _metrics(v, alpha) for k, v in gb]
    df = pd.DataFrame(data, columns=keys + metric_columns)
    return df
