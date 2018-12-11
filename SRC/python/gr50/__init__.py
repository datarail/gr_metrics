from __future__ import division
import math
import logging

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


__all__ = ['compute_gr', 'compute_gr_single', 'gr_metrics', 'logistic',
           'compute_gr_static_toxic']
logger = logging.getLogger(__name__)


def _normalize_log2(n, n_0_0):
    normalized = max(n / n_0_0, 1e-6) # avoiding negative and null values
    return math.log(normalized, 2)

def compute_gr_single(record):
    """Compute Growth Response value for a single sample.

    The input is a namedtuple or pandas Series with at least the following
    numeric fields:

    * cell_count: Number of cells detected in this sample.
    * cell_count__time0: Number of cells in the treatment_duration=0 control for this sample.
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

def assign_ctrls(data, keys):
    '''Compute and assign the controls based on keys

    The input dataframe must contain the column cell_count
    The input keys is a list of column names that will be used to assign controls

    The column 'treatment_duration' is used as a default key and necessary
    Matching of conditions is based on the labels in the column 'role':
    conditions labeled 'treatmemt' are normalized to conditions 'negative_control'

    Parameters
    ----------
    data : pandas.DataFrame
        Input data on which to assign the control values.

    Returns
    -------
    pandas.DataFrame
        Deep copy of input data with columns 'cell_count__ctrl' and
        'cell_count__time0'  appended.

    Example
    -------
    df_out = assign_ctrl(df_in, ['cell_line', 'serum_pct'])

    '''


    # add 'treatment_duration' from the set as it is a default one
    keys = list(set(keys) - set(['treatment_duration']))

    df_keys = data[keys].drop_duplicates()
    df_keys.reset_index(inplace=True)
    dfout = pd.DataFrame(columns=data.columns).copy()

    for i in range(len(df_keys)):
        idx = np.ones(len(data))==1
        for k in keys:
            idx &= data[k] == df_keys[k][i]

        x0 = np.mean(data.loc[idx & (data.treatment_duration==0), 'cell_count'].values)
        times = data['treatment_duration'][idx & (data.treatment_duration>0)].drop_duplicates().values
        for t in times:
            idx_t = idx & (data.treatment_duration==t) & (data.role=='treatment')
            x_ctrl = np.mean(data.loc[idx & (data.treatment_duration==t) & \
                                      (data.role=='negative_control'), 'cell_count'].values)
            df_ctrl = pd.DataFrame(np.repeat([[x0, x_ctrl]],sum(idx_t),axis=0),
                                   columns=['cell_count__time0', 'cell_count__ctrl'])
            dfout = dfout.append(pd.concat([data.loc[idx_t, :].reset_index(drop=True), df_ctrl],
                                           axis=1),
                                 ignore_index=True)
    return dfout


def compute_gr(data):
    """Compute Growth Response value for an entire dataset.

    The input dataframe must contain at least the following numeric fields:

    * cell_count: Number of cells detected per sample.
    * cell_count__time0: Number of cells in the treatment_duration=0 control for each sample.
    * cell_count__ctrl: Number of cells in the no-perturbation control.

    The input must not already contain a column named 'GRvalue'.

    A new dataframe will be returned with the GR values stored in a new
    'GRvalue' column.

    Parameters
    ----------
    data : pandas.DataFrame
        Input data on which to compute the metrics.

    Returns
    -------
    pandas.DataFrame
        Shallow copy of input data with a 'GRvalue' column appended.

    Example
    -------

    """
    if 'GRvalue' in data:
        raise ValueError("Data already contains a 'GRvalue' column; aborting")
    result = data.copy(deep=False)
    result['GRvalue'] = data.apply(compute_gr_single, axis=1)
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

def _metrics(df, alpha, gr_value='GRvalue'):
    df = df.sort_values(by='concentration')
    conc_min = df.concentration.min() / 100
    conc_max = df.concentration.max() * 100
    bounds = np.array([[-1, 1], np.log10([conc_min, conc_max]), [0.1, 5]])
    prior = np.array([0.1, np.log10(np.median(df.concentration)), 2])
    logistic_result = _fit(logistic, df.concentration, df[gr_value],
                           prior, bounds)
    flat_result = _fit(_flat, df.concentration, df[gr_value],
                       prior[[0]], bounds[[0]])
    pval = _calculate_pval(logistic_result, flat_result, len(df.concentration))
    if pval > alpha or not logistic_result.success:
        # Return values for the metrics such that the logistic function will
        # produce the same curve as the flat fit.
        if flat_result.x[0] > 0.5:
            gr50 = np.inf
        else:
            gr50 = -np.inf
        inf = flat_result.x[0]
        gec50 = 0.0
        # Must be non-zero or the logistic function will error.
        slope = 0.01
        r2 = _rsquare(flat_result.x, _flat, df.concentration, df[gr_value])
    else:
        gr50 = _logistic_inv(0.5, logistic_result.x)
        inf = logistic_result.x[0]
        gec50 = 10 ** logistic_result.x[1]
        slope = logistic_result.x[2]
        r2 = _rsquare(logistic_result.x, logistic, df.concentration, df[gr_value])
    # Take the minimum across the highest 2 doses to minimize the effect of
    # outliers (robust minimum).
    max_ = min(df[gr_value][-2:])
    log_conc = np.log10(df.concentration)
    # Normalize AOC by concentration range (width of curve).
    aoc_width = log_conc.max() - log_conc.min()
    aoc = np.trapz(1 - df[gr_value], log_conc) / aoc_width
    return [gr50, max_, aoc, gec50, inf, slope, r2, pval]

def gr_metrics(data, alpha=0.05, gr_value='GRvalue',
               keys=['cell_line', 'agent', 'timepoint']):
    """Compute Growth Response metrics for an entire dataset.

    The input dataframe must contain a column named 'concentration' with the
    dose values of the perturbing agent and a column named 'GRvalue' with the
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
    * GR_AOC: Area over the curve of observed data points. Mathematically this
      is calculated as 1-GR so that increasing GR_AOC values correspond to
      increasing agent effect. Also note the x-axis (concentration) values are
      log10-transformed, and the entire area is normalized by the width
      (difference between maximum and minimum concentration).
    * GEC50: Dose at which GR is halfway between 1 and GRinf.
    * GRinf: Extrapolated GR value at infinite dose.
    * h_GR: Hill slope of fitted GR logistic curve.
    * r2_GR: R squared of fitted GR logistic curve.
    * pval_GR: P-value from the F-test for the GR logistic curve (see below).

    The input data for each experiment are fitted with a logistic curve. An
    F-test is then performed with the null hypothesis being that there is no
    dose response, i.e. the data can be fitted well with a straight horizontal
    line. If the null hypothesis is not rejected, the returned metrics are
    chosen such that the logistic curve determined by GRinf, GEC50 and h_GR is
    equivalent to the horizontal line fit, and GR50 is infinite (potentially
    positive *or* negative).

    Parameters
    ----------
    data : pandas.DataFrame
        Input data on which to compute the metrics.
    alpha : Optional[float]
        Significance level for the F-test.
    gr_value : Optional[str]
        Name of column containing GR values on which metrics are to be computed.
        Default is 'GRvalue'. Additional options include 'GR_static' and 'GR_toxic'.
    keys : list of str
        Name of columns by which to group the input dataframe.
        Default is ['cell_line', 'agent', 'timepoint'].

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
    ...            columns=['drug', 'concentration', 'GRvalue'])
    >>> print gr_metrics(data)
      drug      GR50   GRmax    GR_AOC     GEC50     GRinf      h_GR  \\
    0    A  0.114025  0.0188  0.481125  0.110411  0.018109  1.145262
    1    B       inf  0.9360  0.026375  0.000000  0.971000  0.010000
    <BLANKLINE>
              r2_GR   pval_GR
    0  9.985790e-01  0.000054
    1 -1.243450e-14  1.000000
    """
    if not _packages_available:
        raise RuntimeError("Please install numpy, scipy and pandas in order "
                           "to use this function")
    
    metric_columns = ['GR50', 'GRmax', 'GR_AOC', 'GEC50', 'GRinf', 'h_GR', 'r2_GR',
                      'pval_GR']
    
    gb = data.groupby(keys)
    data = [_mklist(k) + _metrics(v, alpha, gr_value) for k, v in gb]
    df = pd.DataFrame(data, columns=keys + metric_columns)
    return df


def compute_gr_static_toxic(data, time_col='timepoint'):
    """
    Computes gr_static and gr_toxic values for a given dataframe.

    The input dataframe must contain at least the following numeric fields:
    * cell_count: Total number of live cells detected per sample.
    * cell_count__time0: Total number of live cells in the treatment_duration=0 control for each sample.
    * cell_count__ctrl: Total number of cells in the no-perturbation control.
    * dead_count: Total number of dead cells detected per sample.
    * dead_count__time0: Total number of dead cells in the treatment_duration=0 control for each sample.
    * dead_count__ctrl: Total number of dead cells in the no-perturbation control.
    * timepoint : duration of drug treatment per column. The column name can be passed as an argument to
                 allow flexibility of passing time in any units (typically hours or days).
    * role : column that specifies if the sample is 'negative_control', 'positive_control', or 'treatment'.

    Parameters
    ----------
    data : pandas.DataFrame
       Input counts table on which to compute GR static and toxic values.
    time_col : Optional[str]
       Name of column in input dataframe containing duration of drug treatment. Default is 'timepoint'.

    Returns
    -------
    x : pandas.DataFrame
       GR_static and GR_toxic values appended to input dataframe
    
    """
    x = data.copy()
    # If total number of cells (live+dead) is less than 95% of time 0,
    # then increase estimate of dead cell count such that total number of cells is
    # equal to 95% of time 0 control.
    mc = ((x.dead_count + x.cell_count) <
          0.95 * (x.cell_count__time0 + x.dead_count__time0))  
    
    x.loc[mc, 'dead_count'] = np.maximum(x.loc[mc, 'dead_count'],
                                         (x.loc[mc, 'cell_count__time0'] +
                                          x.loc[mc, 'dead_count__time0'] -
                                          np.floor(.95 * x.loc[mc, 'cell_count']) +
                                          1)
                                         )
    logger.warning("%d wells or conditions have 5%% fewer cells than time0 control,"
                   " estimate of dead_count has been increased to compensate." % np.count_nonzero(mc))

    # If total number of cells (live+dead) is more than 115% of untreated control,
    # then reduce estimate of dead cell count such that total cells is equal to
    # 115% of untreated control.
    hd = ((x.role != 'negative_control') &
          ((x.dead_count + x.cell_count) >
           1.15 * (x.dead_count__ctrl + x.cell_count__ctrl))
          )
    
    x.loc[hd, 'dead_count'] = np. minimum(x.loc[hd, 'dead_count'],
                                          (x.loc[hd, 'cell_count__ctrl'] +
                                           x.loc[hd, 'dead_count__ctrl'] -
                                           np.ceil(1.15 * x.loc[hd, 'cell_count']) -
                                           1)
                                          )
    logger.warning("%d wells or conditions have too many cells relative to untreated control,"
                   " estimate of dead_count has been reduced to compensate." % np.count_nonzero(hd)
                   )
    
    d_ratio = np.maximum(x.dead_count - x.dead_count__time0, 1)/\
        (x.cell_count - x.cell_count__time0)
    d_ratio__ctrl = np.maximum(x.dead_count__ctrl - x.dead_count__time0, 1)/\
        (x.cell_count__ctrl - x.cell_count__time0)
    gr = np.log2(x.cell_count/x.cell_count__time0)
    gr__ctrl = np.log2(x.cell_count__ctrl/x.cell_count__time0)

    gr_static = 2 ** (
        ((1 + d_ratio) * gr)/
        ((1 + d_ratio__ctrl) * gr__ctrl)
        ) - 1

    gr_toxic = 2 ** (
        (d_ratio__ctrl * gr__ctrl - d_ratio * gr)/
         x[time_col]
        ) - 1

    x['GR_static'] = gr_static
    x['GR_toxic'] = gr_toxic

    # If number of live cells post treatment ~= time 0 control, analytical solution is not permissive.
    # Hence GR values are computed numerically using Taylor expansion
    fe = (np.abs(x.cell_count - x.cell_count__time0)/x.cell_count) < 1e-10
    if np.any(fe):
         logger.warning("%d wells or conditions have live cell counts approximately equal to time0 control,"
                        " therefore GR static and toxic values computed numerically using Taylor expansion." %
                        np.count_nonzero(fe)
                        )
         x_gr = pd.DataFrame(list(zip(gr, gr__ctrl , d_ratio__ctrl)),
                             columns=['gr', 'gr__ctrl', 'd_ratio__ctrl'])
         x_gr['d_ratio__gr'] = 0
         a = np.array(x.loc[fe, 'cell_count'])
         a = a.reshape(len(a), 1)
         b = np.array(x.loc[fe, 'cell_count__time0'])
         b = b.reshape(len(b), 1)
         nterms = 35 # Number of terms in Taylor expansion series
         
         counts_delta_rep_nterms = np.matlib.repmat((-1) * (a-b), 1, nterms)
         ts_exponents_rep_a =  np.matlib.repmat(range(1, nterms+1), len(a), 1)
         counts_time0_rep_nterms = np.matlib.repmat(b, 1, nterms)
         ts_exponents_rep_b =  np.matlib.repmat(range(1, nterms+1), len(b), 1)
         
         x_gr.loc[fe, 'd_ratio__gr'] = (
             np.maximum(x.loc[fe, 'dead_count'] - x.loc[fe, 'dead_count__time0'], 1).multiply(
                 (1/b + np.diag(
                     np.dot(
                         np.power(counts_delta_rep_nterms, ts_exponents_rep_a),
                         np.power(counts_time0_rep_nterms, ts_exponents_rep_b).T
                         )
                     ).reshape(len(b), 1)
                  ).flatten()
                 )
             )

         x.loc[fe, 'GR_static'] = 2 ** (
             (x_gr.loc[fe, 'gr'] + x_gr.loc[fe, 'd_ratio__gr'])/
             ((1 + x_gr.loc[fe, 'd_ratio__ctrl']) * x_gr.loc[fe, 'gr__ctrl'])
             ) - 1

         x.loc[fe, 'GR_toxic'] = 2 ** (
             ((x_gr.loc[fe, 'd_ratio__ctrl'] * x_gr.loc[fe, 'gr__ctrl'] - x_gr.loc[fe, 'd_ratio__gr'])/
              x.loc[fe, time_col])
             ) - 1      
    
    return x
