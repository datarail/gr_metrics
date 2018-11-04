import numpy as np

def compute_gr_static_death(data, d_ct='Deadcount', x_ct='Cellcount', d_0='Day0DeadCnt', x_0='Day0Cnt',
                            x_0t='Ctrlcount', d_0t='Ctrl_Deadcount', t='Time', t_unit='hour'):
    '''
    Calculate GR static and GR dead parameters from DDD experiment outputs.
    Parameters
    ========
    data: pandas.DataFrame, table of DDD experiment output. Should be a fixed format but just in case of future changes, critical columns are parameterized.
    d_ct: str, column name for treated dead cell count at concentration c.
    x_ct: str, column name for treated live cell count at concentration c.
    d_0: str, column name for t0 control dead cell count.
    x_0: str, column name for t0 control live cell count.
    d_0t: str, column name for control dead cell count.
    x_0t: str, column name for control live cell count.
    t_: str, column name for treatment time.

    '''
    Time = data[t]
    if t_unit == 'day':
        print('Time transformed in days for normalization of GR_d')
        Time = Time / 24

    dead_ratio = np.maximum(data[d_ct] - data[d_0], 1) / (data[x_ct] - data[x_0])
    dead_ratio_ctrl = np.maximum(data[d_0t] - data[d_0], 1) / (data[x_0t] - data[x_0])
    gr = np.log2(data[x_ct] / data[x_0])
    gr_ctrl = np.log2(data[x_0t] / data[x_0])

    def exp_2_func(x):
        return np.power(2,x)-1
    GR_s = exp_2_func((1 + dead_ratio) * gr / ((1 + dead_ratio_ctrl) * gr_ctrl))
    GR_d = exp_2_func(2, (((dead_ratio_ctrl) * gr_ctrl - (dead_ratio) * gr) / Time))
    processed_data = data.copy()
    processed_data['GR_s'] = GR_s
    processed_data['GR_d'] = GR_d
    return processed_data
