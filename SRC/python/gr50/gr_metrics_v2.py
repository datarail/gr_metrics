def GR_static_plus_death(data, d_ct='Deadcount', x_ct='Cellcount', d_0='Day0DeadCnt', x_0='Day0Cnt', x_0t='Ctrlcount', d_0t='Ctrl_Deadcount', t='Time'):
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
    t: str, column name for treatment time.

    '''
    Time = data[t]
    if any(Time > 14):
        print('Time transformed in days for normalization of GR_d')
        Time = Time / 24

    Dratio = (data[d_ct] - data[d_0]).apply(lambda x: max(x, 1)
                                            ) / (data[x_ct] - data[x_0])
    Dratio_ctrl = (data[d_0t] - data[d_0]
                   ).apply(lambda x: max(x, 1)) / (data[x_0t] - data[x_0])
    gr = np.log2(data[x_ct] / data[x_0])
    gr_ctrl = np.log2(data[x_0t] / data[x_0])
    GR_s = np.power(
        2, ((1 + Dratio) * gr / ((1 + Dratio_ctrl) * gr_ctrl))) - 1
    GR_d = np.power(2, (((Dratio_ctrl) * gr_ctrl - (Dratio) * gr) / Time)) - 1

    return GR_s, GR_d
