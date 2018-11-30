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
    processed_data = data.copy()
    Time = processed_data[t]
    if t_unit == 'day':
        print('Time transformed in days for normalization of GR_d')
        Time = Time / 24
    # corrections for missing cells.
    # Test existance of missing cells.
    r_w_miss_cs = processed_data.Deadcount+processed_data.Cellcount < .95*(processed_data.Day0Cnt + processed_data.Day0DeadCnt)
    r_w_miss_cs = r_w_miss_cs[r_w_miss_cs.values].index
    num_missing_cells = len(r_w_miss_cs)
    if num_missing_cells>0:
        print("{} rows had more than 5% missing cells, which will be assigned to DeadCount".format(num_missing_cells))
        missing_cells_count = processed_data.Day0Cnt[r_w_miss_cs]+processed_data.Day0DeadCnt[r_w_miss_cs]-np.floor(0.95*processed_data.Cellcount[r_w_miss_cs])+1
        # new function in pandas for modifying processed_data inplace.
        missing_cells_count.name = 'Deadcount'
        processed_data.Deadcount.update(missing_cells_count)

    # corrections for too many dead cells.
    # Test existance of abnormally large amount of dead cells.
    r_w_dd_cs = (processed_data.pert_type!='ctl_vehicle') & (processed_data.Deadcount+processed_data.Cellcount > 1.15*(processed_data.Ctrlcount + processed_data.Ctrl_Deadcount))
    r_w_dd_cs = r_w_dd_cs[r_w_dd_cs.values].index
    num_dd_cs = len(r_w_dd_cs)

    if num_dd_cs>0:
        print("{} rows too many cells, which will corrected as negative deadcounts".format(num_dd_cs))
        # getting corrected dead counts, these values will be negative as it reflects the treatment actually increased cell numbers.
        # However I do feel this does not make since what is really needed is a thrid GR metric for stimulatory effect which promotes cell growth/division
        # Should not be a problem to implement. 
        updated_deadcount = processed_data.Ctrlcount[r_w_dd_cs] + processed_data.Ctrl_Deadcount[r_w_dd_cs] - np.ceil(1.15*processed_data.Cellcount[r_w_dd_cs])-1
        updated_deadcount.name='Deadcount'
        processed_data.Deadcount.update(updated_deadcount)

    dead_ratio = np.maximum(processed_data[d_ct] - processed_data[d_0], 1) / (processed_data[x_ct] - processed_data[x_0])
    dead_ratio_ctrl = np.maximum(processed_data[d_0t] - processed_data[d_0], 1) / (processed_data[x_0t] - processed_data[x_0])
    gr = np.log2(processed_data[x_ct] / processed_data[x_0])
    gr_ctrl = np.log2(processed_data[x_0t] / processed_data[x_0])

    def exp_2_func(x):
        return np.power(2,x)-1
    GR_s = exp_2_func((1 + dead_ratio) * gr / ((1 + dead_ratio_ctrl) * gr_ctrl))
    GR_d = exp_2_func(((dead_ratio_ctrl) * gr_ctrl - (dead_ratio) * gr) / Time)
    processed_data['GR_s'] = GR_s
    processed_data['GR_d'] = GR_d
    return processed_data