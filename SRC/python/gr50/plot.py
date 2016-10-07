import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from . import logistic

def plot_curves(gr_metrics, gr_values, colorvar=None, colvar=None, rowvar=None):
    '''Produce a trellis plot showing the fitted curves and some of the metrics
     across the different cell lines and drugs.'''

    # padding with dummy columns to have general variables
    grm = gr_metrics.copy()
    grv = gr_values.copy()
    if rowvar is None:
        grm['_r'] = 0
        grv['_r'] = 0
        rowvar = '_r'
    if colvar is None:
        grm['_c'] = 0
        grv['_c'] = 0
        colvar = '_c'
    if colorvar is None:
        grm['_cl'] = 0
        grv['_cl'] = 0
        colorvar = '_cl'

    # display parameters
    colorkeys = grv[colorvar].drop_duplicates().values
    cmap = matplotlib.cm.get_cmap('jet')
    colors = cmap(np.array(range(len(colorkeys)))/(1.*len(colorkeys)))
    colors = np.maximum(colors**1.3-.1,0)
    x_min = grv.concentration.min() / 10
    x_max = grv.concentration.max() * 10
    fit_x = np.logspace(np.log10(x_min), np.log10(x_max))

    #### construction of the plot grid
    sns.set(style="ticks")
    # two cases based on rowvar
    if rowvar is '_r':
        grid = sns.FacetGrid(grv, col_wrap=np.ceil(np.sqrt(len(grv[colvar].drop_duplicates()))),
                             col=colvar, margin_titles=True)
    else:
        grid = sns.FacetGrid(grv, row=rowvar, col=colvar, margin_titles=True)
    # format the axes
    grid.set(xscale="log")
    grid.set_axis_labels(x_var='Concentration (uM)', y_var='GR value')
    grid.set(ylim=(-1, 1.15))
    grid.set(xlim=(x_min, x_max*3))
    grid.fig.tight_layout(w_pad=1)
    # set the positions of the axis labels
    if rowvar is not '_r':
        if colvar is '_c':
            grid.set_titles(row_template='{row_var}={row_name}', col_template='')
        else:
            grid.set_titles(template='{row_var}={row_name}, {col_var}={col_name}',
                            row_template='{row_var}={row_name}',
                            col_template='{col_var}={col_name}')
    else:
        if colvar is not '_c':
            grid.set_titles(template='{col_var}={col_name}',row_template='')
    # small hack to get uniform structure if rowvar was not specified
    if rowvar is '_r':
        grid.row_names = [0]
        grid.axes[0] = grid.axes.copy()

    #### loop to plot the data
    for rowval, row_axes in zip(grid.row_names, grid.axes):
        for colval, ax in zip(grid.col_names, row_axes):
            ax.hlines(0, x_min, x_max, '#707070', lw=0.5)
            h = []
            for i,colorval in enumerate(colorkeys):
                htemp, = ax.plot(np.nan, np.nan, 'o-', ms=5, color=colors[i,], lw=1, label=colorval)
                h.append(htemp)
                v = grv[(grv[colvar] == colval) & (grv[colorvar] == colorval) &
                        (grv[rowvar] == rowval)]
                ax.plot(x_max*np.array([1.5, 3]), [v.GRvalue.min()]*2, '-', color=colors[i,], lw=2)
                ax.plot(v.concentration, v.GRvalue, marker='o', ms=5, color=colors[i,], ls='None')
                for m in grm[(grm[colvar] == colval) & (grm[colorvar] == colorval) &
                                (grm[rowvar] == rowval)].itertuples():
                    fit_y = logistic(fit_x, [m.GRinf, np.log10(m.GEC50), m.h_GR])
                    ax.plot([m.GR50]*2, [-1, -.9], '-', color=colors[i,], lw=2)
                    ax.plot(fit_x, fit_y, lw=1, color=colors[i,])

    # legend (needs to be better positioned)
    plt.legend(tuple(h), tuple(colorkeys), loc='lower left')

    return plt.gcf()
