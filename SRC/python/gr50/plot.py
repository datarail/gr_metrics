import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from . import logistic

def plot_curves(gr_metrics, df, rowvar, colvar):
    '''Produce a trellis plot showing the fitted curves and some of the metrics
     across the different cell lines and drugs.'''

    sns.set(style="ticks")
    grid = sns.FacetGrid(df, row=rowvar, col=colvar, margin_titles=True)
    grid.set(xscale="log")
    grid.map(plt.plot, "concentration", "GRvalue", lw=0, marker='o', ms=3)
    x_min = df.concentration.min() / 10
    x_max = df.concentration.max() * 10
    fit_x = np.logspace(np.log10(x_min), np.log10(x_max))
    for rowval, row_axes in zip(grid.row_names, grid.axes):
        for colval, ax in zip(grid.col_names, row_axes):
            for m in gr_metrics[(gr_metrics[colvar] == colval) &
                                (gr_metrics[rowvar] == rowval)].itertuples():
                fit_y = logistic(fit_x, [m.GRinf, np.log10(m.GEC50), m.h_GR])
                ax.hlines(0, x_min, x_max, '#707070', lw=0.5)
                ax.hlines(m.GRinf, x_min, x_max, '#ff00ff', linestyles='dashed',
                          lw=0.5)
                ax.vlines(m.GEC50, -1, 1, 'b', linestyles='dashed', lw=0.5)
                ax.vlines(m.GR50, -1, 1, 'k', lw=0.5)
                ax.plot(fit_x, fit_y, 'r', lw=1)
    grid.set(ylim=(-1, 1.1))
    grid.fig.tight_layout(w_pad=1)

    return plt.gcf()
