import matplotlib.pyplot as plt
import numpy as np
import scipy.stats


def spearman_plot(
        vals1,
        vals2,
        ax=None,
        xlabel=None,
        ylabel=None,
        title=None,
        fig_fpath=None,
        show=False,
        ):
    """Plot spearman ranking plot, with spearman r and pvalue in title."""
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    if title is None:
        r_spear, pval = scipy.stats.spearmanr(vals1, vals2)
        title = 'Spearman r: %f, P-value: %g' % (r_spear, pval)

    ranks1 = scipy.stats.rankdata(vals1)
    ranks2 = scipy.stats.rankdata(vals2)

    ax.plot(ranks1, ranks2, 'o', alpha=0.1)
    ax.set_aspect(1.0)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if fig_fpath:
        plt.savefig(fig_fpath)
    if show:
        plt.show()
