import scipy.stats
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from copy import deepcopy

def optional_ax(original_function):
    def possibly_new_ax(*args, **kwargs):
        ax_given = kwargs.get('ax')
        if not ax_given:
            fig, ax = plt.subplots(figsize=(12, 8))
            kwargs['ax'] = ax
        
        result = original_function(*args, **kwargs)
        
        figure_file_name = kwargs.get('save_as')
        if figure_file_name:
            assert not ax_given, 'Can\'t give ax and figure_file_name'
            fig.savefig(figure_file_name, bbox_inches='tight')
            plt.close(fig)
        
        return result
    
    return possibly_new_ax

def add_commas_to_yticks(ax):
    def commas_formatter(x, pos):
        return '{0:,}'.format(int(x))
    tick_formatter = matplotlib.ticker.FuncFormatter(commas_formatter)
    ax.yaxis.set_major_formatter(tick_formatter)

def enhanced_scatter(xs, ys, ax,
                     color_by_density=True,
                     do_fit=True,
                     show_p_value=True,
                     hists_height=0,
                    ):
    same_lists = np.allclose(xs, ys)

    if color_by_density and not same_lists:
        indices = np.arange(len(xs))
        np.random.shuffle(indices)
        random_indices = indices[:1000]

        sampled_points = np.vstack([xs[random_indices], ys[random_indices]])
        points = np.vstack([xs, ys])
        kernel = scipy.stats.gaussian_kde(sampled_points)
        colors = kernel(points)
    else:
        colors = np.ones_like(xs)

    if same_lists:
        do_fit = False

    kwargs = {'cmap': matplotlib.cm.jet,
              's': 4,
              'linewidths' : (0.1,),
             }

    ax.scatter(xs, ys, c=colors, **kwargs)

    if do_fit:
        fit = np.polyfit(xs, ys, 1)
        beta, _ = fit
        fit_function = np.poly1d(fit)
        x_lims = ax.get_xlim()
        ax.plot(x_lims, fit_function(x_lims), color='black', alpha=0.5)
        ax.set_xlim(min(xs), max(xs))
        
        ax.annotate(r'$\beta$ = {:0.2f}'.format(beta),
                            xy=(1, 0),
                            xycoords='axes fraction',
                            xytext=(-10, 30),
                            textcoords='offset points',
                            horizontalalignment='right',
                           )
    
    r, p = scipy.stats.pearsonr(xs, ys)
    if show_p_value:
        text = 'r = {:0.2f}, p={:0.2e}'.format(r, p)
    else:
        text = 'r = {:0.2f}'.format(r)

    ax.annotate(text,
                xy=(1, 0),
                xycoords='axes fraction',
                xytext=(-10, 15),
                textcoords='offset points',
                horizontalalignment='right',
               )

    if hists_height > 0:
        ax_x = ax.twinx()
        ax_x.hist(xs, alpha=0.3, histtype='step', bins=100)
        y_min, y_max = ax_x.get_ylim()
        ax_x.set_ylim(ymax=y_max / hists_height)

        ax_y = ax.twiny()
        ax_y.hist(ys, alpha=0.3, histtype='step', bins=100, orientation='horizontal')
        x_min, x_max = ax_y.get_xlim()
        ax_y.set_xlim(xmax=x_max / hists_height)

        ax_x.set_yticks([])
        ax_y.set_xticks([])
        
def draw_diagonal(ax):
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='black', alpha=0.5)


def plot_cdf(ax, data, **kw_args):
    data_copy = deepcopy(data)
    data_copy.sort()
    
    # Large data, many duplicate points -> Slow to graph. Dedup here.
    x = list(set(data_copy))
    x.sort()
    y = []
    xiter = iter(x)
    xx = next(xiter)
    for i, dd in enumerate(data_copy):
        if dd > xx:
            xx = next(xiter)
            y.append(i/float(len(data_copy)))  # would be i-1 for one-based arrays
    y.append(1.0)
    
    # Add start and double points strategically for nice plotting
    x = [x[0]] + 2*x
    x.sort()
    y = [0.0] + [yy for tup in zip(y[:-1], y[1:]) for yy in tup] + [y[-1]]*2
    
    ax.plot(x, y, **kw_args)


