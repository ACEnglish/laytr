"""
Helper functions for plotting a som
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colorbar
from matplotlib.lines import Line2D
import matplotlib.ticker as mticker
from matplotlib.patches import RegularPolygon, Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable


def make_hex_plot(som, hue=None, hue_label="Distance", color_map=cm.Blues, hue_count_ticks=False):
    """
    Creates the hex plot for a som. Hue is a matrix of the same shape as the SOM.
    By default, hue will be a som.distance_map()
    """
    xx, yy = som.get_euclidean_coordinates()
    if hue is None:
        hue = som.distance_map()

    f = plt.figure(figsize=(10,10))
    ax = f.add_subplot(111)

    ax.set_aspect('equal')

    # iteratively add hexagons
    for i in range(hue.shape[0]):
        for j in range(hue.shape[1]):
            wy = yy[(i, j)] * np.sqrt(3) / 2
            hex = RegularPolygon((xx[(i, j)], wy),
                                 numVertices=6,
                                 radius=.95 / np.sqrt(3),
                                 facecolor=color_map(hue[i, j]),
                                 alpha=.4,
                                 edgecolor='gray')
            ax.add_patch(hex)
    xrange = np.arange(hue.shape[0])
    yrange = np.arange(hue.shape[1])
    plt.xticks(xrange - 0.5, xrange)
    plt.yticks(yrange * np.sqrt(3) / 2, yrange)

    x1, x2 = ax.get_xlim()
    plt.xlim(x1-0.5, x2 + 1)
    y1, y2 = ax.get_ylim()
    plt.ylim(y1 -0.5, y2 + 0.5)

    divider = make_axes_locatable(plt.gca())
    ax_cb = divider.new_horizontal(size="5%", pad=0.05)
    cb1 = colorbar.ColorbarBase(ax_cb, cmap=color_map,
                                orientation='vertical', alpha=.4)
    cb1.ax.get_yaxis().labelpad = 16
    cb1.ax.set_ylabel(hue_label,
                      rotation=270, fontsize=16)
    # Make the ticks correspond to values
    if hue_count_ticks:
        ticks = cb1.get_ticks().tolist()
        cb1.set_ticks(ticks, labels=[int(hue.max() * _) for _ in ticks])

    plt.gcf().add_axes(ax_cb)
    return f


def add_places(som, fig, neurons, labels,
             marker_map={'patho':'1', 'codis':'+', 'inter':'x', 'other':'o'},
             default_colors=['C3', 'C1', 'C2', 'C0']):
    """
    Adds markers to a hex plot
    """
    colors = {}
    for key, col in zip(marker_map.keys(), default_colors):
        colors[key] = col

    ax = fig.get_axes()[0]
    # Adding plots to it
    for cnt, w in enumerate(neurons):
        # place a marker on the winning position for the sample xx
        wx, wy = som.convert_map_to_euclidean(tuple(w))
        wy = wy * np.sqrt(3) / 2
        ax.plot(wx, wy,
                marker_map[labels[cnt]],
                markerfacecolor='None',
                markeredgecolor=colors[labels[cnt]],
                markersize=12,
                markeredgewidth=2)

    legend_elements = []
    for key in marker_map:
        legend_elements.append(Line2D([0], [0], marker=marker_map[key], color=colors[key], label=key,
                       markerfacecolor='w', markersize=14, linestyle='None', markeredgewidth=2))

    ax.legend(handles=legend_elements, bbox_to_anchor=(0.1, 1.08), loc='upper left',
              borderaxespad=0., ncol=len(marker_map), fontsize=14)
    return fig
