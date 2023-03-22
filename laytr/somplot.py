import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon, Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm, colorbar
from matplotlib.lines import Line2D

def plot_hex(som, data, target,
             marker_map={'patho':'1', 'codis':'+', 'inter':'x', 'other':'o'},
             default_colors=['C3', 'C1', 'C2', 'C0'],
             umatrix = None,
             umatrix_label='distance from neurons in the neighbourhood',
             color_map=cm.Blues
             ):
    """
    Groups need to map to markers and colors
    data - the features we're analyzing
    target - names for marker_map
    """
    colors = {}
    for key, col in zip(marker_map.keys(), default_colors):
        colors[key] = col

    xx, yy = som.get_euclidean_coordinates()
    if umatrix is None:
        umatrix = som.distance_map()
    weights = som.get_weights()

    f = plt.figure(figsize=(10,10))
    ax = f.add_subplot(111)

    ax.set_aspect('equal')

    # iteratively add hexagons
    for i in range(weights.shape[0]):
        for j in range(weights.shape[1]):
            wy = yy[(i, j)] * np.sqrt(3) / 2
            hex = RegularPolygon((xx[(i, j)], wy),
                                 numVertices=6,
                                 radius=.95 / np.sqrt(3),
                                 facecolor=color_map(umatrix[i, j]),
                                 alpha=.4,
                                 edgecolor='gray')
            ax.add_patch(hex)

    ret_winners = []
    for cnt, x in enumerate(data):
        # getting the winner
        w = som.winner(x)
        ret_winners.append(w)
        # place a marker on the winning position for the sample xx
        wx, wy = som.convert_map_to_euclidean(w)
        wy = wy * np.sqrt(3) / 2
        plt.plot(wx, wy,
                 marker_map[target[cnt]],
                 markerfacecolor='None',
                 markeredgecolor=colors[target[cnt]],
                 markersize=12,
                 markeredgewidth=2)

    xrange = np.arange(weights.shape[0])
    yrange = np.arange(weights.shape[1])
    plt.xticks(xrange-.5, xrange)
    plt.yticks(yrange * np.sqrt(3) / 2, yrange)

    divider = make_axes_locatable(plt.gca())
    ax_cb = divider.new_horizontal(size="5%", pad=0.05)
    cb1 = colorbar.ColorbarBase(ax_cb, cmap=color_map,
                                orientation='vertical', alpha=.4)
    cb1.ax.get_yaxis().labelpad = 16
    cb1.ax.set_ylabel(umatrix_label,
                      rotation=270, fontsize=16)
    plt.gcf().add_axes(ax_cb)

    legend_elements = []
    for key in marker_map:
        legend_elements.append(Line2D([0], [0], marker=marker_map[key], color=colors[key], label=key,
                       markerfacecolor='w', markersize=14, linestyle='None', markeredgewidth=2))

    ax.legend(handles=legend_elements, bbox_to_anchor=(0.1, 1.08), loc='upper left',
              borderaxespad=0., ncol=len(marker_map), fontsize=14)
    return ret_winners

