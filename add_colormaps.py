def add_colormaps(colormaps, data_limits=[0.0,1.0], data_break=0.5,
                             colormap_name="JohnDoe"):
    r"""Concatenate colormaps in *colormaps* list and return Normalize object.
    Appends the colormaps in the list *colormaps* adjusting the mapping to the
    colormap such that it maps the data space limits *data_limits* and puts the
    break in the colormaps at *data_break* which is again in data space.  The
    argument *colormap_name* labels the colormap and is passed to the init
    method of `matplotlib.colors.LinearSegmentedColormap`.
    Originally contributed by Damon McDougall
    returns `matplotlib.colors.LinearSegmentedColormap`,
            `matplotlib.colors.Normalize`

    Taken from visclaw by KMT 6/26/17.

    Example
    -------
    This example takes two colormaps whose data ranges from [-1.0, 5.0] where
    the break in the colormaps occurs at 1.0.
    >>> import matplotlib.pyplot as plt
    >>> import clawpack.visclaw.colormaps as colormaps
    >>> import numpy
    >>> cmaps = (plt.get_cmap("BuGn_r"), plt.get_cmap("Blues_r"))
    >>> new_cmap, norm = colormaps.add_colormaps(cmaps, data_limits=[-1.0, 5.0],
        data_break=1.0)
    >>> x = numpy.linspace(-1,1,100)
    >>> y = x
    >>> X, Y = numpy.meshgrid(x, y)
    >>> fig = plt.figure()
    >>> axes = fig.add_subplot(1,1,1)
    >>> plot = axes.pcolor(X, Y, 3.0 * X + 2.0, cmap=new_cmap, norm=norm)
    >>> fig.colorbar(plot)
    >>> plt.show()
    """

    lhs_dict = colormaps[0]._segmentdata
    rhs_dict = colormaps[1]._segmentdata
    new_dict = dict(red=[], green=[], blue=[], alpha=[])

    # Add first colorbar
    for key in rhs_dict:
        val_list = rhs_dict[key]
        for val in val_list:
            new_dict[key].append((val[0] * 0.5, val[1], val[2]))

    if 'alpha' not in list(rhs_dict.keys()):
        new_dict['alpha'].append((0.0,1.0,1.0))

    # Add second colorbar
    for key in lhs_dict:
        val_list = lhs_dict[key]
        for val in val_list:
            new_dict[key].append(((val[0] + 1.0) * 0.5, val[1], val[2]))

    if 'alpha' not in list(lhs_dict.keys()):
        new_dict['alpha'].append((1.0,1.0,1.0))

    N = 256
    gamma = 1.0

    cmap = colors.LinearSegmentedColormap(colormap_name, new_dict, N, gamma)

    # Compute new norm object
    bounds = numpy.empty(N)
    bounds[:int(N / 2)] = numpy.linspace(data_limits[0], data_break, int(N / 2))
    bounds[int(N / 2):] = numpy.linspace(data_break, data_limits[1],
                                                             int(N / 2) + N % 2)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=N)

    return cmap, norm
