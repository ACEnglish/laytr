"""

laytr - Library for variant benchmarking stratification

See `help()` of specific functions / objects for details

VariantRecord methods:
:meth:`entry_boundaries`
"""
__version__ = '0.0.1'

from laytr.kfeat import (
    kfeat,
    kfeats,
    get_features,
    iter_regions,
    regions_to_kmers,
)

from laytr.somplot import (
    make_hex_plot,
    add_places,
)

