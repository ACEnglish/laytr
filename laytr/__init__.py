"""

laytr - Library for variant benchmarking stratification

See `help()` of specific functions / objects for details

VariantRecord methods:
:meth:`entry_boundaries`
"""
__version__ = '0.0.1'

from laytr.kfeat import (
    kfeat,
    kfeat_seqs,
    get_features,
    RegionIter,
    regions_to_kmers,
)

from laytr.map import (
    map_to_som
)

from laytr.somplot import (
    make_hex_plot,
    add_places,
)

