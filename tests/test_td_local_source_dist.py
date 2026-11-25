"""Z4b: the local-source exclusion must cover a site's whole flank (--flklen, default 3000)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "xtea"))

import global_values  # noqa: E402


def test_local_source_distance_covers_default_flank_plus_insert():
    assert global_values.MIN_POLYMORPHIC_SOURCE_DIST >= 3000 + 1000
