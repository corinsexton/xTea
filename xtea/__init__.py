"""xTea2 — TE insertion calling from BAM/CRAM in one process.

WHY THIS FILE MUTATES sys.path, WHICH IS NORMALLY A SMELL
---------------------------------------------------------
The 29 modules in this directory import each other FLATLY -- `import global_values`, not
`from xtea import global_values`. That works when `run_xtea.py` is executed as a script, because
Python puts the script's own directory on `sys.path`. It does NOT work when the same modules are
imported as `xtea.<name>` from an installed package, and it is the only thing standing between this
tree and a normal `pip install`.

There were two ways to fix it. Rewriting all 29 modules to package-relative imports is the textbook
answer and was NOT taken: this pipeline is currently bit-reproducible (three unpinned WGS
replicates produced byte-identical call sets, array 50062602), that property is the strongest thing
it has, and a 29-module import refactor would put every byte-identity claim back behind a WGS
re-gate for a packaging convenience. The cost is not worth it.

So instead the package restores the ONE condition those imports already depend on -- this directory
being importable -- and changes no module's import statements at all. `import xtea` and running
`python xtea/run_xtea.py` then resolve every sibling identically, which is precisely the point: the
installed entry point and the in-repo script execute the same code down to the import graph.

The insert is idempotent and appends rather than prepends, so a caller who has deliberately shadowed
one of these module names keeps their override.

If those imports are ever converted properly, DELETE THIS BLOCK -- it is scaffolding for a specific
layout, not a design.
"""
import os
import sys

__version__ = "2.0.0"

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.append(_HERE)
