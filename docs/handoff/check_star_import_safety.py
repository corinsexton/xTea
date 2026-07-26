#!/usr/bin/env python
"""Is it safe to delete a `from <mod> import *` line?

Written for trim step 4 (2026-07-25), where `x_TEA_main.py` had 20 star-imports and
one of them -- `from x_TEI_locator import *` -- was to be deleted along with the -C/-D
branches. pyflakes cannot help here: it sees any star-import and reports "unable to
detect undefined names", then stops looking (verified: 19 such lines for x_TEA_main.py,
zero actual undefined-name findings). So a NameError in a branch no test exercises
would not be caught by linting, and `import x_TEA_main` succeeding proves nothing --
the failure would only appear at the moment that branch runs.

What this does instead, entirely from the module namespaces as Python actually builds
them:

  1. import the target module, and separately import every module it star-imports;
  2. compute UNIQUE(victim) = names the victim contributes that NO other star-import
     of the target also contributes (the names that would genuinely disappear);
  3. walk the target's AST for every global name it *reads*, and report the
     intersection.

An empty intersection means deleting the import cannot introduce a NameError. A
non-empty one lists exactly which references have to be re-homed first, with the line
numbers where they are read.

The AST walk deliberately over-collects: it takes every ast.Name in Load context plus
the root of every attribute/call chain, at any nesting depth, and subtracts only names
bound in the same function scope (params, assignments, imports, comprehension targets)
plus module-level bindings and builtins. Over-collecting is the safe direction -- a
false positive costs a manual look, a false negative ships a NameError.

Usage:
    cd xtea && python ../docs/handoff/check_star_import_safety.py x_TEA_main x_TEI_locator

Exit status: 0 if the deletion is safe, 1 if any reference would break.
"""
import ast
import builtins
import importlib
import os
import sys

# xtea/ modules import each other by bare name, so the package dir must be on the
# path. sys.path[0] is this script's own directory (docs/handoff), not the cwd.
sys.path.insert(0, os.getcwd())


def star_imported_modules(path):
    """Every module named in a top-level `from X import *`."""
    tree = ast.parse(open(path).read())
    mods = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and any(a.name == "*" for a in node.names):
            mods.append(node.module)
    return mods


def exported_names(modname):
    """The names `from modname import *` would actually bind."""
    mod = importlib.import_module(modname)
    if hasattr(mod, "__all__"):
        return set(mod.__all__)
    return {n for n in vars(mod) if not n.startswith("_")}


def read_names(path):
    """Global names the module reads, as {name: [lineno, ...]}.

    Over-collects on purpose (see module docstring): only names bound in the
    enclosing function scope, at module level, or in builtins are excluded.
    """
    tree = ast.parse(open(path).read())

    module_bound = set()
    for node in tree.body:
        for sub in ast.walk(node):
            if isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                module_bound.add(sub.name)
            elif isinstance(sub, ast.Name) and isinstance(sub.ctx, (ast.Store, ast.Del)):
                module_bound.add(sub.id)
            elif isinstance(sub, ast.alias):
                module_bound.add((sub.asname or sub.name).split(".")[0])

    reads = {}

    def scope_bound(fn):
        bound = set()
        for a in list(fn.args.posonlyargs) + list(fn.args.args) + list(fn.args.kwonlyargs):
            bound.add(a.arg)
        if fn.args.vararg:
            bound.add(fn.args.vararg.arg)
        if fn.args.kwarg:
            bound.add(fn.args.kwarg.arg)
        for sub in ast.walk(fn):
            if isinstance(sub, ast.Name) and isinstance(sub.ctx, (ast.Store, ast.Del)):
                bound.add(sub.id)
            elif isinstance(sub, ast.alias):
                bound.add((sub.asname or sub.name).split(".")[0])
            elif isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                bound.add(sub.name)
            elif isinstance(sub, ast.ExceptHandler) and sub.name:
                bound.add(sub.name)
        return bound

    functions = [n for n in ast.walk(tree)
                 if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
    bound_by_fn = {id(fn): scope_bound(fn) for fn in functions}

    # map each Name node to its innermost enclosing function, so the right
    # local-binding set is subtracted
    owner = {}
    for fn in functions:
        for sub in ast.walk(fn):
            if isinstance(sub, ast.Name):
                owner[id(sub)] = fn      # walk order: innermost assignment wins last

    for node in ast.walk(tree):
        if not (isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)):
            continue
        name = node.id
        if name in module_bound or hasattr(builtins, name):
            continue
        fn = owner.get(id(node))
        if fn is not None and name in bound_by_fn[id(fn)]:
            continue
        reads.setdefault(name, []).append(node.lineno)
    return reads


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    target, victim = sys.argv[1], sys.argv[2]
    path = target + ".py"

    stars = star_imported_modules(path)
    if victim not in stars:
        sys.exit(f"{target}.py does not star-import {victim} (found: {stars})")

    victim_names = exported_names(victim)
    survivors = set()
    for m in stars:
        if m != victim:
            survivors |= exported_names(m)

    unique = victim_names - survivors
    reads = read_names(path)
    at_risk = {n: sorted(set(ls)) for n, ls in reads.items() if n in unique}

    print(f"{victim} exports {len(victim_names)} names; "
          f"{len(unique)} are unique to it (not re-provided by the other "
          f"{len(stars) - 1} star-imports)")
    print(f"{target}.py reads {len(reads)} global names")
    if not at_risk:
        print(f"\nSAFE: deleting `from {victim} import *` breaks no reference in {target}.py")
        return 0
    print(f"\nNOT SAFE -- {len(at_risk)} reference(s) would become undefined:")
    for n in sorted(at_risk):
        print(f"  {n}: line(s) {', '.join(map(str, at_risk[n]))}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
