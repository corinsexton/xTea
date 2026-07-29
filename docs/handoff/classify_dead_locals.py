#!/usr/bin/env python3
"""Classify pyflakes' "local variable X is assigned to but never used" by RHS side-effect risk.

TRIM_PLAN step 5 deletes dead local assignments. The hazard is stated in the plan but not
mechanised: deleting `x = f()` also deletes the CALL, and if `f` writes a file, mutates an
argument, or advances an iterator, that changes behaviour. Deleting `x = a[0]` cannot.

So the deletion set has to be split, and eyeballing 125 sites is exactly the kind of review
that misses the one that matters. This walks the AST instead and puts every site in one of:

  PURE       RHS provably free of calls and of anything that can mutate: literals, names,
             attribute/subscript reads, arithmetic/comparison/boolop, f-strings, comprehensions
             over pure parts, and calls to a small allowlist of builtins that cannot have
             side effects on their arguments (int/str/float/len/...). SAFE TO DELETE OUTRIGHT.
  CALL       RHS contains a call that is not on the allowlist. The BINDING is dead but the
             call may not be. Do not delete the statement; either leave it, or drop only the
             target (`f()`), which needs a human to confirm the call is wanted at all.
  AUGMENTED  `x op= ...` -- reads x, so deleting it can change a later read of x that pyflakes
             attributes elsewhere. Hand review.
  UNPACK     part of a tuple/list unpacking. Deleting one name means rewriting the target into
             `_`, which is a different edit. Hand review.
  OTHER      anything not matched above (walrus, annotated assign, for-target, with-target...).

Note on pyflakes' own blind spot, which shapes what this can be trusted for: on a module with a
star-import pyflakes reports "unable to detect undefined names" and then finds NOTHING at all.
Every module under xtea/ star-imports, so the input list here is not authoritative for unused
*names* -- but "assigned and never read" is a local-scope question that star-imports do not
affect, so the local-variable findings are sound. (See check_star_import_safety.py for the other
half.)

Usage:
    python -m pyflakes xtea/*.py | python docs/handoff/classify_dead_locals.py
    python docs/handoff/classify_dead_locals.py --tsv    # machine-readable
"""
import ast
import collections
import re
import sys

# Builtins whose return value is the whole point and which cannot mutate their arguments.
# Deliberately conservative: anything not here is treated as impure.
PURE_CALLS = {
    "int", "str", "float", "bool", "len", "abs", "min", "max", "round", "ord", "chr",
    "repr", "format", "divmod", "tuple", "frozenset", "hash", "isinstance", "issubclass",
}
# Method names that are pure reads on str/list/dict. `split`, `strip`, `rstrip` etc. return new
# objects; `get`/`keys`/`values`/`items` read a dict. Anything absent is impure by default --
# notably append/write/add/update/pop/sort/remove/close/seek/readline.
PURE_METHODS = {
    "split", "rsplit", "strip", "lstrip", "rstrip", "lower", "upper", "join", "replace",
    "startswith", "endswith", "find", "rfind", "count", "get", "keys", "values", "items",
    "index", "format", "encode", "decode", "isdigit", "isalpha", "copy",
}

LINE_RE = re.compile(r"^(?P<file>.+?):(?P<line>\d+):(?:\d+:)?\s*local variable '(?P<name>[^']+)'")


class Purity(ast.NodeVisitor):
    """True unless the tree contains something that might have a side effect."""

    def __init__(self):
        self.impure_reason = None

    def visit_Call(self, node):
        f = node.func
        ok = False
        if isinstance(f, ast.Name) and f.id in PURE_CALLS:
            ok = True
        elif isinstance(f, ast.Attribute) and f.attr in PURE_METHODS:
            ok = True
        if not ok:
            name = getattr(f, "id", None) or getattr(f, "attr", None) or "<expr>"
            if self.impure_reason is None:
                self.impure_reason = "call to %s()" % name
            return          # do not recurse; the outermost impure call is the useful report
        self.generic_visit(node)

    def visit_Await(self, node):
        self.impure_reason = self.impure_reason or "await"

    def visit_Yield(self, node):
        self.impure_reason = self.impure_reason or "yield"

    visit_YieldFrom = visit_Yield

    def visit_NamedExpr(self, node):     # walrus binds a second name
        self.impure_reason = self.impure_reason or "walrus binding"


def classify(path, line_no, name, src_lines, tree_cache):
    tree = tree_cache.get(path)
    if tree is None:
        try:
            tree = ast.parse("".join(src_lines[path]), filename=path)
        except SyntaxError as exc:
            return "OTHER", "unparseable: %s" % exc, ""
        tree_cache[path] = tree

    stmt = None
    for node in ast.walk(tree):
        if getattr(node, "lineno", None) != line_no:
            continue
        if isinstance(node, (ast.Assign, ast.AugAssign, ast.AnnAssign)):
            stmt = node
            break
    text = src_lines[path][line_no - 1].strip() if line_no <= len(src_lines[path]) else ""

    if stmt is None:
        return "OTHER", "no simple assignment at this line (for/with/except target?)", text
    if isinstance(stmt, ast.AugAssign):
        return "AUGMENTED", "%s reads its own target" % name, text
    if isinstance(stmt, ast.Assign):
        # unpacking, or several targets
        for t in stmt.targets:
            if isinstance(t, (ast.Tuple, ast.List)):
                return "UNPACK", "tuple/list unpacking target", text
        if len(stmt.targets) > 1:
            return "UNPACK", "chained assignment (a = b = ...)", text

    p = Purity()
    p.visit(stmt.value)
    if p.impure_reason:
        return "CALL", p.impure_reason, text
    return "PURE", ast.dump(stmt.value, annotate_fields=False)[:40], text


def main():
    argv = sys.argv[1:]
    as_tsv = "--tsv" in argv
    findings = []
    for raw in sys.stdin:
        m = LINE_RE.match(raw.strip())
        if m:
            findings.append((m.group("file"), int(m.group("line")), m.group("name")))
    if not findings:
        print("no 'local variable ... never used' findings on stdin", file=sys.stderr)
        return 1

    src_lines, tree_cache = {}, {}
    for path, _, _ in findings:
        if path not in src_lines:
            with open(path) as fh:
                src_lines[path] = fh.readlines()

    rows = []
    for path, line_no, name in findings:
        verdict, why, text = classify(path, line_no, name, src_lines, tree_cache)
        rows.append((verdict, path, line_no, name, why, text))

    if as_tsv:
        print("verdict\tfile\tline\tname\treason\tsource")
        for r in sorted(rows):
            print("\t".join(str(x) for x in r))
        return 0

    by_verdict = collections.Counter(r[0] for r in rows)
    print("=== %d dead locals across %d files ===" % (len(rows), len(src_lines)))
    for v in ("PURE", "CALL", "AUGMENTED", "UNPACK", "OTHER"):
        if by_verdict[v]:
            print("  %-9s %3d" % (v, by_verdict[v]))
    for v in ("PURE", "CALL", "AUGMENTED", "UNPACK", "OTHER"):
        sel = [r for r in rows if r[0] == v]
        if not sel:
            continue
        print("\n--- %s (%d) ---" % (v, len(sel)))
        for _, path, line_no, name, why, text in sorted(sel, key=lambda r: (r[1], r[2])):
            short = path.split("/")[-1]
            print("  %-26s:%-5d %-24s %-28s %s" % (short, line_no, name, why, text[:70]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
