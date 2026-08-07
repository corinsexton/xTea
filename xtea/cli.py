"""Console-script entry point: the `xtea2` command.

`pyproject.toml` maps `xtea2` to `xtea.cli:main`. Importing this module imports the `xtea` package
first, whose `__init__` puts this directory on `sys.path` -- so the flat `import run_xtea` below
resolves exactly as it does when `xtea/run_xtea.py` is run as a script, and the two invocations
execute the same code.

Nothing else belongs here. Argument parsing, defaults and the whole pipeline stay in `run_xtea`,
so `xtea2 --help` and `python xtea/run_xtea.py --help` cannot drift apart.
"""
import sys


def main():
    import run_xtea
    return run_xtea.main()


if __name__ == "__main__":
    sys.exit(main())
