#!/usr/bin/env python3
"""
laytr main entrypoint
"""
import sys
import argparse
from importlib.metadata import version

from laytr import __version__
from laytr.kfeat import kfeat_main
from laytr.map import map_main

def flat_version(args):
    """Print the version"""
    if len(args) and args[0].count("-v"):
        print(f"laytr {version('truvari')}")
    else:
        print(f"laytr v{__version__}")

TOOLS = {"kfeat": kfeat_main,
         "map": map_main,
}

USAGE = f"""\
laytr v{__version__} Library for variant benchmarking stratification

Available commands:
    kfeat   Create kmer featuration of genomic regions
    map     Map kfeats to a SOM and report neurons
"""

def main():
    """
    Main entrypoint for laytr
    """
    parser = argparse.ArgumentParser(prog="laytr", description=USAGE,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="CMD", choices=TOOLS.keys(), type=str, default=None,
                        help="Command to execute")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()
    args = parser.parse_args()

    TOOLS[args.cmd](args.options)

if __name__ == '__main__':
    main()
