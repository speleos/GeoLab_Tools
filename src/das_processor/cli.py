import argparse
from .config import init_config


def main():
    parser = argparse.ArgumentParser(prog="das-processor")
    sub = parser.add_subparsers(dest="command", required=True)

    p_init = sub.add_parser("init", help="Copy config templates to a folder")
    p_init.add_argument("dest", nargs="?", default=".")
    p_init.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()
    if args.command == "init":
        init_config(args.dest, args.overwrite)
