import argparse
import importlib
import sys


SUBCOMMANDS = ('clean', 'split', 'truncate')


def main(args=sys.argv[1:]):
    """
    Parses arguments, passes execution to the appropriate subcommand.
    """
    parser = argparse.ArgumentParser(description="""Prepare data for use with
AmpiclonNoise""")
    subparsers = parser.add_subparsers(title='Commands', help="Valid commands")

    for command_name in SUBCOMMANDS:
        p_module = importlib.import_module('.' + command_name, __package__)
        p = p_module.build_parser(subparsers)
        p.set_defaults(func=p_module.main)

    parsed_args = parser.parse_args(args)
    parsed_args.func(parsed_args)
