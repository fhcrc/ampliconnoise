import argparse
import importlib
import sys

import clean
import split
import truncate

SUBCOMMANDS = ('clean', 'split', 'truncate')


def main(args=sys.argv[1:]):
    """
    Parses arguments, passes execution to the appropriate subcommand.
    """
    parser = argparse.ArgumentParser(description="""Prepare data for use with
AmpiclonNoise""")
    subparsers = parser.add_subparsers(title='Commands')

    for command_name in SUBCOMMANDS:
        p = subparsers.add_parser(command_name)
        p_module = importlib.import_module('.' + command_name, __package__)
        p_module.build_parser(p)
        p.set_defaults(func=p_module.main)

    parsed_args = parser.parse_args(args)
    parsed_args.func(parsed_args)

