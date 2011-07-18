import argparse
import sys


SUBCOMMANDS = ('clean', 'redup', 'split', 'pyronoise', 'seqnoise', 'truncate', 'wfasta', 'sff2raw')


def main(args=sys.argv[1:]):
    """
    Parses arguments, passes execution to the appropriate subcommand.
    """
    parser = argparse.ArgumentParser(description="""Prepare data for use with
AmpliconNoise""")
    subparsers = parser.add_subparsers(title='Commands', help="Valid commands")

    modules = __import__(__package__, fromlist=SUBCOMMANDS)

    for command_name in SUBCOMMANDS:
        p_module = getattr(modules, command_name)
        p = p_module.build_parser(subparsers)
        p.set_defaults(func=p_module.main)

    parsed_args = parser.parse_args(args)
    parsed_args.func(parsed_args)
