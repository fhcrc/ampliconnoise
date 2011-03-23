import argparse
import importlib
import sys

import clean
import split
import truncate

SUBCOMMANDS = ('clean', 'split', 'truncate')
#PARSER_BUILDERS = {
#    'clean': (clean.build_parser, clean.main),
#    'split': (split.build_parser, split.main),
#    'truncate': (truncate.build_parser, truncate.main),
#}

def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="""Prepare data for use with
AmpiclonNoise""")
    subparsers = parser.add_subparsers(title='Commands')

    for command_name in SUBCOMMANDS:
        p = subparsers.add_parser(command_name)
        p_module = importlib.import_module('.' + command_name, __package__)
        p_module.build_parser(p)
        p.set_defaults(func=p_module.main)

#    for command_name, (builder, func) in PARSER_BUILDERS.items():
#
#        p = subparsers.add_parser(command_name)
#        builder(p)
#        p.set_default(func=func)

    parsed_args = parser.parse_args(args)
    parsed_args.func(parsed_args)

