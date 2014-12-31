# import sys
# sys.path.pop(0)
# sys.path.pop(0)
import ensembler
import ensembler.modeling

helpstring_header = """\
Pairwise alignment of target sequences onto template sequences.

Options:"""

helpstring_unique_options = [
    """\
  --targets <target>...       Define one or more target IDs to work on (e.g. "--targets ABL1_HUMAN_D0 --targets SRC_HUMAN_D0") (default: all targets)""",

    """\
  --templates <template>...   Define one or more template IDs to work on (e.g. "--templates ABL1_HUMAN_D0_1OPL_A") (default: all templates)""",

    """\
  -v --verbose                """,
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_unique_options)])
docopt_helpstring = '\n\n'.join(helpstring_unique_options)


def dispatch(args):
    if args['--verbose']:
        loglevel = 'debug'
    else:
        loglevel = 'info'
    ensembler.modeling.align_targets_and_templates(process_only_these_targets=args['--targets'], process_only_these_templates=args['--templates'], loglevel=loglevel)