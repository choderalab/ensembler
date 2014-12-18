import sys
sys.path.pop(0)
sys.path.pop(0)
import ensembler
import ensembler.modeling

helpstring_header = """\
Models a set of target sequences onto a set of template structures using Modeller.

Options:"""

helpstring_unique_options = [
    """\
  --templates <template>...   Define one or more template IDs to work on (e.g. "--templates ABL1_HUMAN_D0_1OPL_A") (default: all templates)""",
]

helpstring_nonunique_options = [
    """\
  --targets <target>...       Define one or more target IDs to work on (e.g. "--targets ABL1_HUMAN_D0 --targets SRC_HUMAN_D0") (default: all targets)""",

    """\
  -v --verbose                """,
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_unique_options + helpstring_nonunique_options)])
docopt_helpstring = '\n\n'.join(helpstring_unique_options)


def dispatch(args):
    if args['--help']:
        print helpstring
        return
    if args['--verbose']:
        loglevel = 'debug'
    else:
        loglevel = 'info'
    ensembler.modeling.build_models(process_only_these_targets=args['--targets'], process_only_these_templates=args['--templates'], loglevel=loglevel)