# import sys
# sys.path.pop(0)
# sys.path.pop(0)
import ensembler
import ensembler.modeling

helpstring_header = """\
Models a set of target sequences onto a set of template structures using Modeller.

Options:"""

helpstring_nonunique_options = [
    """\
  --targets <target>       Define one or more comma-separated target IDs to work on (e.g. "--targets ABL1_HUMAN_D0,SRC_HUMAN_D0") (default: all targets)""",

    """\
  --templates <template>   Define one or more comma-separated template IDs to work on (e.g. "--templates ABL1_HUMAN_D0_1OPL_A") (default: all templates)""",

    """\
  -v --verbose                """,
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_nonunique_options)])
docopt_helpstring = ''


def dispatch(args):
    if args['--targets']:
        targets = args['--targets'].split(',')
    else:
        targets = False
    if args['--templates']:
        templates = args['--templates'].split(',')
    else:
        templates = False

    if args['--verbose']:
        loglevel = 'debug'
    else:
        loglevel = 'info'

    ensembler.modeling.build_models(process_only_these_targets=targets, process_only_these_templates=templates, loglevel=loglevel)