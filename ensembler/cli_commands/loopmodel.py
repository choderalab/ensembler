import ensembler
import ensembler.modeling

helpstring_header = """\
Use Rosetta loopmodel to model missing loops in template structures.

Options:"""

helpstring_unique_options = [
    """\
  --templates <template>       Define one or more comma-separated template IDs to work on (e.g. "--templates ABL1_HUMAN_D0_1OPL_A") (default: all templates)""",

    """\
  --overwrite_structures       Overwrite structure files""",

    """\
  -v --verbose                 """,
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_unique_options)])
docopt_helpstring = '\n\n'.join(helpstring_unique_options)


def dispatch(args):
    if args['--verbose']:
        loglevel = 'debug'
    else:
        loglevel = 'info'

    if args['--templates']:
        templates = args['--templates'].split(',')
    else:
        templates = False

    ensembler.modeling.model_template_loops(process_only_these_templates=templates, overwrite_structures=args['--overwrite_structures'], loglevel=loglevel)