import ensembler
import ensembler.refinement

helpstring_header = """\
Solvates models.

Options."""

helpstring_nonunique_options = [
    """\
  --targets <target>...           Define one or more target IDs to work on (e.g. "--targets ABL1_HUMAN_D0 --targets SRC_HUMAN_D0") (default: all targets)""",
    """\
  --templates <template>...       Define one or more template IDs to work on (e.g. "--templates ABL1_HUMAN_D0_1OPL_A") (default: all templates)""",
    """\
  -v --verbose                    """,
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_nonunique_options)])
docopt_helpstring = ''

def dispatch(args):
    ensembler.refinement.solvate_models(process_only_these_targets=args['--targets'], process_only_these_templates=args['--templates'])
    ensembler.refinement.determine_nwaters(process_only_these_targets=args['--targets'], process_only_these_templates=args['--templates'])