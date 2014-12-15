import ensembler
import ensembler.modeling

helpstring_header = """\
Sorts models by sequence identity, then performs clustering to filter out non-unique models.

Options."""

helpstring_nonunique_options = [
    """\
  --targets <target>...       Define one or more target IDs to work on (e.g. "--targets ABL1_HUMAN_D0 --targets SRC_HUMAN_D0") (default: all targets)""",
    """\
  -v --verbose                """,
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_nonunique_options)])
docopt_helpstring = ''

def dispatch(args):
    if args['--help']:
        print helpstring
        return
    if args['--verbose']:
        loglevel = 'debug'
    else:
        loglevel = 'info'
    ensembler.modeling.sort_by_sequence_identity(process_only_these_targets=args['--targets'], loglevel=loglevel)
    ensembler.modeling.cluster_models(process_only_these_targets=args['--targets'], verbose=args['--verbose'])