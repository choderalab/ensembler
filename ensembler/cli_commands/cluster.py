import ensembler
import ensembler.modeling

helpstring_header = """\
Filter out non-unique models by clustering on RMSD.

Unique models are designated by writing an empty file named "unique_by_clustering" in their model
directory.

Runs serially.

Options:"""

helpstring_unique_options = [
    """\
  --cutoff <cutoff>               Minimum distance cutoff for RMSD-based clustering (nm)
                                  (default: 0.06)""",
]

helpstring_nonunique_options = [
    """\
  --targetsfile <targetsfile>  File containing a list of target IDs to work on (newline-separated).
                               Comment targets out with "#".""",

    """\
  --targets <target>           Define one or more target IDs to work on (comma-separated), e.g.
                               "--targets ABL1_HUMAN_D0,SRC_HUMAN_D0" (default: all targets)""",

    """\
  -v --verbose                 """,
]

helpstring = '\n\n'.join([helpstring_header,  '\n\n'.join(helpstring_unique_options), '\n\n'.join(helpstring_nonunique_options)])
docopt_helpstring = '\n\n'.join(helpstring_unique_options)

def dispatch(args):
    if args['--targetsfile']:
        with open(args['--targetsfile'], 'r') as targetsfile:
            targets = [line.strip() for line in targetsfile.readlines() if line[0] != '#']
    elif args['--targets']:
        targets = args['--targets'].split(',')
    else:
        targets = False

    cutoff = ensembler.utils.set_arg_with_default(args['--cutoff'], default_arg=0.06)

    if args['--verbose']:
        loglevel = 'debug'
    else:
        loglevel = 'info'

    ensembler.modeling.cluster_models(process_only_these_targets=targets, verbose=args['--verbose'], cutoff=cutoff)