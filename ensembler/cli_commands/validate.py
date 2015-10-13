import ensembler
import ensembler.validation

helpstring_header = """\
Calculate model quality.

For each target, this outputs a text file named
``models/[targetid]/validation_scores_sorted-[method]-[ensembler_stage]`` which contains a list of
model IDs sorted by validation score. This can be used by the subsequent ``package_models`` command
to filter out models below a specified quality threshold.

Typically, this should be run after models have been refined to the desired extent (e.g. after
implicit or explicit MD refinement)

More detailed validation results are written to the individual model directories.

MPI-enabled.

Options."""

helpstring_unique_options = [
    """\
  --modeling_stage <stage>     Define the Ensembler modeling stage at which validation should be run
                               Options:
                                 auto - select most advanced stage for which models have been built
                                 build_models
                                 refine_implicit_md
                                 refine_explicit_md
                               Default: auto""",
]

helpstring_nonunique_options = [
    """\
  --targetsfile <targetsfile>  File containing a list of target IDs to work on (newline-separated).
                               Comment targets out with "#".""",

    """\
  --targets <target>           Define one or more target IDs to work on (comma-separated), e.g.
                               "--targets ABL1_HUMAN_D0,SRC_HUMAN_D0" (default: all targets)""",

    """\
  --method <method>            Validation method to use
                               Options:
                                 molprobity
                               Default: molprobity""",

    """\
  -v --verbose                 """,
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_unique_options), '\n\n'.join(helpstring_nonunique_options)])
docopt_helpstring = '\n\n'.join(helpstring_unique_options)


def dispatch(args):
    if args['--targetsfile']:
        with open(args['--targetsfile'], 'r') as targetsfile:
            targets = [line.strip() for line in targetsfile.readlines() if line[0] != '#']
    elif args['--targets']:
        targets = args['--targets'].split(',')
    else:
        targets = False

    if args['--method']:
        method = args['--method']
    else:
        method = 'molprobity'

    if args['--modeling_stage']:
        if args['--modeling_stage'] == 'auto':
            modeling_stage = None
        else:
            modeling_stage = args['--modeling_stage']
    else:
        modeling_stage = None

    if args['--verbose']:
        loglevel = 'debug'
    else:
        loglevel = 'info'

    if method == 'molprobity':
        ensembler.validation.molprobity_validation_multiple_targets(
            targetids=targets,
            modeling_stage=modeling_stage,
            loglevel=loglevel,
        )
