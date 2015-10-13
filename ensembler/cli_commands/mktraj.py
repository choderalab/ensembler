from ensembler.tools.mktraj import MkTraj
from ensembler.core import get_targets

helpstring_header = """\
Write an XTC trajectory containing Ensembler models.
For each target, this function also outputs a PDB topology file and a CSV file containing model IDs
and other data such as target-template sequence identity.

Options."""

helpstring_nonunique_options = [
    """\
  --targetsfile <targetsfile>  File containing a list of target IDs to work on (newline-separated).
                               Comment targets out with "#".""",

    """\
  --targets <target>           Define one or more target IDs to work on (comma-separated), e.g.
                               "--targets ABL1_HUMAN_D0,SRC_HUMAN_D0" (default: all targets)""",

    """\
  --modeling_stage <stage>     Define the Ensembler modeling stage for which a trajectory should be
                               made.
                               Options:
                                 auto - select most advanced stage for which models have been generated
                                 build_models
                                 refine_implicit_md
                                 refine_explicit_md
                               Default: auto""",
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_nonunique_options)])
docopt_helpstring = ''


def dispatch(args):
    if args['--targetsfile']:
        with open(args['--targetsfile'], 'r') as targetsfile:
            targets = [line.strip() for line in targetsfile.readlines() if line[0] != '#']
    elif args['--targets']:
        targets = args['--targets'].split(',')
    else:
        targets = [target.id for target in get_targets()]

    if args['--modeling_stage'] == 'auto':
        modeling_stage = None
    elif args['--modeling_stage']:
        modeling_stage = args['--modeling_stage']
    else:
        modeling_stage = None

    for target in targets:
        MkTraj(targetid=target, ensembler_stage=modeling_stage)
