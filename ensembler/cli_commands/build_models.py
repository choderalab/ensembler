import ensembler
import ensembler.modeling

helpstring_header = """\
Generate models by mapping target sequences onto template structures, using Modeller
(salilab.org/modeller).

MPI-enabled.

Options:"""

helpstring_nonunique_options = [
    """\
  --targetsfile <targetsfile>  File containing a list of target IDs to work on (newline-separated).
                               Comment targets out with "#".""",

    """\
  --targets <target>           Define one or more target IDs to work on (comma-separated), e.g.
                               "--targets ABL1_HUMAN_D0,SRC_HUMAN_D0" (default: all targets)""",

    """\
  --templates <template>       Define one or more template IDs to work on (comma-separated), e.g.
                               "--templates ABL1_HUMAN_D0_1OPL_A" (default: all templates)""",

    """\
  -v --verbose                 """,
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