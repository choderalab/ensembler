import ensembler
import ensembler.refinement

helpstring_header = """\
Conducts explicit-solvent MD refinement on a set of models.

Options."""

helpstring_nonunique_options = [
    """\
  --openmm_platform <platform>    Choose the OpenMM Platform to use (choices: CUDA, OpenCL, CPU, Reference) [default: CUDA].""",
    """\
  --gpupn <gpupn>                  If using GPUs, specify how many are available per node [default: 1].""",
  #   """\
  # --targetsfile <targetsfile>     File containing a list of newline-separated target IDs to work on. Comment targets out with "#".""",
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
    if args['--verbose']:
        loglevel = 'debug'
    else:
        loglevel = 'info'

    ensembler.refinement.refine_implicit_md(openmm_platform=args['--openmm_platform'], gpupn=int(args['--gpupn']), process_only_these_targets=args['--targets'], process_only_these_templates=args['--templates'], verbose=args['--verbose'])