import ensembler
import ensembler.refinement
from ensembler.param_parsers import parse_api_params_string, eval_quantity_string
import simtk.unit as unit

helpstring_header = """\
Refine models by molecular dynamics simulation with implicit solvent (generalized Born surface
area), using OpenMM.

Each MD simulation is preceded by a steepest-descent energy minimzation.

The final refined structure is written as "implicit-refined.pdb.gz" in the model directory.

MPI-enabled.

Options:"""

helpstring_unique_options = [
    """\
  --openmm_platform <platform>      Specify the OpenMM Platform to use {CUDA|OpenCL|CPU|Reference}
                                    (default is to auto-select the fastest platform availble).""",

    """\
  --gpupn <gpupn>                   If using GPUs, specify how many are available per node
                                    [default: 1].""",

    """\
  --simlength <simlength>           Simulation length (ps) [default: 100.0].""",

    """\
  --retry_failed_runs               Retry simulation runs which previously failed, e.g. due to bad
                                    inter-atom contacts.""",

    """\
  --ff <ffname>                     OpenMM force field name [default: amber99sbildn]
                                    See OpenMM documentation for other ff options""",

    """\
  --water_model <modelname>         OpenMM water model name [default: tip3p]
                                    See OpenMM documentation for other water model options""",

    """\
  --api_params <params>             See API documentation for
                                    ensembler.refinement.refine_implicit_md""",
]

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
  --template_seqid_cutoff <cutoff>  Select only templates with sequence identity (percentage)
                                    greater than the given cutoff.""",

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

    if args['--templates']:
        templates = args['--templates'].split(',')
    else:
        templates = False

    if args['--template_seqid_cutoff']:
        template_seqid_cutoff = float(args['--template_seqid_cutoff'])
    else:
        template_seqid_cutoff = False

    if args['--gpupn']:
        gpupn = int(args['--gpupn'])
    else:
        gpupn = 1

    if args['--simlength']:
        sim_length = eval_quantity_string(args['--simlength'])
    else:
        sim_length = 100.0 * unit.picoseconds

    if args['--verbose']:
        loglevel = 'debug'
    else:
        loglevel = 'info'

    if args['--api_params']:
        api_params = parse_api_params_string(args['--api_params'])
    else:
        api_params = {}

    ensembler.refinement.refine_implicit_md(
        openmm_platform=args['--openmm_platform'],
        gpupn=gpupn,
        sim_length=sim_length,
        process_only_these_targets=targets,
        process_only_these_templates=templates,
        template_seqid_cutoff=template_seqid_cutoff,
        retry_failed_runs=args['--retry_failed_runs'],
        ff=args['--ff'],
        implicit_water_model=args['--water_model'],
        verbose=args['--verbose'],
        **api_params
    )