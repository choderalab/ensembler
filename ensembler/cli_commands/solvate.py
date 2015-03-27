import ensembler
import ensembler.refinement

helpstring_header = """\
Determines the number of waters to add when solvating models with explicit water molecules.

The models for each target are to be given the same number of waters.

The function proceeds by first solvating each model individually, given a padding distance
(default: 1 nm). A list of the number of waters added for each model is written to a file
"nwaters.txt" in the "models/[target_id]" directory. A percentile value from the distribution of
the number of waters is selected as the number to use for all models, and this number is written to
the file nwaters-use.txt.

MPI-enabled

Options:"""

helpstring_unique_options = [
    """\
  --padding <padding>          Padding distance for solvation (Angstroms) (default: 10 Angstroms).""",
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
  --ff <ffname>                     OpenMM force field name [default: amber99sbildn]
                                    See OpenMM documentation for other ff options""",

    """\
  --water_model <modelname>         OpenMM water model name [default: tip3p]
                                    See OpenMM documentation for other water model options""",

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

    padding = ensembler.utils.set_arg_with_default(args['--padding'], default_arg=10.0)

    if args['--verbose']:
        loglevel = 'debug'
    else:
        loglevel = 'info'

    ensembler.refinement.solvate_models(
        process_only_these_targets=targets,
        process_only_these_templates=templates,
        template_seqid_cutoff=template_seqid_cutoff,
        padding=padding,
        ff=args['--ff'],
        water_model=args['--water_model'],
        verbose=args['--verbose'],
    )
    ensembler.refinement.determine_nwaters(
        process_only_these_targets=targets,
        process_only_these_templates=templates,
        template_seqid_cutoff=template_seqid_cutoff,
        verbose=args['--verbose']
    )