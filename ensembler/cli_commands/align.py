import ensembler
import ensembler.modeling

helpstring_header = """\
Pairwise alignment of target sequences onto template sequences.

These alignments are used to guide the subsequent modeling step, and are stored as follows:
"models/[target id]/[template id]/alignment.pir"
The ".pir" alignment format is an ascii-based format required by Modeller.

If the loopmodel function was used previously, then templates which have been successfully remodeled
will be selected for this alignment and the subsequent modeling and refinement steps. Otherwise,
Ensembler defaults to using the template structures which contain only resolved residues.

MPI-enabled.

Options:"""

helpstring_unique_options = [
    """\
  --targetsfile <targetsfile>       File containing a list of target IDs to work on (newline-separated).
                                    Comment targets out with "#".""",

    """\
  --targets <target>                Define one or more target IDs to work on (comma-separated), e.g.
                                    "--targets ABL1_HUMAN_D0,SRC_HUMAN_D0" (default: all targets)""",
]

helpstring_nonunique_options = [
    """\
  --templates <template>            Define one or more template IDs to work on (comma-separated), e.g.
                                    "--templates ABL1_HUMAN_D0_1OPL_A" (default: all templates)""",

    """\
  --templatesfile <templatesfile>   File containing a list of template IDs to work on (newline-separated).
                                    Comment templates out with "#".""",

    """\
  -v --verbose                      """,
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_unique_options), '\n\n'.join(helpstring_nonunique_options)])
docopt_helpstring = '\n\n'.join(helpstring_unique_options)


def dispatch(args):
    if args['--verbose']:
        loglevel = 'debug'
    else:
        loglevel = 'info'

    if args['--targetsfile']:
        with open(args['--targetsfile'], 'r') as targetsfile:
            targets = [line.strip() for line in targetsfile.readlines() if line[0] != '#']
    elif args['--targets']:
        targets = args['--targets'].split(',')
    else:
        targets = False

    if args['--templatesfile']:
        with open(args['--templatesfile'], 'r') as templatesfile:
            templates = [line.strip() for line in templatesfile.readlines() if line[0] != '#']
    elif args['--templates']:
        templates = args['--templates'].split(',')
    else:
        templates = False

    ensembler.modeling.align_targets_and_templates(process_only_these_targets=targets, process_only_these_templates=templates, loglevel=loglevel)