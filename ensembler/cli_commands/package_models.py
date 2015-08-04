import ensembler
import ensembler.packaging

helpstring_header = """\
Package models for transfer or for set-up as a Folding@Home project.

Packaged models are written in directories of the form "packaged_models/fah-projects/[target id]".

MPI-enabled.

Options."""

helpstring_unique_options = [
    """\
  --package_for <choice>                Specify which packaging method to use (required).
                                        - transfer: compress results into a single .tgz file
                                        - FAH: set up the input files and directory structure
                                          necessary to start a Folding@Home project.""",

    """\
  --nfahclones <n>                      If packaging for Folding@Home, select the number of clones
                                        to use for each model [default: 1].""",

    """\
  --compressruns                        If packaging for Folding@Home, choose whether to compress
                                        each RUN into a .tgz file.""",
]

helpstring_nonunique_options = [
    """\
  --targetsfile <targetsfile>  File containing a list of target IDs to work on (newline-separated).
                               Comment targets out with "#".""",

    """\
  --targets <target>           Define one or more target IDs to work on (comma-separated), e.g.
                               "--targets ABL1_HUMAN_D0,SRC_HUMAN_D0" (default: all targets)""",

    """\
  --templates <template>            Define one or more template IDs to work on (comma-separated), e.g.
                                    "--templates ABL1_HUMAN_D0_1OPL_A" (default: all templates)""",

    """\
  --templatesfile <templatesfile>   File containing a list of template IDs to work on (newline-separated).
                                    Comment targets out with "#".""",

    """\
  --template_seqid_cutoff <cutoff>  Select only templates with sequence identity (percentage)
                                    greater than the given cutoff.""",

    """\
  -v --verbose                 """,
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_unique_options), '\n\n'.join(helpstring_nonunique_options)])
docopt_helpstring = '\n\n'.join(helpstring_unique_options)

def dispatch(args):
    if args['--package_for']:
        package_for = args['--package_for']
    else:
        raise KeyError('--package_for choice must be defined')

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

    if args['--template_seqid_cutoff']:
        template_seqid_cutoff = float(args['--template_seqid_cutoff'])
    else:
        template_seqid_cutoff = False

    if args['--nfahclones']:
        n_fah_clones = int(args['--nfahclones'])
    else:
        n_fah_clones = 1

    if args['--compressruns']:
        archive = True
    else:
        archive = False

    if args['--verbose']:
        loglevel = 'debug'
    else:
        loglevel = 'info'

    if package_for.lower() == 'transfer':
        ensembler.packaging.package_for_transfer(
            process_only_these_targets=targets,
            process_only_these_templates=templates,
        )

    elif package_for.lower() == 'fah':
        ensembler.packaging.package_for_fah(
            process_only_these_targets=targets,
            process_only_these_templates=templates,
            template_seqid_cutoff=template_seqid_cutoff,
            nclones=n_fah_clones,
            archive=archive,
            loglevel=loglevel,
        )