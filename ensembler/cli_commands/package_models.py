import ensembler
import ensembler.packaging

helpstring_header = """\
Packages models for transfer or for set-up as a Folding@Home project. Packaged models are stored in the packaged-models/ directory.

Options."""

helpstring_unique_options = [
    """\
  --package_for <choice>                Specify which packaging method to use (required).
                                        - transfer: compress results into a single .tgz file
                                        - FAH: set-up the input files and directory structure necessary to start a Folding@Home project.""",
    """\
  --nfahclones <n>                      If packaging for Folding@Home, select the number of clones to use for each model [default: 1].""",
    """\
  --archivefahproject                If packaging for Folding@Home, choose whether to compress the results into a .tgz file.""",
]

helpstring_nonunique_options = [
    """\
  --targetsfile <targetsfile>           File containing a list of newline-separated target IDs to work on. Comment targets out with "#".""",
    """\
  --targets <target>                    Define one or more target IDs to work on (e.g. "--targets ABL1_HUMAN_D0 --targets SRC_HUMAN_D0") (default: all targets)""",
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

    if args['--nfahclones']:
        n_fah_clones = int(args['--nfahclones'])
    else:
        n_fah_clones = 1

    if args['--archivefahproject']:
        archive = True
    else:
        archive = False

    if package_for.lower() == 'transfer':
        ensembler.packaging.package_for_transfer(process_only_these_targets=targets)

    elif package_for.lower() == 'fah':
        ensembler.packaging.package_for_fah(process_only_these_targets=targets, nclones=n_fah_clones, archive=archive)