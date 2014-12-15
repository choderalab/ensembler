import ensembler
import ensembler.packaging

helpstring_header = """\
Packages models for transfer or for set-up as a Folding@Home project. Packaged models are stored in the packaged-models/ directory.

Options."""

helpstring_unique_options = [
    """\
  --package_for <choice>                Specify which packaging method to use.
                                        - transfer: compress results into a single .tgz file
                                        - FAH: set-up the input files and directory structure necessary to start a Folding@Home project.""",
    """\
  --nFAHclones <n>                      If packaging for Folding@Home, select the number of clones to use for each model [default: 10].""",
    """\
  --no-archiveFAHproject                If packaging for Folding@Home, choose whether to compress the results into a .tgz file.""",
]

helpstring_nonunique_options = [
    """\
  --targets <target>...           Define one or more target IDs to work on (e.g. "--targets ABL1_HUMAN_D0 --targets SRC_HUMAN_D0") (default: all targets)""",
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_unique_options), '\n\n'.join(helpstring_nonunique_options)])
docopt_helpstring = '\n\n'.join(helpstring_unique_options)

def dispatch(args):
    if args['--help']:
        print helpstring
        return

    if args['--package_for'] == 'transfer':
        ensembler.packaging.package_for_transfer(process_only_these_targets=args['--targets'])

    elif args['--package_for'] == 'FAH':
        ensembler.packaging.package_for_fah(process_only_these_targets=args['--targets'], nclones=int(args['--nFAHclones']), archive=args['--no-archiveFAHproject'])