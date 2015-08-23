from ensembler.tools.renumber_residues import RenumberResidues

helpstring_header = """Renumber residues using the canonical UniProt sequence coordinates.
Target IDs must start with the UniProt mnemonic, e.g. 'ABL1_HUMAN'
"""

helpstring_unique_options = [
    """\
  --target <targetid>          ID for target to work, e.g. 'ABL1_HUMAN_D0'""",
]

helpstring_nonunique_options = [
    """\
  -v --verbose                 """,
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_unique_options), '\n\n'.join(helpstring_nonunique_options)])
docopt_helpstring = '\n\n'.join(helpstring_unique_options)


def dispatch(args):
    if args['--verbose']:
        log_level = 'debug'
    else:
        log_level = 'info'

    RenumberResidues(
        targetid=args['--target'],
        log_level=log_level
    )
