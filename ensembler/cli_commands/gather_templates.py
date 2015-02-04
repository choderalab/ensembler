import ast
import ensembler
import ensembler.initproject

helpstring_header = """\
Gather template protein data from a specified resource, such as UniProt or a TargetExplorer
database.

Template sequences are written in two forms, to "templates/templates-resolved-seq.fa" and
"templates/templates-full-seq.fa". The former contains only residues resolved in the structure; the
latter contains the complete UniProt sequence contained within the span of the structure, including
unresolved residues. Template structures (containing only resolved residues) are written to
"templates/structures-resolved". Template structures containing the full sequences can optionally be
generated with a subsequent step - the loopmodel function.

Templates are given IDs of the form [UniProt mnemonic]_D[domain id]_[PDB id]_[chain id], where the
final two elements represent the PBD ID and chain identifier. Example: EGFR_HUMAN_D0_2GS7_B.

MPI-enabled.

Options:"""

gather_from_options = ['targetexplorer', 'uniprot', 'pdb']

helpstring_nonunique_options = [
    """\
  --gather_from <method>         Choose a source from which to gather data.
                                 {targetexplorer|uniprot|pdb} [default: targetexplorer]
                                   - "uniprot": UniProt (www.uniprot.org). Requires a query string
                                     defined by the --query flag, plus an optional regular
                                     expression for selecting domains (--uniprot_domain_regex).
                                   - "targetexplorer": a TargetExplorer database
                                     (https://github.com/choderalab/targetexplorer). Requires a
                                     query string defined by the --query flag, and a URI specified
                                     by the --dbapi_uri flag.
                                   - "pdb": PDB (www.pdb.org). Requires a query string defined by
                                   the --query falg, plus an optional regular expression for
                                   selecting domains (--uniprot_domain_regex).""",

    """\
  --query <query>                Query string for selecting templates
                                   if --gather_from="uniprot":
                                     Use the same syntax as used on the UniProt website. Note that
                                     *all* domains contained within the returned UniProt entries
                                     will be selected as templates, unless the
                                     --uniprot_domain_regex option is used to select a subset. The
                                     script will print some information on the set of unique domain
                                     names returned by the initial UniProt search, which can help
                                     with constructing a suitable regex.

                                     Example: 'domain:"Protein kinase" AND taxonomy:9606 AND
                                     reviewed:yes' - this will return reviewed UniProt entries for
                                     human (taxonomy ID: 9606) proteins containing "Protein kinase"
                                     domain annotations.

                                   if --gather_from="targetexplorer":
                                     Use syntax as described in the TargetExplorer documentation.
                                   if --gather_from="pdb":
                                     List of comma-separated PDB IDs.""",

    """\
  --dbapi_uri <uri>              TargetExplorer database API URI, e.g.
                                 "http://plfah2.mskcc.org/kinomeDBAPI\"""",

    """\
  --uniprot_domain_regex <regex> Optional regular expression (case-sensitive) for selecting domains
                                 contained within UniProt entries. If not provided, all domains
                                 contained within returned UniProt entries will be selected as
                                 templates.

                                 Example: '^Protein kinase(?!; truncated)(?!; inactive)' - matches
                                 "Protein kinase", "Protein kinase; 1" and "Protein kinase; 2",
                                 and excludes "Protein kinase; truncated" and "Protein kinase;
                                 inactive\"""",
]

helpstring_unique_options = [
    """\
  --structure_paths <path>       Local directories within which to search for PDB and SIFTS files
                                 (comma-separated)""",

    """\
  --chainids <chainids>          if --gather_from="pdb":
                                   Optionally specify which PDB chain IDs to parse. Use a Python
                                   list of lists syntax (one list for each PDB ID), e.g.
                                   '[["A", "D"], ["A", "B", "C"]]'""",
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_nonunique_options), '\n\n'.join(helpstring_unique_options)])
docopt_helpstring = '\n\n'.join(helpstring_unique_options)


def dispatch(args):
    if args['--structure_paths']:
        structure_paths = args['--structure_paths'].split(',')
    else:
        structure_paths = False

    if args['--gather_from'].lower() == 'targetexplorer':
        required_args = ['--dbapi_uri']
        ensembler.cli.validate_args(args, required_args)
        ensembler.initproject.gather_templates_from_targetexplorer(args['--dbapi_uri'], search_string=args['--query'], structure_dirs=structure_paths)

    elif args['--gather_from'].lower() == 'uniprot':
        required_args = ['--query']
        ensembler.cli.validate_args(args, required_args)
        ensembler.initproject.gather_templates_from_uniprot(args['--query'], uniprot_domain_regex=args['--uniprot_domain_regex'], structure_dirs=structure_paths)

    elif args['--gather_from'].lower() == 'pdb':
        required_args = ['--query', '--uniprot_domain_regex']
        ensembler.cli.validate_args(args, required_args)

        pdbids = args['--query'].split(',')

        if args['--chainids']:
            chainids_list = ast.literal_eval(args['--chainids'])
            if len(chainids_list) != len(pdbids):
                raise Exception('The number of lists passed to --chainids must be the same as the number of PDB IDs.')
            chainids = {pdbids[i]: chainids_list[i] for i in range(len(pdbids))}
        else:
            chainids = None

        ensembler.initproject.gather_templates_from_pdb(pdbids, uniprot_domain_regex=args['--uniprot_domain_regex'], chainids=chainids, structure_dirs=structure_paths)

    else:
        raise Exception('--gather_from flag must be set to any of %r' % gather_from_options)