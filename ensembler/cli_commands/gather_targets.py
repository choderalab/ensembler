import ensembler
import ensembler.initproject

helpstring_header = """\
Gather target protein data from a specified resource, such as UniProt or a TargetExplorer database.

Sequences are written to "targets/targets.fa".

Targets are given IDs of the form [UniProt mnemonic]_D[domain id], which consist of the UniProt
name for the target and an identifier for the domain (since a single target protein contain
multiple domains of interest). Example: EGFR_HUMAN_D0.

If the specified resource is UniProt, it is recommended to use the --uniprot_domain_regex flag to
specify exactly which domains will be selected.

Runs serially.

Options:"""

helpstring_unique_options = [
    """\
  --gather_from <method>         Resource from which to gather data {uniprot|targetexplorer}
                                 [default: uniprot]
                                   - "uniprot": UniProt (www.uniprot.org). Requires a query string
                                     defined by the --query flag, plus an optional regular
                                     expression for selecting domains (--uniprot_domain_regex).
                                   - "targetexplorer": a TargetExplorer database
                                     (https://github.com/choderalab/targetexplorer). Requires a
                                     query string defined by the --query flag, and a URI specified
                                     by the --dbapi_uri flag.""",

    """\
  --query <query>                Query string for selecting targets
                                   if --gather_from="uniprot":
                                     Use the same syntax as used on the UniProt website. Note that
                                     *all* domains contained within the returned UniProt entries
                                     will be selected as targets, unless the --uniprot_domain_regex
                                     option is used to select a subset. The script will print some
                                     information on the set of unique domain names returned by the
                                     initial UniProt search, which can help with constructing a
                                     suitable regex.

                                     Example: 'domain:"Protein kinase" AND taxonomy:9606 AND
                                     reviewed:yes' - this will return reviewed UniProt entries for
                                     human (taxonomy ID: 9606) proteins containing "Protein kinase"
                                     domain annotations.

                                   if --gather_from="targetexplorer":
                                     Use syntax as described in the TargetExplorer documentation.""",

    """\
  --dbapi_uri <uri>              TargetExplorer database API URI, e.g.
                                 "http://plfah2.mskcc.org/kinomeDBAPI\"""",

    """\
  --uniprot_domain_regex <regex> Optional regular expression (case-sensitive) for selecting domains
                                 contained within UniProt entries. If not provided, all domains
                                 contained within returned UniProt entries will be selected as
                                 targets.

                                 Example: '^Protein kinase(?!; truncated)(?!; inactive)' - matches
                                 "Protein kinase", "Protein kinase; 1" and "Protein kinase; 2",
                                 and excludes "Protein kinase; truncated" and "Protein kinase;
                                 inactive\"""",

    """\
  -v --verbose                   """,
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_unique_options)])
docopt_helpstring = '\n\n'.join(helpstring_unique_options)


def dispatch(args):
    if args['--verbose']:
        loglevel = 'debug'
    else:
        loglevel = 'info'

    if args['--gather_from'].lower() == 'targetexplorer':
        required_args = ['--dbapi_uri']
        ensembler.cli.validate_args(args, required_args)
        ensembler.initproject.GatherTargetsFromTargetExplorer(args['--dbapi_uri'], search_string=args['--query'], loglevel=loglevel)

    elif args['--gather_from'].lower() == 'uniprot':
        required_args = ['--query']
        ensembler.cli.validate_args(args, required_args)
        ensembler.initproject.GatherTargetsFromUniProt(args['--query'], uniprot_domain_regex=args['--uniprot_domain_regex'], loglevel=loglevel)

    else:
        raise Exception('--gather_from flag must be set to either "uniprot" or "targetexplorer"')