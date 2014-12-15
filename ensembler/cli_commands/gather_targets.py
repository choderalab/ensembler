import ensembler
import ensembler.initproject

helpstring_header = """\
Gather target protein data.

Options:"""

helpstring_unique_options = [
    """\
  --gather_from <method>            Choose a source from which to gather data. {targetexplorer|uniprot} [default: targetexplorer]
                                      - "targetexplorer": An existing TargetExplorer database, specified via the --dbapi_uri argument.
                                      - "uniprot": UniProt (www.uniprot.org). Requires a user-defined query
                                      string, plus an optional regex string for further filtering.""",

    """\
  --query <query>                   Query string
                                      if --gather_from="uniprot":
                                        Specify a UniProt search string, using the same syntax as on the UniProt
                                        site (note: not all syntax may be supported, but most basic searches will
                                        work). *All* domains contained within the returned UniProt entries will be
                                        selected as targets, unless the --uniprot_domain_regex option is used to
                                        select a subset. The script will print some information on the set of unique
                                        domain names returned by the initial UniProt search, which can help with
                                        constructing a suitable string for --uniprot_domain_regex.

                                        Example: 'domain:"Protein kinase" AND taxonomy:9606 AND reviewed:yes' - this
                                        will return reviewed UniProt entries for human (taxonomy ID: 9606) proteins
                                        containing "Protein kinase" domain annotations. Note that all domains contained
                                        with those entries (including domains which are not "Protein kinase") will be
                                        selected as targets, unless the --uniprot_domain_regex flag is also set.""",

    """\
  --dbapi_uri <uri>                 TargetExplorer database API URI, e.g. "http://plfah2.mskcc.org/kinomeDBAPI\"""",

    """\
  --uniprot_domain_regex <regex>    Optional regular expression for subselecting domains from within UniProt
                                    entries (case-sensitive). If not provided, all domains contained within
                                    returned UniProt entries will be selected as targets (this will often not be
                                    the desired behavior).

                                    Example: '^Protein kinase(?!; truncated)(?!; inactive)' - matches "Protein
                                    kinase" as well as "Protein kinase; 1" and "Protein kinase; 2\"""",
]

helpstring = '\n\n'.join([helpstring_header, '\n\n'.join(helpstring_unique_options)])
docopt_helpstring = helpstring


def dispatch(args):
    if args['--help']:
        print helpstring
        return
    if args['--gather_from'] == 'targetexplorer':
        required_args = ['--dbapi_uri']
        ensembler.cli.validate_args(args, required_args)
        ensembler.initproject.GatherTargetsFromTargetExplorer(args['--dbapi_uri'], search_string=args['--query'])

    elif args['--gather_from']== 'uniprot':
        required_args = ['--query']
        ensembler.cli.validate_args(args, required_args)
        ensembler.initproject.GatherTargetsFromUniProt(args['--query'], uniprot_domain_regex=args['--uniprot_domain_regex'])