helpstring_header = """\
Ensembler
"""

ensembler_helpstring = """\
Usage:
  ensembler -h | --help
  ensembler init [-h | --help] [--project_dir <dir>]
  ensembler gather_targets [-h | --help] [--gather_from <method>] [--query <query>] [--dbapi_uri <uri>] [--uniprot_domain_regex <regex>]
  ensembler gather_templates [-h | --help] [--gather_from <method>] [--query <query>] [--dbapi_uri <uri>] [--uniprot_domain_regex <regex>] [--structure_path <path> ...] [--no-loopmodel] [--overwrite_structures]
  ensembler build_models [-h | --help] [--targets <target>...] [--templates <template>...] [-v | --verbose]
  ensembler sort_and_cluster [-h | --help] [--targets <target>...] [-v | --verbose]
  ensembler refine_implicit [-h | --help] [--targets <target>...] [--targetsfile <targetsfile>] [--templates <template>...] [--gpupn <gpupn>] [--openmm_platform <platform>] [-v | --verbose]
  ensembler refine_explicit [-h | --help] [--targets <target>...] [--templates <template>...] [--gpupn <gpupn>] [--openmm_platform <platform>] [-v | --verbose]
  ensembler package_models [-h | --help] [--targets <target>...] [--package_for <choice>...] [--nFAHclones <n>] [--no-archiveFAHproject]

Commands:
  init                          Initialize a new Ensembler project
  gather_targets                Gather targets

General options:
  -h --help                     Print command line help
"""

gather_from_helpstring = """\
  --gather_from <method>            Choose a source from which to gather data. {targetexplorer|uniprot} [default: targetexplorer]
                                      - "targetexplorer": An existing TargetExplorer database, specified via the --dbapi_uri argument.

                                      - "uniprot": UniProt (www.uniprot.org). Requires a user-defined query
                                      string, plus an optional regex string for further filtering.
"""