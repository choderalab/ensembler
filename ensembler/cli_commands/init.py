import ensembler.utils
import ensembler.initproject

helpstring_header = """\
Initialize Ensembler project by creating necessary subdirectories and a project metadata .yaml file.

Options:"""

helpstring_unique_options = """\
  --project_dir <dir>           Directory within which to initialize new project (default: ".")"""

helpstring = '\n\n'.join([helpstring_header, helpstring_unique_options])
docopt_helpstring = helpstring_unique_options

required_args = []


def dispatch(args):
    project_dir = ensembler.utils.set_arg_with_default(args['--project_dir'], default_arg='.')
    ensembler.cli.validate_args(args, required_args)
    ensembler.initproject.InitProject(project_dir)