import ensembler.initproject

helpstring = """\
Initialize Ensembler project by creating necessary subdirectories and a project metadata .yaml file.

Options:
  --project_dir <dir>           Directory within which to initialize new project [default: .]
"""
docopt_helpstring = helpstring

required_args = []


def dispatch(args):
    ensembler.cli.validate_args(args, required_args)
    ensembler.initproject.InitProject(args['--project_dir'])