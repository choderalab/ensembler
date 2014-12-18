from docopt import docopt
import ensembler.cli_commands

docopt_elements = [ensembler.cli_commands.general.helpstring_header, ensembler.cli_commands.general.ensembler_helpstring]
for command_str in ensembler.cli_commands.command_list:
    command = getattr(ensembler.cli_commands, command_str)
    docopt_elements.append(command.docopt_helpstring)

docopt_full_helpstring = '\n'.join(docopt_elements)


def validate_args(args, required_args):
    for required_arg in required_args:
        if args[required_arg] is None:
            raise Exception('Required argument %s is missing.' % required_arg)


def main():
    args = docopt(docopt_full_helpstring, help=False)

    # print args

    command_dispatched = False

    for command_str in ensembler.cli_commands.command_list:
        if args[command_str]:
            if args['--help']:
                command = getattr(ensembler.cli_commands, command_str)
                print command.helpstring
            else:
                if not args['init']:
                    ensembler.core.check_project_toplevel_dir()
                command = getattr(ensembler.cli_commands, command_str)
                command.dispatch(args)
            command_dispatched = True

    if not command_dispatched and args['--help']:
        print '\n'.join([ensembler.cli_commands.general.helpstring_header, ensembler.cli_commands.general.ensembler_helpstring])
        pass