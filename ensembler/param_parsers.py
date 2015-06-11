import ast
import re
import operator as op
import simtk.unit

# needs to be sorted in reverse order so that 'nanoseconds' is matched before 'nano'
unit_membernames = sorted([name for name in simtk.unit.__dict__], key=len, reverse=True)
safe_names = {'None': None, 'True': True, 'False': False}

quantity_without_operator_regex = re.compile(
    '([0-9.]+) ?({0})'.format('|'.join(unit_membernames))
)   # e.g. "2 picoseconds"

valid_operators = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul,
             ast.Div: op.truediv, ast.Pow: op.pow, ast.BitXor: op.xor,
             ast.USub: op.neg}


def parse_api_params_string(params_string):
    """
    Safely parse a string representing a dict of kwargs to be passed to an API function.

    Parameters
    ----------
    params_string: str

    Returns
    -------
    dict

    Examples
    --------
    >>> parse_api_params_string('{"a": 3 / picoseconds, "b": "x", "c": 2.4}')
    """

    expr = ast.parse(params_string, mode='eval')
    if not isinstance(expr.body, ast.Dict):
        raise TypeError('Parsed string - {0} - should return a dict'.format(params_string))
    return safe_eval(expr.body)


def safe_eval(node):
    if isinstance(node, ast.Dict):   # <dict>
        return {safe_eval(key): safe_eval(value) for key, value in zip(node.keys, node.values)}
    elif isinstance(node, ast.Num):   # <number>
        return node.n
    elif isinstance(node, ast.Str):   # <str>
        return node.s
    elif isinstance(node, ast.Name) and node.id in unit_membernames:   # <member of simtk.unit>
        return getattr(simtk.unit, node.id)
    elif isinstance(node, ast.Name) and node.id in safe_names:   # None, True, False
        return safe_names[node.id]
    elif isinstance(node, ast.BinOp):   # <left> <operator> <right>
        return valid_operators[type(node.op)](safe_eval(node.left), safe_eval(node.right))
    elif isinstance(node, ast.UnaryOp):   # <operator> <operand> e.g., -1
        return valid_operators[type(node.op)](safe_eval(node.operand))
    else:
        raise TypeError(node)


def eval_quantity_string(param_value_string):
    """
    Safely evaluate simtk quantities passed from CLI, using either Python expression syntax
    ('2 * picoseconds' or '2 / picoseconds') or a more natural syntax ('2 picoseconds').

    Parameters
    ----------
    param_value_string: str

    Examples
    --------
    >>> eval_quantity_string('2 picoseconds')
    >>> eval_quantity_string('2 / picoseconds')
    >>> eval_quantity_string('2')
    """

    quantity_as_number_space_unit_match = re.match(quantity_without_operator_regex, param_value_string)
    if quantity_as_number_space_unit_match:
        number, unit_name = quantity_as_number_space_unit_match.groups()
        number = ast.literal_eval(number)
        unit_obj = getattr(simtk.unit, unit_name)
        return number * unit_obj

    else:
        expr = ast.parse(param_value_string, mode='eval')
        return safe_eval(expr.body)