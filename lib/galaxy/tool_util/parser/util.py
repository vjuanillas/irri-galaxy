from collections import OrderedDict

DEFAULT_DELTA = 10000
DEFAULT_DELTA_FRAC = None

DEFAULT_METRIC = "mae"
DEFAULT_EPS = 0.01
DEFAULT_PIN_LABELS = None


def is_dict(item):
    return isinstance(item, dict) or isinstance(item, OrderedDict)


def _parse_name(name, argument):
    """Determine name of an input source from name and argument
    returns the name or if absent the argument property
    In the latter case, leading dashes are stripped and
    all remaining dashes are replaced by underscores.
    """
    if name is None:
        if argument is None:
            raise ValueError("parameter must specify a 'name' or 'argument'.")
        name = argument.lstrip("-").replace("-", "_")
    return name
