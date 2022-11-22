import re
import os


alleles = {'ref': 'alt', 'alt': 'ref'}
starting_columns = ['#chr', 'start', 'end', 'ID', 'ref', 'alt']

def make_bad_dir(fit_dir, bad):
    BAD = convert_frac_to_float(bad)
    bad_dir = os.path.join(fit_dir, 'BAD{:.2f}'.format(BAD))
    if not os.path.exists(bad_dir):
        os.mkdir(bad_dir)
    return BAD, bad_dir


def convert_frac_to_float(string):
    if re.match(r"^[1-9]+[0-9]*/[1-9]+[0-9]*$", string):
        num, denom = string.split('/')
        if int(denom) <= 0:
            return False
        else:
            value = int(num) / int(denom)
    elif re.match(r"^[1-9]+[0-9]*\.[1-9]+[0-9]*$", string):
        try:
            value = float(string)
        except ValueError:
            return False
    elif re.match(r"^[1-9]+[0-9]*$", string):
        try:
            value = int(string)
        except ValueError:
            return False
    else:
        return False
    if value >= 1:
        return value
    else:
        return False


def check_states(string):
    if not string:
        raise ValueError
    string = string.strip().split(',')
    ret_val = list(map(convert_frac_to_float, string))
    if not all(ret_val):
        raise ValueError
    else:
        return ret_val

def get_field_by_ftype(allele, ftype='pval'):
    if ftype == 'pval':
        return f'pval_{allele}'
    elif ftype == 'es':
        return f'es_{allele}'
    elif ftype == 'w':
        return f'w_{allele}'
    elif ftype == 'pval-ag':
        return f'logit_pval_{allele}'
    elif ftype == 'es-weighted-mean':
        return f'aggregated_es_weighted_{allele}'
    elif ftype == 'es-mean':
        return f'aggregated_es_{allele}'
    else:
        raise ValueError