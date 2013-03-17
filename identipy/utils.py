import re
from pyteomics import mass

def neutral_masses(spectrum):
    if 'params' in spectrum:
        exp_mass = spectrum['params']['pepmass'][0]
        charge = spectrum['params'].get('charge')
    else:
        exp_mass = spectrum['base peak']
        charge = spectrum['charge']
    if isinstance(charge, str):
        states = re.findall(r'(\d+)(\+|-)?', charge)
        states = [int(num)*(-1 if sign == '-' else 1)
                for num, sign in states]
    else: states = [charge]
    states.sort()
    return zip((exp_mass*ch - ch*mass.nist_mass['H+'][0][0]
            for ch in states), states)


