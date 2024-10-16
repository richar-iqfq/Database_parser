import re

energie_regex = re.compile(r'\s+([+-]?\d+.\d*)')
alpha_homo_regex = re.compile(r'Alpha\s+occ.')
alpha_lumo_regex = re.compile(r'Alpha\s+virt.')
beta_homo_regex = re.compile(r'Beta\s+occ.')
beta_lumo_regex = re.compile(r'Beta\s+virt.')
alpha_electrons_regex = re.compile(r'\s+(\d+)\s+alpha\s+electrons')
beta_electrons_regex = re.compile(r'\s+(\d+)\s+beta\s+electrons')
molar_volume_regex = re.compile(r'Molar volume')
molar_volume_value_regex = re.compile(r'(\d+.\d*)\s+cm\*\*3\/mol')
point_group_line_regex = re.compile(r'Full point group')
point_group_value_regex = re.compile(r"group\s+([A-Z]+\*?\d*\w*)\s+")
dipole_line_regex = re.compile(r'^\s+Dipole\s+moment')
quadrupole_line_regex = re.compile(r'^\s+Quadrupole\s+moment')
traceless_line_regex = re.compile(r'^\s+Traceless\s+Quadrupole')
coordinates_values_regex = (r'\s+([XYZ]{1,2})=\s+([-+]?\d+.\d*)')
equals_to_pattern = r'\s+{}=\s*([+-]?[0-9.]+)'