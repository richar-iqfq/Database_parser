import re
import os
import pandas as pd
from tqdm import tqdm
from .RegexConstants import *

'''
Extract the information from the log files in the input_path given to the class
'''

class MainLogReader():
    '''
    Attributes
    ----------
    input_path (`str`):
        Folder with the log files to calculate

    output_path (`str`):
        Folder where the results will be stored

    Methods
    -------
    search(get=False):
        Compute data extraction for the molecules inside the input_path folder,
        if get is True, returns all the found values.

    save():
        Save the info into csv files to the output_path:  data.csv and energies.csv
    '''

    def __init__(self, input_path, output_path):
        self.input_path = input_path
        self.output_path = output_path
        self.energies = {}
        self.values = {}
        self.Properties = [
            'ID', 'NAtoms', 'Ne', 'Homo-4', 'Homo-3',
            'Homo-2', 'Homo-1', 'Homo', 'Lumo', 
            'Lumo+1', 'Lumo+2', 'Lumo+3', 'Lumo+4',
            'CIe', 'HF', 'ET', 'EV', 'EJ', 'EK', 
            'ENuc', 'Dipole', 'Volume', 'Pgroup',
            'Factor', 'X', 'Y', 'Z', 'XX', 'YY', 'ZZ',
            'XY', 'XZ', 'YZ'
        ]
        self.sep = os.sep
                
    def get_by_index(self, list, index):
        try:
            return list[index]
        except:
            return 0

    def __orbitals(self, file_content):
        orbital_energies = {
            'a_homo' : [],
            'a_lumo' : [],
            'b_homo' : [], 
            'b_lumo' : []
        }

        is_betha = any(beta_homo_regex.search(line) for line in file_content)
        factor = 1 if is_betha else 2

        Ea_lines_homo = []
        Ea_lines_lumo = []
        
        Eb_lines_homo = []
        Eb_lines_lumo = []

        for line in file_content:
            matcha_homo = alpha_homo_regex.search(line)
            matcha_lumo = alpha_lumo_regex.search(line)
            
            Ea_lines_homo.append(line) if matcha_homo else None
            Ea_lines_lumo.append(line) if matcha_lumo else None
            
            if is_betha:
                matchb_homo = beta_homo_regex.search(line)
                matchb_lumo = beta_lumo_regex.search(line)

                Eb_lines_homo.append(line) if matchb_homo else None
                Eb_lines_lumo.append(line) if matchb_lumo else None
    
        # Extract all the orbital energies
        for line in Ea_lines_homo:
            e = energie_regex.findall(line)
            orbital_energies['a_homo'] += e

        for line in Ea_lines_lumo:
            e = energie_regex.findall(line)
            orbital_energies['a_lumo'] += e

        if is_betha:
            for line in Eb_lines_homo:
                e = energie_regex.findall(line)
                orbital_energies['b_homo'] += e

            for line in Eb_lines_lumo:
                e = energie_regex.findall(line)
                orbital_energies['b_lumo'] += e
        
        energies = []
        for list in [orbital_energies[key] for key in orbital_energies.keys()]:
            energies += list

        # Append both lists (alpha and beta)
        orbital_homo = orbital_energies['a_homo'] + orbital_energies['b_homo']
        orbital_homo = [float(val) for val in orbital_homo]

        orbital_lumo = orbital_energies['a_lumo'] + orbital_energies['b_lumo']
        orbital_lumo = [float(val) for val in orbital_lumo]

        if is_betha:
            orbital_homo = sorted(orbital_homo)
            orbital_lumo = sorted(orbital_lumo)

        # Assign the homo and lumo value
        E_homo4 = self.get_by_index(orbital_homo, -5)
        E_homo3 = self.get_by_index(orbital_homo, -4)
        E_homo2 = self.get_by_index(orbital_homo, -3)
        E_homo1 = self.get_by_index(orbital_homo, -2)

        E_homo = self.get_by_index(orbital_homo, -1)
        E_lumo = self.get_by_index(orbital_lumo, 0)

        E_lumo1 = self.get_by_index(orbital_lumo, 1)
        E_lumo2 = self.get_by_index(orbital_lumo, 2)
        E_lumo3 = self.get_by_index(orbital_lumo, 3)
        E_lumo4 = self.get_by_index(orbital_lumo, 4)
            
        return (E_homo4, E_homo3, E_homo2, E_homo1, E_homo, E_lumo, E_lumo1, E_lumo2, E_lumo3, E_lumo4, energies, factor)
    
    def __Equals2(self, keyword, file_content):

        RE_equalsto = equals_to_pattern.format(keyword)
        
        if keyword == 'DE\(Corr\)':
            matches = []
            for line in file_content:
                match = re.search(RE_equalsto, line)
                if match != None:
                    matches.append(line)
            
            value = re.search(RE_equalsto, matches[-1])
            return value.group(1)
        
        else:
            for line in file_content:
                value = re.search(RE_equalsto, line)
                if value != None:
                    return value.group(1)

    def __electron_number(self, file_content):
        for line in file_content:
            matcha = alpha_electrons_regex.search(line)
            matchb = beta_electrons_regex.search(line)

            ae = matcha.group(1) if matcha else 0
            be = matchb.group(1) if matchb else 0

            if matcha and matchb:
                break

        return int(ae) + int(be)
    
    def __Volume(self, file_content):
        result = ''
        volume = None

        for line in file_content:
            match = molar_volume_regex.search(line)

            if match:
                result = line
                break
        
        volume_match = molar_volume_value_regex.search(result)

        if volume_match:
            volume = volume_match.group(1)

        return volume
    
    def __point_group(self, file_content):
        result = ''
        point_group = None

        for line in file_content:
            match = point_group_line_regex.search(line)

            if match:
                result = line
                break
        
        match = point_group_value_regex.search(result)

        if match:
            point_group = match.group(1)
        
        return point_group

    def __coordinates(self, file_content):
        dip_index = None
        quad_index = None
        trace_index = None

        # Get the lines index
        for i, line in enumerate(file_content):
            if dipole_line_regex.findall(line):
                dip_index = i

            if quadrupole_line_regex.findall(line):
                quad_index = i

            if traceless_line_regex.findall(line):
                trace_index = i

            if dip_index and quad_index and trace_index:
                break
        
        dipole_values = file_content[dip_index+1:quad_index]
        quadrupole_values = file_content[quad_index+1:trace_index]

        dipole = []
        quadrupole = []

        for found in dipole_values:
            res = coordinates_values_regex.findall(found)
            dipole += res

        for found in quadrupole_values:
            res = coordinates_values_regex.findall(found)
            quadrupole += res

        return dipole, quadrupole

    def __build_structure(self):
        for key in self.Properties:
            self.values[key] = []

    def __get_files(self):
        main_path = os.path.join(os.getcwd(), self.input_path)
        CIt = [] # Provisional path
        HFt = [] # Provisional path

        for Dir, SubDir, Files in os.walk(main_path):
            if '80-e' in Dir or 'incomplete' in Dir:
                continue
            elif 'HF' in Dir:
                for file in Files:
                    if 'log' in file:
                        HFt.append(os.path.join(Dir, file))
            elif 'CI' in Dir:
                for file in Files:
                    if 'log' in file:
                        CIt.append(os.path.join(Dir, file))
        
        CI = []
        HF = []
        print(f"Preparing files...")
        for fileHF in set(HFt):
            for fileCI in set(CIt):
                if fileCI[fileCI.rfind(self.sep)+1::] in fileHF:
                    CI.append(fileCI)
                    HF.append(fileHF)
                    break

        return HF, CI
            
    def search(self, get=False):
        '''
        Compute data extraction for the molecules inside the input_path folder,
        if get is True, returns all the found values.

        Parameters
        ----------
        get (`bool`)
            if get return self.values, default is False
        '''
        print('Initializing...')
        self.__build_structure()

        HF, CI = self.__get_files()

        pbar_HF = tqdm(total=len(HF), desc='HF files', colour='cyan')

        try:
            for file in HF:
                with open(file, 'r') as log_file:
                    file_content = log_file.readlines()
                
                # Identifiers
                index = file.rfind(self.sep) + 1
                t_file = file[index::]

                ID = t_file.replace('.log', '')
                self.values['ID'].append(ID)

                # NAtoms
                NAtoms = self.__Equals2('NAtoms', file_content)
                self.values['NAtoms'].append(NAtoms)

                # Ne
                Ne = self.__electron_number(file_content)
                self.values['Ne'].append(Ne)

                # Homo, Lumo, orbital_energies
                (Homo4, Homo3, Homo2, Homo1, Homo, Lumo, Lumo1, Lumo2, Lumo3, Lumo4, energies, factor) = self.__orbitals(file_content)
                self.values['Homo-4'].append(Homo4)
                self.values['Homo-3'].append(Homo3)
                self.values['Homo-2'].append(Homo2)
                self.values['Homo-1'].append(Homo1)
                self.values['Homo'].append(Homo)
                self.values['Lumo'].append(Lumo)
                self.values['Lumo+1'].append(Lumo1)
                self.values['Lumo+2'].append(Lumo2)
                self.values['Lumo+3'].append(Lumo3)
                self.values['Lumo+4'].append(Lumo4)
                self.values['Factor'].append(factor)

                self.energies[ID] = pd.Series(energies)

                # HF
                HFe = self.__Equals2('ETot', file_content)
                if HFe == None:
                    HFe = self.__Equals2('HF', file_content)
                    self.values['HF'].append(HFe)
                else:
                    self.values['HF'].append(HFe)

                # ET
                ET = self.__Equals2('ET', file_content)
                if ET != None:
                    self.values['ET'].append(float(ET)/float(Ne))
                else:
                    self.values['ET'].append('NaN')

                # EV
                EV = self.__Equals2('EV', file_content)
                if EV != None:
                    self.values['EV'].append(float(EV)/float(Ne))
                else:
                    self.values['EV'].append('NaN')
                # EJ
                EJ = self.__Equals2('EJ', file_content)
                if EJ != None:
                    self.values['EJ'].append(float(EJ)/float(Ne))
                else:
                    self.values['EJ'].append('NaN')

                # EK
                EK = self.__Equals2('EK', file_content)
                if EK != None:
                    self.values['EK'].append(float(EK)/float(Ne))
                else:
                    self.values['EK'].append('NaN')

                # ENuc
                ENuc = self.__Equals2('ENuc', file_content)
                if ENuc != None:
                    self.values['ENuc'].append(float(ENuc)/float(Ne))
                else:
                    self.values['ENuc'].append('NaN')

                # Dipole moment
                Dm = self.__Equals2('Tot', file_content)
                if Dm != None:
                    self.values['Dipole'].append(Dm)
                else:
                    self.values['Dipole'].append('NaN')

                # Volume
                Vol = self.__Volume(file_content)
                if Vol != None:
                    self.values['Volume'].append(Vol)
                else:
                    self.values['Volume'].append('NaN')

                # Point group
                group = self.__point_group(file_content)
                if group != None:
                    self.values['Pgroup'].append(group)
                else:
                    self.values['Pgroup'].append('NaN')

                # Coordinates
                dipole, quadrupole = self.__coordinates(file_content)
                
                for coordinates in (dipole, quadrupole):
                    for value in coordinates:
                        name = value[0]
                        self.values[name].append(float(value[1]))
                
                pbar_HF.update()
        except:
            print(file)
            raise

        print('\n')
        pbar_CI = tqdm(total=len(CI), desc='CI files', colour='cyan')

        for file in CI:
            with open(file, 'r') as log_file:
                file_content = log_file.readlines()

            CIE = self.__Equals2('DE\(Corr\)', file_content)

            self.values['CIe'].append(CIE)

            pbar_CI.update()

        if get:
            return self.values

    def save(self):
        '''
        Save the info into csv files to the output_path:  data.csv and energies.csv
        '''
        if not os.path.isdir(self.output_path):
            os.mkdir(self.output_path)
        
        data = pd.DataFrame(self.values)
        data = data.sort_values(by=['ID'])

        energies_DF = pd.DataFrame(self.energies)

        data.to_csv(self.output_path + '/data.csv', index=False)
        energies_DF.to_feather(self.output_path + '/energies.feather')
        
        print('Exported')