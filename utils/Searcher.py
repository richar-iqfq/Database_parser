import re
import os
import pandas as pd

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
        self.__Internal_file_selected = None
        self.energies = {}
        self.values = {}
        self.Properties = ['ID', 'NAtoms', 'Ne', 'Homo-4', 'Homo-3',
                            'Homo-2', 'Homo-1', 'Homo', 'Lumo', 
                            'Lumo+1', 'Lumo+2', 'Lumo+3', 'Lumo+4',
                            'CIe', 'HF', 'ET', 'EV', 'EJ', 'EK', 
                            'ENuc', 'Dipole', 'Factor']
        
    def __orbitals(self):
        file = self.__Internal_file_selected
        orbital_energies = {'a_homo' : [],
                             'a_lumo' : [],
                             'b_homo' : [], 
                             'b_lumo' : []}

        REa_homo = r'Alpha\s+occ.'
        REa_lumo = r'Alpha\s+virt.'
        REb_homo = r'Beta\s+occ.'
        REb_lumo = r'Beta\s+virt.'

        RE_values = re.compile(r'\s+([+-]?\d+.\d*)')
        
        is_betha = False

        Ea_lines_homo = []
        Ea_lines_lumo = []
        
        Eb_lines_homo = []
        Eb_lines_lumo = []

        with open(file, 'r') as log:
            for line in log:
                if re.search(REb_homo, line):
                    is_betha = True
                    break
        
        if is_betha:
            factor = 1
        else:
            factor = 2

        with open(file, 'r') as log:
            if is_betha:
                for line in log:
                    matcha_homo = re.search(REa_homo, line)
                    matcha_lumo = re.search(REa_lumo, line)
                    
                    matchb_homo = re.search(REb_homo, line)
                    matchb_lumo = re.search(REb_lumo, line)

                    if matcha_homo != None:
                        Ea_lines_homo.append(line)
                    if matcha_lumo != None:
                        Ea_lines_lumo.append(line)
                    
                    if matchb_homo != None:
                        Eb_lines_homo.append(line)
                    if matchb_lumo != None:
                        Eb_lines_lumo.append(line)
            else:
                for line in log:
                    matcha_homo = re.search(REa_homo, line)
                    matcha_lumo = re.search(REa_lumo, line)
                    if matcha_homo != None:
                        Ea_lines_homo.append(line)
                    if matcha_lumo != None:
                        Ea_lines_lumo.append(line)
        
        # Extract all the orbital energies
        energies = []
        for line in Ea_lines_homo:
            e = RE_values.findall(line)
            orbital_energies['a_homo'] += e

        for line in Ea_lines_lumo:
            e = RE_values.findall(line)
            orbital_energies['a_lumo'] += e

        if is_betha:
            for line in Eb_lines_homo:
                e = RE_values.findall(line)
                orbital_energies['b_homo'] += e

            for line in Eb_lines_lumo:
                e = RE_values.findall(line)
                orbital_energies['b_lumo'] += e
        
        for list in [orbital_energies[key] for key in orbital_energies.keys()]:
            energies += list

        if is_betha:
            orbital_homo = orbital_energies['a_homo'] + orbital_energies['b_homo']
            orbital_homo = [float(val) for val in orbital_homo]
            orbital_homo = sorted(orbital_homo)

            orbital_lumo = orbital_energies['a_lumo'] + orbital_energies['b_lumo']
            orbital_lumo = [float(val) for val in orbital_lumo]
            orbital_lumo = sorted(orbital_lumo)
            
            try:
                E_homo4 = orbital_homo[-5]
            except:
                E_homo4 = 0
            
            try:
                E_homo3 = orbital_homo[-4]
            except:
                E_homo3 = 0

            try:
                E_homo2 = orbital_homo[-3]
            except:
                E_homo2 = 0
            
            E_homo1 = orbital_homo[-2]
            E_homo = orbital_homo[-1]
            E_lumo = orbital_lumo[0]
            E_lumo1 = orbital_lumo[1]

            try:
                E_lumo2 = orbital_lumo[2]
            except:
                E_lumo2 = 0

            try:
                E_lumo3 = orbital_lumo[3]
            except:
                E_lumo3 = 0
            
            try:
                E_lumo4 = orbital_lumo[4]
            except:
                E_lumo4 = 0

        else:
            # Assign the homo and lumo value
            try:
                E_homo4 = orbital_energies['a_homo'][-5]
            except:
                E_homo4 = 0

            try:
                E_homo3 = orbital_energies['a_homo'][-4]
            except:
                E_homo3 = 0
            
            try:
                E_homo2 = orbital_energies['a_homo'][-3]
            except:
                E_homo2 = 0
            
            E_homo1 = orbital_energies['a_homo'][-2]
            E_homo = orbital_energies['a_homo'][-1]
            E_lumo = orbital_energies['a_lumo'][0]
            E_lumo1 = orbital_energies['a_lumo'][1]
            
            try:
                E_lumo2 = orbital_energies['a_lumo'][2]
            except:
                E_lumo2 = 0
            
            try:
                E_lumo3 = orbital_energies['a_lumo'][3]
            except:
                E_lumo3 = 0
            
            try:
                E_lumo4 = orbital_energies['a_lumo'][4]
            except:
                E_lumo4 = 0
        
        return (E_homo4, E_homo3, E_homo2, E_homo1, E_homo, E_lumo, E_lumo1, E_lumo2, E_lumo3, E_lumo4, energies, factor)
    
    def __Equals2(self, keyword):
        file = self.__Internal_file_selected

        RE_equalsto = r'\s+{}=\s*([+-]?[0-9.]+)'.format(keyword)
        
        if keyword == 'DE\(Corr\)':
            matches = []
            with open(file, 'r') as log:
                for line in log:
                    match = re.search(RE_equalsto, line)
                    if match != None:
                        matches.append(line)
            value = re.search(RE_equalsto, matches[-1])
            return value.group(1)
        else:
            with open(file, 'r') as log:
                for line in log:
                    value = re.search(RE_equalsto, line)
                    if value != None:
                        return value.group(1)

    def __electron_number(self):
        file = self.__Internal_file_selected

        RE_ae = r'\s+(\d+)\s+alpha\s+electrons'
        RE_be = r'\s+(\d+)\s+beta\s+electrons'

        with open(file, 'r') as log:
            for line in log:
                matcha = re.search(RE_ae, line)
                matchb = re.search(RE_be, line)

                if matcha != None:
                    ae = matcha.group(1)
                if matchb != None:
                    be = matchb.group(1)

        return int(ae) + int(be)

    def __build_structure(self):
        for key in self.Properties:
            self.values[key] = []

    def search(self, get=False):
        '''
        Compute data extraction for the molecules inside the input_path folder,
        if get is True, returns all the found values.

        Parameters
        ----------
        get (`bool`)
            if get return self.values, default is False
        '''
        self.__build_structure()

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
        for fileHF in HFt:
            for fileCI in CIt:
                if fileCI[fileCI.rfind('/')+1::] in fileHF:
                    CI.append(fileCI)
                    HF.append(fileHF)

        for file in HF:
            self.__Internal_file_selected = file
            # print(file)

            # Identifiers
            index = file.rfind('/') + 1
            t_file = file[index::]
            ID = re.search(r'[/\\]?([A-Z0-9]+)[_a-z]*.log', t_file).group(1)
            self.values['ID'].append(ID)

            # NAtoms
            NAtoms = self.__Equals2('NAtoms')
            self.values['NAtoms'].append(NAtoms)

            # Ne
            Ne = self.__electron_number()
            self.values['Ne'].append(Ne)

            # Homo, Lumo, orbital_energies
            (Homo4, Homo3, Homo2, Homo1, Homo, Lumo, Lumo1, Lumo2, Lumo3, Lumo4, energies, factor) = self.__orbitals()
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
            HFe = self.__Equals2('ETot')
            if HFe == None:
                HFe = self.__Equals2('HF')
                self.values['HF'].append(HFe)
            else:
                self.values['HF'].append(HFe)

            # ET
            ET = self.__Equals2('ET')
            self.values['ET'].append(float(ET)/float(Ne))

            # EV
            EV = self.__Equals2('EV')
            if EV != None:
                self.values['EV'].append(float(EV)/float(Ne))
            else:
                self.values['EV'].append('NaN')
            # EJ
            EJ = self.__Equals2('EJ')
            self.values['EJ'].append(float(EJ)/float(Ne))

            # EK
            EK = self.__Equals2('EK')
            self.values['EK'].append(float(EK)/float(Ne))

            # ENuc
            ENuc = self.__Equals2('ENuc')
            self.values['ENuc'].append(float(ENuc)/float(Ne))

            # Dipole moment
            Dm = self.__Equals2('Tot')
            self.values['Dipole'].append(Dm)

        for file in CI:
            self.__Internal_file_selected = file
            CIE = self.__Equals2('DE\(Corr\)')

            self.values['CIe'].append(CIE)

        if get:
            return self.values

    def Save(self):
        '''
        Save the info into csv files to the output_path:  data.csv and energies.csv
        '''
        main_path = os.path.join(os.getcwd(), self.input_path)
        exit_path = os.path.join(main_path, self.output_path)
        
        if not os.path.isdir(exit_path):
            os.mkdir(exit_path)
        
        data = pd.DataFrame(self.values)
        data = data.sort_values(by=['ID'])

        energies_DF = pd.DataFrame(self.energies)

        data.to_csv(exit_path + '/data.csv', index=False)
        energies_DF.to_csv(exit_path + '/energies.csv', index=False)
        print('Exported')