import pandas as pd
import os

class Dataloader():
    '''
    Dataloader for compute calculations

    Attributes
    ----------
    database (`str`):
        Csv database with all data

    energies (`str`):
        Feather type file with all the molecule's energies.

    Methods
    -------
    load(raw, amount='All'):
        Load the data inside the database and energies and returns a dictionary with
        the info. 

        raw (`bool`):
            True if database do not have optimized values

        amount ('All' or `int`):
            Number of molecules to return in dict, default is 'All'
    '''
    def __init__(self, database, energies):
        self.data = pd.read_csv(database)
        self.energies = pd.read_feather(energies)
        self.root = os.getcwd()

    def load(self, raw, amount='All'):
        '''
        load the data inside the database and energies given to class constructor.

        Parameters
        ----------
        raw (`bool`):
            True for raw dataset (without optimized values like b1 or theta, mu)
        
        amount (`str`) or (`int`):
            refers to the number of molecules to return in dict, default
            is 'All'

        Returns
        -------
        Dictionary with the info.
        '''
        if amount == 'All':
            self.amount = 'All'
            mols = self.data['ID']
        else:
            self.amount = amount
            mols = self.data['ID'][0:amount]

        values = {}
        for i, ID in enumerate(mols):
            if raw:
                values[ID] = {
                    'NAtoms' : int(self.data['NAtoms'][i]),
                    'Ne' : int(self.data['Ne'][i]),
                    'Homo' : float(self.data['Homo'][i]),
                    'Lumo' : float(self.data['Lumo'][i]),
                    'CIe' : float(self.data['CIe'][i]),
                    'HF' : float(self.data['HF'][i]),
                    'Factor': int(self.data['Factor'][i]),
                    'Energies' : [float(energie) for energie in self.energies[ID].dropna()]
                }

            else:
                values[ID] = {
                    'NAtoms' : int(self.data['NAtoms'][i]),
                    'Ne' : int(self.data['Ne'][i]),
                    'Homo' : float(self.data['Homo'][i]),
                    'Lumo' : float(self.data['Lumo'][i]),
                    'CIe' : float(self.data['CIe'][i]),
                    'HF' : float(self.data['HF'][i]),
                    'Factor': int(self.data['Factor'][i]),
                    'Theta' : float(self.data['Theta'][i]),
                    'Mu' : float(self.data['Mu'][i]),
                    'B_opt1' : float(self.data['B_opt1'][i]),
                    'B_opt2' : float(self.data['B_opt2'][i]),
                    'B_opt3' : float(self.data['B_opt3'][i]),
                    'B_opt4' : float(self.data['B_opt4'][i]),
                    'B_opt5' : float(self.data['B_opt5'][i]),
                    'B_opt6' : float(self.data['B_opt6'][i]),
                    'B_opt7' : float(self.data['B_opt7'][i]),
                    'Energies' : [float(energie) for energie in self.energies[ID].dropna()]
                }
        
        return values

    def get_IDS(self):
        if not self.amount:
            print('First run load method!')        
        elif self.amount == 'All':
            IDS = self.data['ID']
        else:
            IDS = self.data['ID'][0:self.amount]

        return IDS