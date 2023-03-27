import pandas as pd
import os

class dataloader():
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
                    'NAtoms' : self.data['NAtoms'][i],
                    'Ne' : self.data['Ne'][i],
                    'Homo' : self.data['Homo'][i],
                    'Lumo' : self.data['Lumo'][i],
                    'CIe' : self.data['CIe'][i],
                    'HF' : self.data['HF'][i],
                    'Factor': self.data['Factor'][i],
                    'Energies' : [energie for energie in self.energies[ID].dropna()]
                }

            else:
                values[ID] = {
                    'NAtoms' : self.data['NAtoms'][i],
                    'Ne' : self.data['Ne'][i],
                    'Homo' : self.data['Homo'][i],
                    'Lumo' : self.data['Lumo'][i],
                    'CIe' : self.data['CIe'][i],
                    'HF' : self.data['HF'][i],
                    'Factor': self.data['Factor'][i],
                    'Theta' : self.data['Theta'][i],
                    'Mu' : self.data['Mu'][i],
                    'B_opt1' : self.data['B_opt1'][i],
                    'B_opt2' : self.data['B_opt2'][i],
                    'B_opt3' : self.data['B_opt3'][i],
                    'B_opt4' : self.data['B_opt4'][i],
                    'B_opt5' : self.data['B_opt5'][i],
                    'Energies' : [energie for energie in self.energies[ID].dropna()]
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