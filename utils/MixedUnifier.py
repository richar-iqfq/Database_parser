import os
import pandas as pd

class MixedDataUnifier ():
    def __init__(self):
        self.database = 'Database_##.csv'
        self.energies = 'Database_Energies_##.feather'
        self.smiles = 'Database_SMILES_##.csv'

        self.names = ['data', 'smile', 'energie']

    def unify(self, save=True, sets = ['anion', 'cation', 'radical']):
        print('Appending data...')

        for name in self.names:

            if name == 'data':
                main_df = pd.read_csv(self.database.replace('##', 'neutro'))
                main_df = self.include_columns(main_df, 'neutro')

                for type in sets:
                    df = self.rename_id(name, type)
                    df = self.include_columns(df, type)
                    main_df = pd.concat([main_df, df], ignore_index=True)

                if save:
                    main_df.to_csv('Database_Full.csv', index=False)

            if name == 'smile':
                main_df = pd.read_csv(self.smiles.replace('##', 'neutro'))

                for type in sets:
                    df = self.rename_id(name, type)
                    main_df = pd.concat([main_df, df], ignore_index=True)

                if save:
                    main_df.to_csv('Database_SMILES_Full.csv', index=False)

            if name == 'energie':
                main_df = pd.read_feather(self.energies.replace('##', 'neutro'))

                for type in sets:
                    df = self.rename_id(name, type)
                    main_df = pd.concat([main_df, df], axis=1)

                if save:
                    main_df = main_df.loc[:,~main_df.columns.duplicated()].copy()

                    main_df.to_feather('Database_Energies_Full.feather')

    def rename_id(self, selection, type):
        keys = {
            'anion' : 'a',
            'cation' : 'c',
            'radical' : 'r'
        }

        if selection == 'data':
            # For data
            df = pd.read_csv(self.database.replace('##', type))
            df['ID'] = df['ID'] + keys[type]

        elif selection == 'smile':
            # For smiles
            df = pd.read_csv(self.smiles.replace('##', type))
            df['ID'] = df['ID'] + keys[type]

        elif selection == 'energie':
            # For energies
            df = pd.read_feather(self.energies.replace('##', type))
            df.columns = df.columns + keys[type]

        return df
    
    def include_columns(self, DataFrame, type):
        charges = {
            'anion' : -1,
            'cation' : 1,
            'radical' : 0,
            'neutro' : 0
        }

        multiplicity = {
            'anion' : 2,
            'cation' : 2,
            'radical' : 3,
            'neutro' : 1
        }

        lenght = len(DataFrame)

        charge = [charges[type]] * lenght
        multiplicity = [multiplicity[type]] * lenght

        DataFrame['charge'] = charge
        DataFrame['multiplicity'] = multiplicity

        return DataFrame
