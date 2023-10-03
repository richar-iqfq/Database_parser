from utils.Searcher import MainLogReader
from utils.Tools import raw_data_cleaner, raw_data_unify
from utils.MolBuilder import MolBuilder
import os
import re
import pandas as pd

class Extractor():
    '''
    Extract the data from the HF and CI files into csv contained inside output_path

    Attributes
    ----------
    output_path (`str`):
        output_path for the parsed data

    Methods
    -------
    extract(unify=True, output_keywords='AA')
        Extract the data from the HF and CI files

    unify (`bool`)
        if True, mix all results inside one csv file

    output_keywords (`str`)
        keywords to append at the result file, default is 'AA'
    '''
    def __init__(self, input_path, output_path):
        self.root = os.getcwd()
        
        for folder in os.listdir(input_path):
            if 'HF' in folder:
                self.HF_path = os.path.join(input_path, folder)
                break
        
        self.input_path = input_path
        self.output_path = os.path.join(input_path, output_path)

    def __save_smiles(self, mols_data, smiles, path):
        df = pd.DataFrame({
            'ID' : mols_data,
            'smiles' : smiles
        })

        df.to_csv(os.path.join(path, 'smiles.csv'), index=False)

    def __write_failed(self, failed, path):
        with open(os.path.join(path, 'failed.txt'), 'w') as txt:
            for fail in failed:
                txt.write(f'{fail}.txt')

    def extract(self, force=False):
        '''
        Extract the data from the HF and CI files

        Parameters
        ----------
        unify (`bool`)
            if True, mix all results inside one csv file, default is True

        output_keywords (`str`)
            keywords to append at the result file, default is 'AA'
        '''

        # Data extraction
        print(f'Writing data files on {self.output_path}')

        if not os.path.isfile(os.path.join(self.output_path, 'data.csv')) or force:
            searcher = MainLogReader(self.input_path, self.output_path)
            searcher.search()
            searcher.save()

            builder = MolBuilder(self.output_path)

            log_files = [file for file in os.listdir(self.HF_path) if '.log' in file]

            for log in log_files:
                log_path = os.path.join(self.HF_path, log)

                builder.build_xyz(log_path)

            mols = []
            smiles = []
            failed = []

            xyz_route = os.path.join(self.output_path, 'xyz_molecules')
            xyz_files = [file for file in os.listdir(xyz_route) if '.xyz' in file]

            charge_keys = {
            '_c' : 1,
            '_a' : -1,
            '_r' : 0
            }

            for xyz in xyz_files:
                charge = 0
                ID = re.search(r'[/\\]?([A-Z0-9]+[_a-z]*).xyz', xyz).group(1)
                xyz_file = os.path.join(xyz_route, xyz)

                for key in charge_keys.keys():
                    if key in xyz:
                        charge = charge_keys[key]

                smile = builder.get_smiles(xyz_file, charge)
                
                mols.append(ID)
                if not smile:
                    failed.append(ID)
                    smiles.append('')
                    print(f'{ID} failed :c')
                else:
                    smiles.append(smile)
                
                builder.get_image(xyz_file, charge)

            self.__write_failed(failed, self.output_path)
            self.__save_smiles(mols, smiles, self.output_path)