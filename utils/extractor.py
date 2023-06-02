from utils.searcher import main_log_reader
from utils.tools import raw_data_cleaner, raw_data_unify
from utils.mol_builder import Mol_builder
import os
import re
import pandas as pd

class extractor():
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
    def __init__(self, output_path):
        self.root = os.getcwd()
        sets = [folder for folder in os.listdir(self.root) if 'set' in folder.lower()]
        
        self.folders = sorted(sets)
        self.output_path = output_path

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

    def extract(self, unify=True, output_keywords='AA', charge=0):
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
        for folder in self.folders:
            folder_path = os.path.join(self.root, folder)
            if not os.path.isfile(os.path.join(folder_path, self.output_path, 'data.csv')):
                searcher = main_log_reader(folder, self.output_path)
                searcher.search()
                searcher.Save()

            HF_path = [os.path.join(folder_path, HF_path) for HF_path in os.listdir(folder_path) if 'HF' in HF_path]
            builder = Mol_builder(os.path.join(folder_path, self.output_path), charge=charge)

            log_files = [file for file in os.listdir(HF_path[0]) if '.log' in file]

            for log in log_files:
                log_path = os.path.join(HF_path[0], log)

                builder.build_xyz(log_path)

            mols = []
            smiles = []
            failed = []

            xyz_route = os.path.join(folder_path, self.output_path, 'xyz_molecules')
            xyz_files = [file for file in os.listdir(xyz_route) if '.xyz' in file]

            for xyz in xyz_files:
                ID = re.search(r'[/\\]?([A-Z0-9]+)[_a-z]*.xyz', xyz).group(1)
                xyz_file = os.path.join(xyz_route, xyz)
                smile = builder.get_smiles(xyz_file)
                
                mols.append(ID)
                if not smile:
                    failed.append(ID)
                    smiles.append('')
                    print(f'{ID} failed :c')
                else:
                    smiles.append(smile)
                
                builder.get_image(xyz_file)

            self.__write_failed(failed, os.path.join(folder_path, self.output_path))
            self.__save_smiles(mols, smiles, os.path.join(folder_path, self.output_path))

        if unify:
            print('Writing and cleaning final dataset')
            file_name = raw_data_unify(self.output_path, output_keywords)
            raw_data_cleaner(file_name, overwrite=True)

    