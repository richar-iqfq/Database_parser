from utils.searcher import main_log_reader
from utils.tools import raw_data_cleaner, raw_data_unify
import os

class extractor():
    '''
    Extract the data from the HF and CI files into csv contained inside output_path

    Parameters
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

    def extract(self, unify=True, output_keywords='AA'):
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
            if not os.path.isfile(os.path.join(self.root, folder, self.output_path, 'data.csv')):
                searcher = main_log_reader(folder, self.output_path)
                searcher.search()
                searcher.Save()

        if unify:
            print('Writing and cleaning final dataset')
            file_name = raw_data_unify(self.output_path, output_keywords)
            raw_data_cleaner(file_name, overwrite=True)