from utils.Extractor import Extractor
from utils.TerminationChecker import TerminationChecker

if __name__=='__main__':
    '''
    Extract and unify data in one csv file
    '''

    input_path = 'log_files'
    output_path = 'extracted_data'

    t_checker = TerminationChecker(input_path)
    t_checker.check(move=True)

    mb = Extractor(input_path, output_path)
    mb.extract(force=True)