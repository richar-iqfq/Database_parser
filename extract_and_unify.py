from utils.Extractor import Extractor
from utils.TerminationChecker import CheckNormalT

if __name__=='__main__':
    '''
    Extract and unify data in one csv file
    '''

    t_checker = CheckNormalT()
    # t_checker.check(move=True)
    
    output_path = 'extracted_data'
    output_keywords = 'AA'

    mb = Extractor(output_path)
    mb.extract(unify=True, output_keywords=output_keywords, charge=0)