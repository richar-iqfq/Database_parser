from utils.extractor import extractor
from utils.termination_checker import check_Normal_T

if __name__=='__main__':
    '''
    Extract and unify data in one csv file
    '''

    t_checker = check_Normal_T()
    # t_checker.check(move=True)
    
    output_path = 'extracted_data'
    output_keywords = 'AA'

    mb = extractor(output_path)
    mb.extract(unify=True, output_keywords=output_keywords)