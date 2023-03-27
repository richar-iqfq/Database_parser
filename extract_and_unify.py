from utils.extractor import extractor

if __name__=='__main__':
    '''
    Extract and unify data in one csv file
    '''
    
    output_path = 'extracted_data'
    output_keywords = 'AA'

    mb = extractor(output_path)
    mb.extract(unify=True, output_keywords=output_keywords)