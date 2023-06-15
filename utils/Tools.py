import pandas as pd
import os
import re

def raw_data_unify(container, output_keywords):
    '''
    Function that mix all the independent result files together

    Parameters
    ----------
    container (`str`):
        Folder name with the independent raw data extracted, wich contains data.csv and energies.csv
    
    output_keywords (`str`):
        Extra words to name the new database

    Returns
    -------
    file_name (`str`):
        Name of the ouput database.
    '''
    SETS = [set for set in os.listdir('./') if 'set' in set.lower()]
    SETS = sorted(SETS)

    print('Appending data...')

    # For data
    main_name = os.path.join(SETS[0], container, 'data.csv')
    main_df = pd.read_csv(main_name)

    path = os.path.join(container, 'data.csv')

    for set in SETS[1::]:
        print(main_df.shape, end='\r')
        
        df = pd.read_csv(os.path.join(set, path))
        main_df = pd.concat([main_df, df], ignore_index=True)

    main_df = main_df.sort_values(by=['ID'], ignore_index=True)

    file_name = f'./Database_{output_keywords}'
    main_df.to_csv(file_name+'.csv', index=False)

    del(main_df)
    del(df)

    print('Appending energies...')

    # For energies
    main_name = os.path.join(SETS[0], container, 'energies.csv')
    main_df = pd.read_csv(main_name)

    path = os.path.join(container, 'energies.csv')

    for set in SETS[1::]:
        print(main_df.shape, end='\r')
        
        df = pd.read_csv(os.path.join(set, path))
        main_df = pd.concat([main_df, df], axis=1)

    main_df = main_df.loc[:,~main_df.columns.duplicated()].copy()

    main_df.to_feather(f'./Database_Energies_{output_keywords}.feather')

    # For SMILES
    main_name = os.path.join(SETS[0], container, 'smiles.csv')
    main_df = pd.read_csv(main_name)

    path = os.path.join(container, 'smiles.csv')

    for set in SETS[1::]:
        print(main_df.shape, end='\r')
        
        df = pd.read_csv(os.path.join(set, path))
        main_df = pd.concat([main_df, df], ignore_index=True)

    main_df = main_df.sort_values(by=['ID'], ignore_index=True)

    smiles_name = f'./Database_SMILES_{output_keywords}'
    main_df.to_csv(smiles_name+'.csv', index=False)

    return file_name

def raw_data_cleaner(input_file, overwrite=False):
    '''
    Function that cleans the database removing the duplicated molecules and creates a criteria csv file

    Parameters
    ----------
    input_file `str`: 
        Path of the input csv database

    overwrite (`bool`):
        If False, new file with _cleaned will be created to save the new database.
    '''

    path = input_file + '.csv'

    # Load csv
    database = pd.read_csv(path)
    print(f'Initial Database: {database.shape}')

    # Create a new criterium with the HF energies rounded to 1e-5
    database['HFr'] = [round(val, 5) for val in database['HF']]

    # Drop duplicate rows
    cleaned_database = database.drop_duplicates(['NAtoms', 'Ne', 'HFr'])

    # Drop the criterium
    cleaned_database = cleaned_database.drop(columns=['HFr'])

    # Save cleaned database
    if overwrite:
        output_file = input_file + '.csv'
    else:
        output_file = input_file + '_cleaned.csv'
    
    cleaned_database.to_csv(output_file, index=False)

    print(f'Final Database: {cleaned_database.shape}')
    print('Saved')

    # Save the criterium database
    criterium_database = pd.DataFrame().assign(ID=cleaned_database['ID'], NAtoms=cleaned_database['NAtoms'], Ne=cleaned_database['Ne'], HF=cleaned_database['HF'])

    criterium_database.to_csv('criteria.csv', index=False)

def rename_cols(file_name, words):
    df = pd.read_csv(file_name)

    names = df['ID']
    names += words

    df['ID'] = names

    df.to_csv(words + file_name, index=False)

    print('File successfully edited')

def concat_df(files, output_file):
    df = pd.concat(files, ignore_index=True)

    df.to_csv(output_file, index=False)

    print(f'{output_file} writed successfully')

def search_chlorides():
    root = os.getcwd()
    sets = [folder for folder in os.listdir(root) if 'set' in folder.lower()]


    folders = sorted(sets)
    chlorides_list = []

    for folder in folders:
        xyz_files = os.path.join(folder, 'extracted_data', 'xyz_molecules')

        for xyz in os.listdir(xyz_files):
            file = os.path.join(xyz_files, xyz)

            with open(file, 'r') as txt:
                lines = txt.readlines()
                
                for line in lines:
                    if re.search('Cl\t', line):
                        chlorides_list.append(xyz)
                        break

    df = pd.DataFrame({'Chlorides' : chlorides_list})
    df.to_csv('Database_Chlorides_AA.csv', index=False)