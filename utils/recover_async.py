from utils.optimizerV3 import b_optimizer
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
from tqdm import tqdm

#============================== For Ecorr prediction =================================
def update_conv(ans):
    conv_values[ans[0]] = ans[1]
    pbar.update()

def conv1(index):
    ID = IDS[index]
    Ecorr_p, theta_p, mu_p, _ = optimizer[ID].Recover_Ecorr(b_test[ID], 1)
    return index, (Ecorr_p, theta_p, mu_p)

def conv2(index):
    ID = IDS[index]
    Ecorr_p, theta_p, mu_p, _ = optimizer[ID].Recover_Ecorr(b_test[ID], 2)
    return index, (Ecorr_p, theta_p, mu_p)

def conv3(index):
    ID = IDS[index]
    Ecorr_p, theta_p, mu_p, _ = optimizer[ID].Recover_Ecorr(b_test[ID], 3)
    return index, (Ecorr_p, theta_p, mu_p)

def conv4(index):
    ID = IDS[index]
    Ecorr_p, theta_p, mu_p, _ = optimizer[ID].Recover_Ecorr(b_test[ID], 4)
    return index, (Ecorr_p, theta_p, mu_p)

def conv5(index):
    ID = IDS[index]
    Ecorr_p, theta_p, mu_p, _ = optimizer[ID].Recover_Ecorr(b_test[ID], 5)
    return index, (Ecorr_p, theta_p, mu_p)

class calculator():
    def __init__(self, output_path, calc_type=1, cores='All', tolerance=None):
        self.output_path = output_path
        self.calc_type = calc_type
        self.cores = cores
        self.tolerance = tolerance

    def run(self, database, energies, alpha):
        
        global optimizer, IDS, pbar, conv_values, b_test

        root = os.getcwd()

        save_path = os.path.join(root, self.output_path, 'saved_values')
        if not os.path.isdir(save_path):
            os.makedirs(save_path)

        optimizer = {}
        b_test = {}
        final_values = {}

        values_dataset = pd.read_csv(database)
        energies_dataset = pd.read_feather(energies)

        IDS = [id for id in values_dataset['ID']]

        dataset = {}
        for i, ID in enumerate(IDS):
            dataset[ID] = {
                    'NAtoms' : values_dataset['NAtoms'][i],
                    'Ne' : values_dataset['Ne'][i],
                    'CIe' : values_dataset['CIe'][i],
                    'Factor': values_dataset['Factor'][i],
                    'Energies' : [energie for energie in energies_dataset[ID].dropna()],
                    'b_pred' : values_dataset['b_pred'][i]
                }

        for ID in IDS:
            b_test[ID] = dataset[ID]['b_pred']
            orbital_energies = dataset[ID]['Energies']
            Ne = dataset[ID]['Ne']
            CIe = dataset[ID]['CIe']
            factor = dataset[ID]['Factor']
            
            b0 = 1
            mu0 = -0.2
            theta0 = 9.5

            optimizer[ID] = b_optimizer(orbital_energies, Ne, CIe, alpha, b0, mu0, theta0, factor, tolerance=self.tolerance)

        #--------------------------------------------- Run 1 ----------------------------------------------
        if self.calc_type == 1:
            conv_values = np.zeros([len(IDS), 3])

            pool = Pool(processes=self.cores if self.cores != 'All' else None)

            print('.'*20)
            print('Running Expansion 1')
            pbar = tqdm(total=len(IDS), desc='Molecules', colour='cyan')
            
            for i in range(len(IDS)):
                pool.apply_async(conv1, args=(i, ), callback=update_conv)

            pool.close()
            pool.join()
            pbar.close()

            file_name = os.path.join(save_path, f'Ecorrb1_a{alpha}')
            np.save(file_name, conv_values) #Save the results in case of breaking

            final_values['Ecorr_1'] = conv_values[:, 0]
            final_values['Theta_1'] = conv_values[:, 1]
            final_values['Mu_1'] = conv_values[:, 2]

            del(pool)

        #--------------------------------------------- Run 2 ----------------------------------------------
        if self.calc_type == 2:
            conv_values = np.zeros([len(IDS), 3])

            pool = Pool(processes=self.cores if self.cores != 'All' else None)

            print('.'*20)
            print('Running Expansion 2')
            pbar = tqdm(total=len(IDS), desc='Molecules', colour='cyan')
            
            for i in range(len(IDS)):
                pool.apply_async(conv2, args=(i, ), callback=update_conv)

            pool.close()
            pool.join()
            pbar.close()

            file_name = os.path.join(save_path, f'Ecorrb2_a{alpha}')
            np.save(file_name, conv_values) #Save the results in case of breaking

            final_values['Ecorr_2'] = conv_values[:, 0]
            final_values['Theta_2'] = conv_values[:, 1]
            final_values['Mu_2'] = conv_values[:, 2]

            del(pool)

        #--------------------------------------------- Run 3 ----------------------------------------------
        if self.calc_type == 3:
            conv_values = np.zeros([len(IDS), 3])

            pool = Pool(processes=self.cores if self.cores != 'All' else None)

            print('.'*20)
            print('Running Expansion 3')
            pbar = tqdm(total=len(IDS), desc='Molecules', colour='cyan')
            
            for i in range(len(IDS)):
                pool.apply_async(conv3, args=(i, ), callback=update_conv)

            pool.close()
            pool.join()
            pbar.close()

            file_name = os.path.join(save_path, f'Ecorrb3_a{alpha}')
            np.save(file_name, conv_values) #Save the results in case of breaking

            final_values['Ecorr_3'] = conv_values[:, 0]
            final_values['Theta_3'] = conv_values[:, 1]
            final_values['Mu_3'] = conv_values[:, 2]

            del(pool)

        #--------------------------------------------- Run 4 ----------------------------------------------
        if self.calc_type == 4:
            conv_values = np.zeros([len(IDS), 3])

            pool = Pool(processes=self.cores if self.cores != 'All' else None)

            print('.'*20)
            print('Running Expansion 4')
            pbar = tqdm(total=len(IDS), desc='Molecules', colour='cyan')
            
            for i in range(len(IDS)):
                pool.apply_async(conv4, args=(i, ), callback=update_conv)

            pool.close()
            pool.join()
            pbar.close()

            file_name = os.path.join(save_path, f'Ecorrb4_a{alpha}')
            np.save(file_name, conv_values) #Save the results in case of breaking

            final_values['Ecorr_4'] = conv_values[:, 0]
            final_values['Theta_4'] = conv_values[:, 1]
            final_values['Mu_4'] = conv_values[:, 2]

            del(pool)
        
        #--------------------------------------------- Run 5 ----------------------------------------------
        if self.calc_type == 5:
            conv_values = np.zeros([len(IDS), 3])

            pool = Pool(processes=self.cores if self.cores != 'All' else None)

            print('.'*20)
            print('Running Expansion 5')
            pbar = tqdm(total=len(IDS), desc='Molecules', colour='cyan')
            
            for i in range(len(IDS)):
                pool.apply_async(conv5, args=(i, ), callback=update_conv)

            pool.close()
            pool.join()
            pbar.close()

            file_name = os.path.join(save_path, f'Ecorrb5_a{alpha}')
            np.save(file_name, conv_values) #Save the results in case of breaking

            final_values['Ecorr_5'] = conv_values[:, 0]
            final_values['Theta_5'] = conv_values[:, 1]
            final_values['Mu_5'] = conv_values[:, 2]

            del(pool)