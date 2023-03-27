from utils.dataloader import dataloader
from utils.optimizerV3 import b_optimizer
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
from tqdm import tqdm

'''
Main Script, extract all the information from the log files and calculate the b, theta and mu for each molecule.
'''
# ========================= For optimization ================================
def update(ans):
    optim_values[ans[0]] = ans[1]
    pbar.update()

def launcher1(index):
    ID = IDS[index]
    optimizer[ID].Run_boptimization(1)
    theta_opt, mu_opt, b_opt, Err = optimizer[ID].get_values(type=1)

    return index, (theta_opt, mu_opt, b_opt, Err)

def launcher2(index):
    ID = IDS[index]
    optimizer[ID].Run_boptimization(2)
    theta_opt, mu_opt, b_opt, Err = optimizer[ID].get_values(type=2)

    return index, (theta_opt, mu_opt, b_opt, Err)

def launcher3(index):
    ID = IDS[index]
    optimizer[ID].Run_boptimization(3)
    theta_opt, mu_opt, b_opt, Err = optimizer[ID].get_values(type=3)

    return index, (theta_opt, mu_opt, b_opt, Err)

def launcher4(index):
    ID = IDS[index]
    optimizer[ID].Run_boptimization(4)
    theta_opt, mu_opt, b_opt, Err = optimizer[ID].get_values(type=4)

    return index, (theta_opt, mu_opt, b_opt, Err)

def launcher5(index):
    ID = IDS[index]
    optimizer[ID].Run_boptimization(5)
    theta_opt, mu_opt, b_opt, Err = optimizer[ID].get_values(type=5)

    return index, (theta_opt, mu_opt, b_opt, Err)

#============================== For Ecorr prediction =================================
def update_conv(ans):
    conv_values[ans[0]] = ans[1]
    pbar.update()

def conv1(index):
    ID = IDS[index]
    Ecorr_p, theta_p, mu_p, _ = optimizer[ID].Recover_Ecorr(b1_test[ID], 1)
    return index, (Ecorr_p, theta_p, mu_p)

def conv2(index):
    ID = IDS[index]
    Ecorr_p, theta_p, mu_p, _ = optimizer[ID].Recover_Ecorr(b2_test[ID], 2)
    return index, (Ecorr_p, theta_p, mu_p)

def conv3(index):
    ID = IDS[index]
    Ecorr_p, theta_p, mu_p, _ = optimizer[ID].Recover_Ecorr(b3_test[ID], 3)
    return index, (Ecorr_p, theta_p, mu_p)

def conv4(index):
    ID = IDS[index]
    Ecorr_p, theta_p, mu_p, _ = optimizer[ID].Recover_Ecorr(b4_test[ID], 4)
    return index, (Ecorr_p, theta_p, mu_p)

def conv5(index):
    ID = IDS[index]
    Ecorr_p, theta_p, mu_p, _ = optimizer[ID].Recover_Ecorr(b5_test[ID], 5)
    return index, (Ecorr_p, theta_p, mu_p)

class main_calculator():
    '''
    Calculates both, the optimized b and the Ecorr recovering

    Attributes
    ----------
    output_path (`str`):
        Output path to storage results, for `run_optimization`: root/results/optimization/output_path, 
        for `run_recovering`: root/results/recovering/output_path

    calc_types (`tuple`) of (`int`):
        Expansion type to compute calculation, valid values: 1, 2, 3, 4, 5, default->(1, 2, 3, 4, 5)
    
    cores ('All') or (`int`):
        Number of cores/threads to use in calculation.

    tolerance (`dict`):
        Dictionary with the tolerance values. Keys->'Ne_calc','Energy_S','Energy_Savrg','Recover'

    Methods
    -------
    run_obtimization(database, energies, alpha)
        Computes the b optimization and generate a csv file with the results and folder with the .npy
        files.
        
        Parameters
        ----------
        database (`str`):
            database.csv path
        
        energies (`str`):
            energies.feather path

        alpha (`float`):
            alpha value

    run_recovering(database, energies, percent, alpha):
        Computes the Ecorr, theta and mu recovering and generate a folder with the .npy files

        Parameters
        ----------
        database (`str`):
            database.csv path
        
        energies (`str`):
            energies.feather path

        percent (`float`):
            percent value of b for recovering, 1 as 100%
            
        alpha (`float`):
            alpha value
    '''
    def __init__(self, output_path, calc_types=(1, 2, 3, 4, 5), cores='All', tolerance=None, options=None):
        self.output_path = output_path
        self.calc_types = calc_types
        self.cores = cores
        self.tolerance = tolerance
        self.options = options

    def run_optimization(self, database, energies, alpha):
        '''
        Computes the b optimization and generate a csv file with the results and folder with the .npy
        files.
        
        Parameters
        ----------
        database (`str`):
            database.csv path
        
        energies (`str`):
            energies.feather path

        alpha (`float`):
            alpha value
        '''
        global optimizer, IDS, pbar, optim_values

        os.system('clear')
        root = os.getcwd()
        save_path = os.path.join(root, 'results', 'optimization', self.output_path)
        if not os.path.isdir(save_path):
            os.makedirs(save_path)

        save_path2 = os.path.join(root, 'results', 'optimization', self.output_path, 'saved_values')
        if not os.path.isdir(save_path2):
            os.makedirs(save_path2)
        
        print(f'Processing Alpha = {alpha}...')
        print(f'\nCores = {self.cores}')
        
        optimized_values = {
            'Theta' : [],
            'Mu' :  [],
            'B_opt1' : [],
            'B_opt2' : [],
            'B_opt3' : [],
            'B_opt4' : [],
            'B_opt5' : [],
            'Err1' : [],
            'Err2' : [],
            'Err3' : [],
            'Err4' : [],
            'Err5' : [],
        }

        optimizer = {}

        loader = dataloader(database, energies)
        dataset = loader.load(raw=True)
        IDS = loader.get_IDS()

        for ID in IDS:
            Ne = dataset[ID]['Ne']
            Ecorr = dataset[ID]['CIe']
            factor = dataset[ID]['Factor']
            orbital_energies = dataset[ID]['Energies']

            if self.options != None:
                b0 = self.options['b0']
            else:
                b0 = 1
            
            mu0 = -0.2
            theta0 = 9.5
            
            optimizer[ID] = b_optimizer(orbital_energies, Ne, Ecorr, alpha, b0, mu0, theta0, factor, tolerance=self.tolerance)

        if 1 in self.calc_types:
            # Calc for b1
            optim_values = np.zeros([len(IDS), 4])

            pool = Pool(processes=self.cores if self.cores != 'All' else None)

            print('.............')
            print('Optimizing b1')
            pbar = tqdm(total=len(IDS), desc='Molecules', colour='cyan')
            for i in range(len(IDS)):
                pool.apply_async(launcher1, args=(i,), callback=update)
            
            pool.close()
            pool.join()
            pbar.close()

            for i in range(len(optim_values)):
                optimized_values['Theta'].append(optim_values[i][0])
                optimized_values['Mu'].append(optim_values[i][1])
                optimized_values['B_opt1'].append(optim_values[i][2])
                optimized_values['Err1'].append(optim_values[i][3])

            file_name = os.path.join(save_path2, f'optim1_a{alpha}')
            np.save(file_name, optim_values) #Save the results in case of breaking


        if 2 in self.calc_types:
            # Calc for b2
            optim_values = np.zeros([len(IDS), 4])

            pool = Pool(processes=self.cores if self.cores != 'All' else None)

            print('.............')
            print(f'Optimizing b2')        
            pbar = tqdm(total=len(IDS), desc='Molecules', colour='cyan')
            for i in range(len(IDS)):
                pool.apply_async(launcher2, args=(i,), callback=update)
            
            pool.close()
            pool.join()
            pbar.close()

            for i in range(len(optim_values)):
                optimized_values['B_opt2'].append(optim_values[i][2])
                optimized_values['Err2'].append(optim_values[i][3])

            file_name = os.path.join(save_path2, f'optim2_a{alpha}')
            np.save(file_name, optim_values) #Save the results in case of breaking

        if 3 in self.calc_types:
            # Calc for b3
            optim_values = np.zeros([len(IDS), 4])

            pool = Pool(processes=self.cores if self.cores != 'All' else None)

            print('.............')
            print(f'Optimizing b3')
            pbar = tqdm(total=len(IDS), desc='Molecules', colour='cyan')
            for i in range(len(IDS)):
                pool.apply_async(launcher3, args=(i,), callback=update)
            
            pool.close()
            pool.join()
            pbar.close()

            for i in range(len(optim_values)):
                optimized_values['B_opt3'].append(optim_values[i][2])
                optimized_values['Err3'].append(optim_values[i][3])

            file_name = os.path.join(save_path2, f'optim3_a{alpha}')
            np.save(file_name, optim_values) #Save the results in case of breaking

        if 4 in self.calc_types:
            # Calc for b4
            optim_values = np.zeros([len(IDS), 4])

            pool = Pool(processes=self.cores if self.cores != 'All' else None)

            print('.............')
            print(f'Optimizing b4')
            pbar = tqdm(total=len(IDS), desc='Molecules', colour='cyan')
            for i in range(len(IDS)):
                pool.apply_async(launcher4, args=(i,), callback=update)
            
            pool.close()
            pool.join()
            pbar.close()

            for i in range(len(optim_values)):
                optimized_values['B_opt4'].append(optim_values[i][2])
                optimized_values['Err4'].append(optim_values[i][3])

            file_name = os.path.join(save_path2, f'optim4_a{alpha}')
            np.save(file_name, optim_values) #Save the results in case of breaking

        if 5 in self.calc_types:
            # Calc for b4
            optim_values = np.zeros([len(IDS), 4])

            pool = Pool(processes=self.cores if self.cores != 'All' else None)

            print('.............')
            print(f'Optimizing b5')
            pbar = tqdm(total=len(IDS), desc='Molecules', colour='cyan')
            for i in range(len(IDS)):
                pool.apply_async(launcher5, args=(i,), callback=update)
            
            pool.close()
            pool.join()
            pbar.close()

            for i in range(len(optim_values)):
                optimized_values['B_opt5'].append(optim_values[i][2])
                optimized_values['Err5'].append(optim_values[i][3])

            file_name = os.path.join(save_path2, f'optim5_a{alpha}')
            np.save(file_name, optim_values) #Save the results in case of breaking

        print('Optimization Completed :D')
        print('Saving data...')

        if os.path.isfile(os.path.join(save_path, f'results_a{alpha}.csv')):
            database_df = pd.read_csv(os.path.join(save_path, f'results_a{alpha}.csv'))
        else:
            database_df = pd.read_csv(database)
        
        # Save values in one dataframe
        for key in optimized_values.keys():
            if len(optimized_values[key]) != 0:
                database_df[key] = optimized_values[key]
        
        if len(self.calc_types) != 0:
            database_df.to_csv(os.path.join(save_path, f'results_a{alpha}.csv'), index=False)

        if not os.path.isfile(os.path.join(save_path, f'results_a{alpha}.csv')):
            print(f'results_a{alpha}.csv not saved, something went wrong or any calc specified!')
            input('Press Enter to continue...')

        print('Completed')

    def run_recovering(self, database, energies, percent, alpha):
        '''
        Computes the Ecorr, theta and mu recovering and generate a folder with the .npy files

        Parameters
        ----------
        database (`str`):
            database.csv path
        
        energies (`str`):
            energies.feather path

        percent (`float`):
            percent value of b for recovering, 1 as 100%
            
        alpha (`float`):
            alpha value
        '''
        global optimizer, IDS, pbar, conv_values, b1_test, b2_test, b3_test, b4_test, b5_test
        
        root = os.getcwd()
        save_path = os.path.join(root, 'results', 'recovering', 'saved_values', self.output_path)
        if not os.path.isdir(save_path):
            os.makedirs(save_path)

        optimizer = {}
        b1_test = {}
        b2_test = {}
        b3_test = {}
        b4_test = {}
        b5_test = {}
        final_values = {}

        loader = dataloader(database, energies)
        dataset = loader.load(raw=False)
        IDS = loader.get_IDS()

        for ID in IDS:
            orbital_energies = dataset[ID]['Energies']
            Ne = dataset[ID]['Ne']
            CIe = dataset[ID]['CIe']
            factor = dataset[ID]['Factor']
            b1_test[ID] = dataset[ID]['B_opt1']*percent
            b2_test[ID] = dataset[ID]['B_opt2']*percent
            b3_test[ID] = dataset[ID]['B_opt3']*percent
            b4_test[ID] = dataset[ID]['B_opt4']*percent
            b5_test[ID] = dataset[ID]['B_opt5']*percent

            b0 = 1
            mu0 = -0.2
            theta0 = 9.5

            optimizer[ID] = b_optimizer(orbital_energies, Ne, CIe, alpha, b0, mu0, theta0, factor, tolerance=self.tolerance)

        # ====================== Run 1 ===================================
        if 1 in self.calc_types:
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

            file_name = os.path.join(save_path, f'Ecorrb1_a{alpha}_perc{percent}')
            np.save(file_name, conv_values) #Save the results in case of breaking

            final_values['Ecorr_1'] = conv_values[:, 0]
            final_values['Theta_1'] = conv_values[:, 1]
            final_values['Mu_1'] = conv_values[:, 2]

            del(pool)

        # ====================== Run 2 ===================================
        if 2 in self.calc_types:
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

            file_name = os.path.join(save_path, f'Ecorrb2_a{alpha}_perc{percent}')
            np.save(file_name, conv_values) #Save the results in case of breaking

            final_values['Ecorr_2'] = conv_values[:, 0]
            final_values['Theta_2'] = conv_values[:, 1]
            final_values['Mu_2'] = conv_values[:, 2]
            
            del(pool)

        # ====================== Run 3 ===================================
        if 3 in self.calc_types:
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

            file_name = os.path.join(save_path, f'Ecorrb3_a{alpha}_perc{percent}')
            np.save(file_name, conv_values) #Save the results in case of breaking

            final_values['Ecorr_3'] = conv_values[:, 0]
            final_values['Theta_3'] = conv_values[:, 1]
            final_values['Mu_3'] = conv_values[:, 2]
            
            del(pool)

        # ====================== Run 4 ===================================
        if 4 in self.calc_types:
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

            file_name = os.path.join(save_path, f'Ecorrb4_a{alpha}_perc{percent}')
            np.save(file_name, conv_values) #Save the results in case of breaking

            final_values['Ecorr_4'] = conv_values[:, 0]
            final_values['Theta_4'] = conv_values[:, 1]
            final_values['Mu_4'] = conv_values[:, 2]
            
            del(pool)

        # ====================== Run 5 ===================================
        if 5 in self.calc_types:
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

            file_name = os.path.join(save_path, f'Ecorrb5_a{alpha}_perc{percent}')
            np.save(file_name, conv_values) #Save the results in case of breaking

            final_values['Ecorr_5'] = conv_values[:, 0]
            final_values['Theta_5'] = conv_values[:, 1]
            final_values['Mu_5'] = conv_values[:, 2]
            
            del(pool)

        print('Finished')

        return final_values