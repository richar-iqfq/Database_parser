from utils.calculator_asyncV6 import main_calculator
from utils.convergence import convergence_checker
import os

'''
Script that unifies the recovering of the data, the calculation of each b and
the corroboration of the convergence for the molecules 
'''

alpha_list = (-0.3, -0.14, -0.34, -0.24)
calc_types = {
    -0.30 : (1, 0, 0, 0, 0),
    -0.14 : (0, 2, 0, 0, 0),
    -0.34 : (0, 0, 3, 0, 0),
    -0.24 : (0, 0, 0, 4, 0)
}

if __name__=='__main__':
    cores = 50
    tolerance = {
        'Ne_calc' : 1E-7,
        'Energy_S' : 1E-9,
        'Energy_Savrg' : 1E-9,
        'Recover' : 1E-9
    }

    database = 'Database_II.csv'
    energies = 'Database_Energies_II.feather'


    for alpha in alpha_list:
        
        output_path = f'a{alpha}_final-results'

        calculator = main_calculator(output_path, calc_types=calc_types[alpha], cores=cores, tolerance=tolerance)
        calculator.run_optimization(database, energies, alpha)

        out_file = os.path.join('results', 'optimization', output_path, f'results_a{alpha}.csv')
        c_checker = convergence_checker(out_file)
        c_checker.run_check()