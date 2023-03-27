from utils.calculator_asyncV6 import main_calculator
from utils.convergence import convergence_checker
from utils.termination_checker import check_Normal_T
import os

'''
Script that unifies the recovering of the data, the calculation of each b and
the corroboration of the convergence for the molecules 
'''

if __name__=='__main__':
    cores = 20
    calc_types= (1, 2, 3, 4, 5)
    tolerance = {
        'Ne_calc' : 1E-4,
        'Energy_S' : 1E-7,
        'Energy_Savrg' : 1E-7,
        'Recover' : 1E-7
    }

    database = 'Database_II.csv'
    energies = 'Database_Energies_II.feather'

    t_checker = check_Normal_T()
    t_checker.check(move=True)

    alpha = -0.5
    output_path = f'a{-0.5}_results'

    calculator = main_calculator(output_path, calc_types=calc_types, cores=cores, tolerance=tolerance)
    calculator.run_optimization(database, energies, alpha)

    out_file = os.path.join('results', 'optimization', output_path, 'results_a{alpha}.csv')
    c_checker = convergence_checker(out_file)
    c_checker.run_check()