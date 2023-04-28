from utils.calculator_asyncV6 import main_calculator

'''
Script that unifies the recovering of the data, the calculation of each b and
the corroboration of the convergence for the molecules 
'''

alpha = (-0.27)
calc_types = (0, 0, 0, 4, 0)

if __name__=='__main__':
    cores = 50
    tolerance = {
        'Ne_calc' : 1E-7,
        'Energy_S' : 1E-9,
        'Energy_Savrg' : 1E-9,
        'Recover' : 1E-9
    }

    database = 'Database_AA.csv'
    energies = 'Database_Energies_AA.feather'

    output_path = f'a{alpha}_final-results'

    calculator = main_calculator(output_path, calc_types=calc_types, cores=cores, tolerance=tolerance)
    calculator.run_optimization(database, energies, alpha)