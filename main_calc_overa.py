from utils.CalculatorAsyncV6 import MainCalculator
# from utils.convergence import convergence_writter
from utils.StatisticAnalyser import StatisticAnalyser
import os

'''
Script that unifies the recovering of the data, the calculation of each b and
the corroboration of the molecules convergence for given alpha values
as well as the statistical analysis for the percent values.
'''

if __name__=='__main__':
    cores = 14
    calc_types=(1, 2, 3, 4, 5)
    alpha_list = (-0.27, 0, 0)
    
    percent = (0.9, 0.95, 0.96, 0.97, 1, 1.03, 1.04, 1.05, 1.1)

    tolerance = {
        'Ne_calc' : 1E-4,
        'Energy_S' : 1E-7,
        'Energy_Savrg' : 1E-7,
        'Recover' : 1E-7
    }

    database = 'Database_BB.csv'
    energies = 'Database_Energies_BB.feather'

    total_calcs = len(alpha_list)
    count = 1

    for alpha in alpha_list:
        print(f'Running: {count}/{total_calcs}\n')
        output_path = f'a{alpha}_results'

        # calculator = MainCalculator(output_path, calc_types=calc_types, cores=cores, tolerance=tolerance)
        # calculator.run_optimization(database, energies, alpha)

        out_file = os.path.join('results', 'optimization', output_path, f'results_a{alpha}.csv')

        # c_checker = convergence_writter(out_file)
        # c_checker.run_write()

        st = StatisticAnalyser(out_file, energies, output_path, alpha=alpha, percent=percent, specific=4)
        st.plot_dispersion(save=True, show=False)
        st.plot_correlation_matrix(save=True, show=False, diag=False)
        #st.plot_Err(cores=cores, save=True, show=False, tolerance=tolerance)

        count += 1
