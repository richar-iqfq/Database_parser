from utils.Convergence import Analyzer
from utils.CalculatorAsyncV6 import MainCalculator

if __name__=='__main__':
    database = 'Database_AA.csv'
    energies = 'Database_Energies_AA.feather'

    tolerance = {
        'Ne_calc' : 1E-4,
        'Energy_S' : 1E-7,
        'Energy_Savrg' : 1E-7,
        'Recover' : 1E-7
    }

    re_path = r'a-[0-9.]+_results'

    analyzer = Analyzer(re_path)
    values = analyzer.analyze()

    for new_try in values:
        alpha, type, b_mean = new_try[0], new_try[1], new_try[2]

        output_path = f'a{alpha}_results'
        options = {
            'b0' : b_mean
        }

        mb = MainCalculator(output_path, calc_types=(type, 0, 0, 0, 0), cores=14, tolerance=tolerance, options=options)
        mb.run_optimization(database, energies, alpha)