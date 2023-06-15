from utils.RecoverAsync import Calculator

calc_type = 4
alpha = -0.27

database = f'Net_4Hlayerb_opt{calc_type}_a27_predictions.csv'
energies = 'Database_Energies_BB.feather'

tolerance = {
    'Ne_calc' : 1E-4,
    'Energy_S' : 1E-7,
    'Energy_Savrg' : 1E-7,
    'Recover' : 1E-7
}

mb = Calculator('predictions', calc_type, cores=5, tolerance=tolerance)
mb.run(database, energies, alpha)