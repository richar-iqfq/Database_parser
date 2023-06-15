from utils.PredictionAnalyzer import PredictionAnalyzer

database = 'predictions/Net_4Hlayerb_opt4_Full_R1_FPredictions.csv'
predictions_file = 'predictions/saved_values/Ecorrb4_a-0.27_R1.npy'
alpha = -0.27
b = 'b4'
out = '_Full_R1'

mb = PredictionAnalyzer(database, predictions_file, alpha, b, out)
mb.make_plots(show=True)