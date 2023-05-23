from utils.prediction_analyzer import prediction_analyzer

database = 'predictions/Net_3Hlayerb_opt4_a27_Full_predictions.csv'
predictions_file = 'predictions/saved_values/Ecorrb4_a-0.27_Full.npy'
alpha = -0.27
b = 'b4'
out = '_Full'

mb = prediction_analyzer(database, predictions_file, alpha, b, out)
mb.make_plots(show=True)