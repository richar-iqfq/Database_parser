import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error, mean_absolute_percentage_error
import os

class PredictionAnalyzer():
    def __init__(self, database, predictions_file, alpha, b, out=None):
        self.alpha = alpha
        self.database = pd.read_csv(database)
        self.energies = [val for val in self.database['CIe']]
        self.predictions = np.load(predictions_file)
        self.energies_pred = self.predictions[:,0]
        self.b = b
        self.out = out

    def make_plots(self, show=False, save=True):
        '''
        Build and save the resulting plots

        Parameters
        ----------
        show `bool`
            If true display the plots in screen, default is False

        save `bool`
            If true save plots in pdf format inside prediction/plots
        '''
        MAE_e = mean_absolute_error(self.energies, self.energies_pred)
        MAE_b = mean_absolute_error(self.database['B_opt4'], self.database['b_pred'])

        PE_e = mean_absolute_percentage_error(self.energies, self.energies_pred)
        PE_b = mean_absolute_percentage_error(self.database['B_opt4'], self.database['b_pred'])

        #------------------------------- Err ------------------------------
        Err = abs(self.energies - self.energies_pred)
        maxErr = max(Err)
        
        fig, ax = plt.subplot_mosaic('AB;CD')
        fig.suptitle(f'{self.b} alpha: {self.alpha}')
        fig.set_size_inches(20, 13)

        # Energie Err plot
        ax['A'].plot(Err, '.r')
        ax['A'].set_xlabel('Molecules')
        ax['A'].set_ylabel('Absolute Error')
        ax['A'].set_title(f'Error over Energie\nMAE = {MAE_e:.4f}      Max Error: {maxErr:.2f}')

        # Correlation energie plot
        ax['B'].plot(self.energies, self.energies, '-b')
        ax['B'].plot(self.energies, self.energies_pred, '.g')
        
        r2 = np.corrcoef(self.energies, self.energies_pred)

        ax['B'].set_xlabel('Energies')
        ax['B'].set_ylabel('Energies Predicted')
        ax['B'].set_title(f'MAPE: {PE_e*100:.2f}      r$^2$: {r2[0,1]:.3f}')

        # B Err plot
        b_Err = abs(self.database['B_opt4'] - self.database['b_pred'])
        maxBErr = max(b_Err)
        ax['C'].plot(b_Err, '.c')
        ax['C'].set_xlabel('Molecules')
        ax['C'].set_ylabel('Absolute Error')
        ax['C'].set_title(f'Error over b_opt\nMAE = {MAE_b:.4f}      Max Error: {maxBErr:.2f}')

        # Correlation b plot
        ax['D'].plot(self.database['B_opt4'], self.database['B_opt4'], '-b')
        ax['D'].plot(self.database['B_opt4'], self.database['b_pred'], '.y')
        
        r2 = np.corrcoef(self.database['B_opt4'], self.database['b_pred'])

        ax['D'].set_xlabel('b_opt')
        ax['D'].set_ylabel('b Predicted')
        ax['D'].set_title(f'MAPE: {PE_b*100:.2f}      r$^2$: {r2[0,1]:.3f}')

        if show:
            plt.show()

        if save:
            path = os.path.join('predictions', 'plots')
            if not os.path.isdir(path):
                os.makedirs(path)
            
            fig.savefig(os.path.join(path, f'a{self.alpha}_{self.b}{self.out}.pdf'), dpi=450, format='pdf')