from pypdf import PdfMerger

import os

if __name__=='__main__':
    calc_types=(1, 2, 3, 4, 5)
    alpha_list = (-0.06, -0.08, -0.1, -0.12, -0.14, -0.16, -0.18,
                  -0.2, -0.22, -0.24, -0.26, -0.28, -0.3, -0.34,
                  -0.38, -0.42, -0.46)
    
    total_calcs = len(alpha_list)
    count = 1

    out_file = os.path.join('results', 'recovering', 'Statistic_Plots', 'b_dispersion_mean.pdf')

    path = os.path.join('results', 'recovering', 'Statistic_Plots')
    
    pdf_list = []

    merger = PdfMerger()

    for alpha in alpha_list:
        print(f'Merging: {count}/{total_calcs}\n')

        pdf = os.path.join(path, f'a{alpha}_results', f'b_dispersion_a{alpha}.pdf')

        merger.append(pdf)
        
        count += 1

    merger.write(out_file)
    merger.close()