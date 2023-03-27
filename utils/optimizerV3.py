from scipy.optimize import minimize
import numpy as np

'''
Optimizer script, calculate the optimization for the parameters given to the class, resturns b, mu and theta
'''

class b_optimizer():
    '''
    Attributes
    ----------
    orbital_energies (`list` or `array`):
        Vector with all the orbital energies

    Ne (`int`):
        Number of electrons

    Ecorr (`float`):
        Ecorr energie

    alpha (`float`):
        Alpha values for the optimization

    b0 (`float`):
        Initial value for b_opt

    mu0 (`float`):
        Initial value for mu

    theta0 (`float`):
        Initial value for theta

    factor (`int`):
        1 or 2, depends on calculation type
    '''

    def __init__(self, orbital_energies, Ne, Ecorr, alpha, b0, mu0, theta0, factor, tolerance=None):
        # Fixed parameters
        self.orbital_energies = orbital_energies
        self.Ne = Ne # Electron number
        self.Ecorr = Ecorr # Correlation energie
        self.alpha = alpha/factor # Alpha
        self.factor = factor # Factor

        # Initial values
        self.b0 = b0
        self.theta0 = theta0
        self.mu0 = mu0

        # Optimized Parameters
        self.Ne_calc = 0 # Electron number calculated (inside class)
        self.b_calc = 0 # b calculated (inside class)
        self.theta_calc = 0 # theta calculated (inside class)
        self.mu_calc = 0 # mu calculated (inside class)

        if tolerance == None:
            # This tolerance values are too high, take so many time (about 30s per molecule)
            self.tolerance = {
                'Ne_calc' : 1E-11,
                'Energy_S' : 1E-15,
                'Energy_Savrg' : 1E-15,
                'Recover' : 1E-11
            }
        else:
            self.tolerance = tolerance

    def Ne_calculation(self, theta=None, mu=None):
        if theta==None and mu==None:
            theta = self.theta_calc
            mu = self.mu_calc
            
        sum = 0
        for energy in self.orbital_energies:
            ni = self.Occupation(energy, theta, mu)
            sum += ni

        ne = self.factor*sum

        return ne

    def Occupation(self, energy, theta, mu):
        # occupation = ni
        # con ni = {1 + exp[theta(energy - mu)]}^-1

        # Calculation protection
        if (theta * (energy - mu)) < -700:
            ni = 1
        elif (theta * (energy - mu)) > 700:
            ni = 0
        else:
            ni = 1 / (1 + np.exp(theta*(energy - mu)))

        return ni

    def Entropy_VN(self, theta, mu):
        # S(theta) = Sum(ni*Ln(ni) + (1 - ni)*Ln(1 - ni))

        S = 0
        for energy in self.orbital_energies:
            ni = self.Occupation(energy, theta, mu)
            # Calculation protection
            if ni == 0:
                a = 0
            else:
                a = ni * np.log(ni)

            if ni == 1:
                b = 0
            else:
                b = (1 - ni) * np.log(1 - ni)

            Si = a + b
            S += Si
        
        return S

    def Energy_S(self, theta, mu):
        # E(S) = S/theta

        S = self.Entropy_VN(theta, mu)
        Energy_S = ((self.factor*S)/theta)

        return Energy_S

    def Entropy_avrg(self, expansion, theta, mu, b):
        # ci = ni*(1 - ni)
        # <S(theta)> = alpha + (theta/2Ne)*Sum(ci*[betha*mu - energy])

        Sum_a = 0
        Sum_b = 0
        Sum_c = 0
        Sum_d = 0
        Sum_e = 0

        for energy in self.orbital_energies:
            ni = self.Occupation(energy, theta, mu)
            ci = ni*(1 - ni)

            if expansion >= 1:
                sum_ai = ci*(b*mu - energy)
                Sum_a += sum_ai
            
            if expansion >= 2:
                sum_bi = ci*(1 - 2*ni)*(b*mu - energy)**2
                Sum_b += sum_bi

            if expansion >= 3:
                sum_ci = ci*(1 - 6*ci)*(b*mu - energy)**3
                Sum_c += sum_ci

            if expansion >= 4:
                sum_di = ci*(1 - 2*ni)*(1 - 12*ci)*(b*mu - energy)**4
                Sum_d += sum_di

            if expansion >= 5:
                sum_ei = ci*(1 - 30*ci*(1 - 4*ci))*(b*mu - energy)**5
                Sum_e += sum_ei

        expansion_terms = 0

        if expansion >= 1:
            term_1 = (theta*Sum_a)/(self.Ne)

            expansion_terms += term_1

        if expansion >= 2:
            term_2 = ((theta**2)*Sum_b)/(2*self.Ne)

            expansion_terms += term_2

        if expansion >= 3:
            term_3 = ((theta**3)*Sum_c)/(6*self.Ne)

            expansion_terms += term_3

        if expansion >= 4:
            term_4 = ((theta**4)*Sum_d)/(24*self.Ne)

            expansion_terms += term_4

        if expansion >= 5:
            term_5 = ((theta**5)*Sum_e)/(120*self.Ne)

            expansion_terms += term_5

        S_avrg = self.alpha + expansion_terms

        return S_avrg

    def Energy_Savrg(self, expansion, theta=None, mu=None, b_test=None):
        # E(<S>) = (Ne*<S>)/theta
        if theta==None and mu==None:
            theta = self.theta_calc
            mu = self.mu_calc
        if b_test==None:
            b_test = self.b_calc

        S = self.Entropy_avrg(expansion, theta, mu, b_test)
        
        Energy_Savrg = ((self.factor*S*self.Ne)/theta)

        return Energy_Savrg

    def __Mu_optimizer(self, theta):
        # Optimizer for mu
        def errorN(mu_test):
            Ne_calc = self.Ne_calculation(theta, mu_test)
            error = abs(self.Ne - Ne_calc)

            # print(f'Ne_calc: {error}')
            return error[0]

        mu_local = self.mu0
        result_mu = minimize(errorN, mu_local, method='Nelder-Mead', tol=self.tolerance['Ne_calc'])

        self.mu_calc = result_mu.x[0]
        self.Ne_calc = self.Ne_calculation(theta, self.mu_calc)
        
        return self.mu_calc
    
    def Mu_theta_optimizer(self):
        # Optimizer for thata and mu
        def errorS(theta_test):
            mu_i = self.__Mu_optimizer(theta_test)
            Energy_Si = self.Energy_S(theta_test, mu_i)

            error = abs(Energy_Si - self.Ecorr)
            # print(f'E_S: {error}')
            return error

        theta_local = self.theta0
        result_theta = minimize(errorS, theta_local, method='Nelder-Mead', tol=self.tolerance['Energy_S'])

        self.theta_calc = result_theta.x[0]
        self.Ne_calc = self.Ne_calculation()

        return self.theta_calc, self.mu_calc

    def Run_boptimization(self, expansion):
        # Optimizer for b
        if abs(self.Ne - self.Ne_calculation()) >= 1E-3:
            _, _ = self.Mu_theta_optimizer()
            
        def error_Savrg(b):
            Energy_Savrg = self.Energy_Savrg(expansion, b_test=b)
            
            error = abs(Energy_Savrg - self.Ecorr)
            # print(f'Type {expansion}: {error}')
            return error
        
        options = {'maxiter': 300}
        def run():
            b_local = self.b0
            result = minimize(error_Savrg, b_local, method='Nelder-Mead', tol=self.tolerance['Energy_Savrg'], options=options)

            return result
        
        result_b = run()
        self.b_calc = result_b.x[0]

        return self.b_calc

    def Recover_Ecorr(self, bopt, type):
        # Error function
        def errorB(theta_test):
            mu_i = self.__Mu_optimizer(theta_test)

            Energy_avrg = self.Energy_Savrg(type, theta=theta_test, mu=mu_i, b_test=bopt)
            Energy_S = self.Energy_S(theta_test, mu_i)

            error = abs(Energy_avrg - Energy_S)

            # print(f'R_MT {type}: {error}')
            return error

        theta_local = self.theta0
        result_theta = minimize(errorB, theta_local, method='Nelder-Mead', tol=self.tolerance['Recover']) #Try to increase the convergence to Nelder

        mu_predict = self.mu_calc
        theta_predict = result_theta.x[0]

        Ecorr_p = self.Energy_Savrg(type, theta=theta_predict, mu=mu_predict, b_test=bopt)
        
        err = abs(Ecorr_p - self.Ecorr)
        # print(f'Final Err: {err}')

        return Ecorr_p, theta_predict, mu_predict, err#, bopt #remove b

    def get_values(self, type):

        theta_opt = self.theta_calc
        mu_opt = self.mu_calc
        b_opt = self.b_calc

        S_avrg = self.Energy_Savrg(expansion=type)
        Err = abs(S_avrg - self.Ecorr)

        return theta_opt, mu_opt, b_opt, Err

    def Convergence(self, theta, mu, b, type=1):
        Energy = self.Energy_Savrg(type, theta=theta, mu=mu, b_test=b)

        convergence = abs(Energy - self.Ecorr)
        if convergence <= 1E-3:
            Flag = True
        else:
            Flag = False
        
        return Flag