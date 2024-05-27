import template as tmp
import sympy as sp
import numpy as np
import scipy as sc
class BRTemplate:
    def __init__(self, Gp, Ga, wo, wp1, wa1):
        self.band_pass_template = None
        self.approximation_function = None
        self.Gp = Gp
        self.Ga = Ga
        self.wp1 = wp1
        self.wa1 = wa1
        self.final_function = None
        self.final_function_expr = None

        self.wp2 = (wo**2)/wp1
        self.wa2 = (wo**2)/wa1
        self.wo = wo
        self.bwp = self.wp2 - self.wp1
        self.bwa = self.wa2 - self.wa1

    def set_approximation_function(self, approximation_function):
        self.approximation_function = approximation_function
        print(self.bwp / self.bwa)
        print(self.bwp)
        print(self.bwa)
        self.low_pass_template = tmp.Template(self.Ga, self.Gp, [1, self.bwp / self.bwa], tmp.TemplateType.LP)
        self.low_pass_template.set_approximation_function(self.approximation_function)
        self.low_pass_template.compute()

    def get_low_pass_filter(self):
        return self.low_pass_template

    def generate_final_tf(self):
        t_br = 1 / ((self.wo / self.bwp) * ((tmp.s / self.wo) + (self.wo / tmp.s)))
        self.final_function_expr = self.low_pass_template.gs_function.subs(tmp.s, t_br)

        num, den = [[float(i) for i in sp.Poly(i, tmp.s).all_coeffs()] for i in
                    sp.simplify(self.final_function_expr).as_numer_denom()]

        zeros, poles, gain = sc.signal.tf2zpk(num, den)

        new_zeros = []
        for zero in zeros:
            if sp.re(zero) < 0:
                new_zeros.append(zero)

        new_poles = []
        for pole in poles:
            if sp.re(pole) < 0:
                new_poles.append(pole)

        if self.approximation_function == tmp.ApproximationFunction.Ch2:
            new_zeros = tmp.remove_duplicates(new_zeros, 3)

        self.final_function = tmp.TransferFunction(sc.signal.ZerosPolesGain(new_zeros, new_poles, np.sqrt(gain)))
