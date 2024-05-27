import template as tmp
import sympy as sp
import numpy as np

class BPTemplate:
    def __init__(self, Gp, Ga, wo, wp1, wa1):
        self.low_pass_template = None
        self.approximation_function = None
        self.Gp = Gp
        self.Ga = Ga
        self.wp1 = wp1
        self.wa1 = wa1
        self.wa2 = wo**2 / wa1
        self.wp2 = wo**2 / wp1
        self.final_function = None
        self.final_function_expr = None

        self.bwp = self.wp2 - wp1
        self.bwa = self.wa2 - wa1
        self.wo = wo

    def set_approximation_function(self, approximation_function):
        self.approximation_function = approximation_function
        self.low_pass_template = tmp.Template(self.Ga, self.Gp, [1, self.bwa / self.bwp], tmp.TemplateType.LP)
        self.low_pass_template.set_approximation_function(self.approximation_function)
        self.low_pass_template.compute()

    def get_low_pass_filter(self):
        return self.low_pass_template

    def generate_final_tf(self):
        t_bp = (self.wo / self.bwp) * ((tmp.s / self.wo) + (self.wo / tmp.s))
        self.final_function_expr = self.low_pass_template.final_function_expr.subs(tmp.s, t_bp)

        num, den = [[float(i) for i in sp.Poly(i, tmp.s).all_coeffs()] for i in
                    sp.simplify(self.final_function_expr).as_numer_denom()]

        self.final_function = tmp.TransferFunction(num, den)
