from sympy import *
from utils import *
from app_func.approximation_function import *
P, s, e, w = symbols('P, s, e, w')


class ButterManager(ApproximationFunctionManager):
    def find_n(self, template):
        temp = (1 - template.normalized_template.Ga ** 2) / ((template.normalized_template.Ga ** 2) * template.xi_val ** 2)
        n_precise = log_base(temp, template.normalized_template.w_array[1]) / 2
        final_n = int(n_precise + 1)
        return final_n

    def get_approximation_function(self, template):
        self.validate_template(template)
        return w ** template.order

manager = ButterManager()
