from sympy import *
from utils import *
from app_func.approximation_function import *

P, s, e, w = symbols('P, s, e, w')
from scipy.signal import *
from scipy.optimize import root_scalar


class CauerManager(ApproximationFunctionManager):
    def find_n(self, template):
        # template parameters
        xi_a = self.find_xi_a(template)
        xi = template.xi_val

        order = ellipord(template.normalized_template.w_array[0],
                         template.normalized_template.w_array[0],
                         -20 * log_base(template.normalized_template.Gp, 10),
                         -20 * log_base(template.normalized_template.Ga, 10),
                         True)


        return final_n

    def get_approximation_function(self, template):
        self.validate_template(template)
        return w ** template.order

    @staticmethod
    def find_xi_a(template):
        return float(sqrt((1 / template.normalized_template.Ga ** 2) - 1))


manager = CauerManager()
