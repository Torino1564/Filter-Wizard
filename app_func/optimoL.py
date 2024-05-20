import numpy as np
from sympy import *
from app_func.approximation_function import *

P, s, e, w = symbols('P, s, e, w')


class OptimoLManager(ApproximationFunctionManager):
    def __init__(self):
        self.C = []
        self.C.append(Poly([1, 0, 0], w))
        self.C.append(Poly([1, 0, 0, 0, 0], w))
        self.C.append(Poly([3, 0, -3, 0, 1, 0, 0], w))
        self.C.append(Poly([20, 0, -40, 0, 28, 0, -8, 0, 1, 0, 0], w))
        self.C.append(Poly([50, 0, -120, 0, 105, 0, -40, 0, 6, 0, 0, 0, 0], w))
        self.C.append(Poly([175, 0, -525, 0, 615, 0, -355, 0, 105, 0, -15, 0, 1, 0, 0], w))
        self.C.append(Poly([490, 0, -1668, 0, 2310, 0, -1624, 0, 615, 0, -120, 0, 10, 0, 0, 0, 0], w))
        self.C.append(Poly([1764, 0, -7056, 0, 11704, 0, -10416, 0, 5376, 0, -1624, 0, 276, 0, -24, 1, 0, 0], w))
        self.C.append(
            Poly([5292, 0, -23520, 0, 44100, 0, -45360, 0, 27860, 0, -10416, 0, 2310, 0, -280, 0, 15, 0, 0, 0, 0], w))
        self.order_counter = self.C.__len__()

    def eval_ol(self, order, value):
        if order > self.order_counter:
            pprint("Maximum order exceeded")
            return np.Inf
        return self.C[order].subs(w, value)

    def find_n(self, template):
        lower_limit = sqrt(((1 / (template.normalized_template.Ga ** 2)) - 1) / (template.xi_val ** 2))

        ol_order = 0
        ol_value = 0
        while ol_value < lower_limit:
            ol_order += 1
            ol_value = self.eval_ol(ol_order, template.normalized_template.w_array[1])
        return ol_order

    def get_approximation_function(self, template):
        self.validate_template(template)
        approximation_function = self.C[template.n].as_expr()
        return sqrt(approximation_function)


manager = OptimoLManager()
