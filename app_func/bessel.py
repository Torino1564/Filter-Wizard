from sympy import *
from app_func.approximation_function import *
P, s, e, w = symbols('P, s, e, w')


class BesselManager(ApproximationFunctionManager):
    def __init__(self):
        self.C = []
        self.C.append(Poly([1], w))
        self.C.append(Poly([1, 0], w))
        self.order_counter = 1

    def add_cheby(self):
        new_poly = self.C[self.order_counter]*w*2 - self.C[self.order_counter - 1]
        self.C.append(new_poly)
        self.order_counter += 1

    def eval_cheby(self, order, value):
        while order > self.order_counter:
            self.add_cheby()
        return self.C[order].subs(w, value)

    def find_n(self, template):
        lower_limit = sqrt(((1 / (template.normalized_template.Ga ** 2)) - 1) / (self.find_xi(template) ** 2))

        ch1_order = 0
        ch1_value = 0
        while ch1_value < lower_limit:
            ch1_order += 1
            ch1_value = self.eval_cheby(ch1_order, template.normalized_template.w_array[1])
        return ch1_order

    def get_approximation_function(self, template):
        self.validate_template(template)
        approximation_function = self.C[template.n].as_expr()
        return approximation_function



manager = BesselManager()