import app_func.cheby as ch1
from numpy import *

class Cheby2Manager(ch1.ChebyManager):

    def find_xi(self, template):
        return 1/sqrt(template.normalized_template.Ga ** 2 / ( 1 + template.normalized_template.Ga ** 2))

    def find_n(self, template):
        return ch1.manager.find_n(template)

    def get_approximation_function(self, template):
        self.validate_template(template)
        approximation_function = 1 / ch1.manager.C[template.n].as_expr().subs(ch1.w, template.normalized_template.w_array[1] / ch1.w)
        return approximation_function


manager = Cheby2Manager()
