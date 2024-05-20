import template as tmp


class HPTemplate:
    def __init__(self, Gp, Ga, wp, wa, ):
        self.low_pass_template = None
        self.approximation_function = None
        self.Gp = Gp
        self.Ga = Ga
        self.wp = wp
        self.wa = wa
        self.final_function = None

    def set_approximation_function(self, approximation_function):
        self.approximation_function = approximation_function
        self.low_pass_template = tmp.Template(self.Ga, self.Gp, [1, self.wp / self.wa], tmp.TemplateType.LP)
        self.low_pass_template.compute()

    def get_low_pass_filter(self):
        return self.low_pass_template

    def generate_final_tf(self):
        self.final_function = self.low_pass_template.final_function_expr.subs(tmp.s, self.wp / tmp.s)
