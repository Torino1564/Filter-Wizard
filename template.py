from enum import Enum
from scipy import signal
import matplotlib.pyplot as plt
from sympy import init_printing
init_printing()
from app_func.full import *
from utils import *

import IPython.display as ipd

class TemplateType(Enum):
    LP = 1
    HP = 2
    BP = 3
    BR = 4


class BasicTemplate:
    Ga = 0
    Gp = 0
    w_array = [0]
    template_type = TemplateType.LP


class Template:
    def __init__(self, ga, gp, w_array, template_type):
        self.fw_function = None
        self.gs_function = None
        self.approximation_function_expr = None
        self.xi_val = 0
        self.order = 0
        self.updated = false
        self.basic_template.Ga = ga
        self.basic_template.Gp = gp
        self.basic_template.w_array = w_array
        self.basic_template.template_type = template_type

        self.normalized_template.Ga = ga
        self.normalized_template.Gp = gp
        self.normalized_template.w_array = normalize_w_array(w_array, template_type)
        self.normalized_template.template_type = template_type

        self.find_xi()

    def find_xi(self):
        manager = managers[self.approximation_function.value]
        self.xi_val = manager.find_xi(self)
        return self.xi_val

    def find_n(self):
        manager = managers[self.approximation_function.value]
        self.order = manager.find_n(self)
        return self.order

    def get_correct_approx_function(self):
        manager = managers[self.approximation_function.value]
        self.approximation_function_expr = manager.get_approximation_function(self)
        ipd.display(self.approximation_function_expr)
        return self.approximation_function_expr

    def compute(self):
        self.find_xi()
        self.find_n()
        self.get_correct_approx_function()
        self.generate_fw_function()
        self.cull_positive_poles()
        self.print_pzm()

    def set_approximation_function(self, approximation_function):
        self.approximation_function = approximation_function
        self.updated = false

    def generate_fw_function(self):
        temp_function = 1 / (1 + (self.xi_val * self.approximation_function_expr) ** 2)
        self.fw_function = temp_function.subs(w, s / I)

    def cull_positive_poles(self):
        numerator, denominator = self.fw_function.as_numer_denom()

        array_denominator = denominator.as_poly().all_coeffs()

        if isinstance(numerator, Number):
            array_numerator = [float(numerator)]
        else:
            array_numerator = numerator.as_poly().all_coeffs()

        builtin_format_array_denominator = list(map(float, array_denominator))
        builtin_format_array_numerator = list(map(float, array_numerator))
        zeros, poles, gain = signal.tf2zpk(builtin_format_array_numerator, builtin_format_array_denominator)

        new_denominator = 1
        new_numerator = 1

        for pole in poles:
            if re(pole) > 0:
                pass
            else:
                new_denominator *= (1 - (s / pole))

        for zero in zeros:
            new_numerator *= (1 - (s / zero))

        self.gs_function = (gain * new_numerator) / new_denominator

    def print_pzm(self):
        assert (self.gs_function is not None)
        numerator, denominator = self.gs_function.as_numer_denom()

        array_denominator = denominator.as_poly().all_coeffs()

        if isinstance(numerator, Number):
            array_numerator = [float(numerator)]
        else:
            array_numerator = numerator.as_poly().all_coeffs()

        builtin_format_array_denominator = list(map(float, array_denominator))
        test_poly = Poly(builtin_format_array_denominator, w)
        #poles = solve(test_poly, w)
        builtin_format_array_numerator = list(map(float, array_numerator))
        zeros, poles, gain = signal.tf2zpk(builtin_format_array_numerator, builtin_format_array_denominator)

        # Plot poles and zeros
        plt.scatter(poles.real, poles.imag, marker='x', color='red', label='Poles')
        plt.scatter(zeros.real, zeros.imag, marker='o', color='blue', label='Zeros')

        # Add labels and legend
        plt.xlabel('Real')
        plt.ylabel('Imaginary')
        plt.title('Poles and Zeros Map')
        plt.legend()

        # Add unit circle (optional)
        unit_circle = plt.Circle((0, 0), 1, fill=False, linestyle='dotted', color='gray')
        plt.gca().add_artist(unit_circle)

        # Set aspect ratio to equal
        plt.gca().set_aspect('equal', adjustable='box')

        # Show plot
        plt.grid(True)
        plt.show()

    approximation_function = ApproximationFunction.Butter
    basic_template = BasicTemplate
    normalized_template = BasicTemplate


def normalize_w_array(w_array, template_type):
    if template_type == TemplateType.LP:
        w_array[1] /= w_array[0]
        w_array[0] = 1
    if template_type == TemplateType.HP:
        w_array[0] /= w_array[1]
        w_array[1] = 1
    # TODO: Complete for band pass and band reject
    return w_array
