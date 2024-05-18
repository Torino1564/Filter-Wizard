from enum import Enum

import scipy.signal
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
        self.final_function = None
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
        return self.approximation_function_expr

    def compute(self):
        self.find_xi()
        self.find_n()
        self.get_correct_approx_function()
        self.generate_fw_function()
        self.generate_gs_function()
        self.generate_final_tf()

    def set_approximation_function(self, approximation_function):
        self.approximation_function = approximation_function
        self.updated = false

    def generate_fw_function(self):
        temp_function = 1 / (1 + (self.xi_val * self.approximation_function_expr) ** 2)
        self.fw_function = temp_function

    def generate_gs_function(self):
        temp_function = self.fw_function.subs(w, s / I)
        self.gs_function = temp_function

    def generate_final_tf(self):
        num, den = [[float(i) for i in Poly(i, s).all_coeffs()] for i in self.gs_function.as_numer_denom()]
        zeros, poles, gain = signal.tf2zpk(num, den)

        new_poles = []
        for pole in poles:
            if re(pole) <= 0:
                new_poles.append(pole)

        self.final_function = signal.TransferFunction(signal.ZerosPolesGain(zeros, new_poles, gain))

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
