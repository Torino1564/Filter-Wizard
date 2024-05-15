# dependencies

from IPython.display import display
from sympy import *
from numpy import *
from sympy.abc import s, w, n, xi
from enum import Enum

from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

# global variables
n = symbols('n')
xi = symbols('xi')
F = symbols('F')
G = symbols('G')
P, s, e, w = symbols('P, s, e, w')


class ChebyManager:
    def __init__(self):
        self.C = []
        self.C.append(Poly([0], w))
        self.C.append(Poly([1, 0], w))
        self.order_counter = 1

    def add_cheby(self):
        self.C.append(self.C[self.order_counter - 1] * w * 2 - self.C[self.order_counter - 2])
        self.order_counter += 1

    def eval_cheby(self, order, value):
        while order > self.order_counter:
            self.add_cheby()
        return self.C[order].subs(w, value)


cheby_manager = ChebyManager()


class TemplateType(Enum):
    LP = 1
    HP = 2
    BP = 3
    BR = 4


class ApproximationFunction(Enum):
    Butter = 1
    Ch1 = 2
    Ch2 = 3


class BasicTemplate:
    Ga = 0
    Gp = 0
    w_array = [0]
    template_type = TemplateType.LP


def log_base(value, base):
    try:
        result = log(value) / log(base)
        return result
    except (ValueError, TypeError) as err:
        print("Error:", err)
        return None


class Template:
    def __init__(self, ga, gp, w_array, template_type):
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
        self.xi_val = ((1 / self.normalized_template.Gp ** 2) - 1)

    def find_n(self, approximation_function):
        assert (self.xi_val != 0)
        if approximation_function == ApproximationFunction.Butter:
            temp = (1 - self.normalized_template.Ga ** 2) / ((self.normalized_template.Ga ** 2) * self.xi_val ** 2)
            n_precise = log_base(temp, self.normalized_template.w_array[1])/2
            final_n = int(n_precise + 1)
            return final_n
        if approximation_function == ApproximationFunction.Ch1:
            lower_limit = sqrt(((1 / (self.normalized_template.Ga ** 2)) - 1) / xi)

            ch1_order = 0
            ch1_value = 0
            while ch1_value < lower_limit:
                ch1_value = cheby_manager.eval_cheby(ch1_order, self.normalized_template.w_array[1])

            return ch1_order

    basic_template = BasicTemplate
    normalized_template = BasicTemplate
    xi_val = 0


def normalize_w_array(w_array, template_type):
    if template_type == TemplateType.LP:
        w_array[1] /= w_array[0]
        w_array[0] = 1
    if template_type == TemplateType.HP:
        w_array[0] /= w_array[1]
        w_array[1] = 1
    # TODO: Complete for band pass and band reject
    return w_array


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    example_template = Template(0.01, 0.707, [1000, 4000], TemplateType.LP)
    example_template.find_n(ApproximationFunction.Butter)
