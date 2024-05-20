from enum import Enum

import scipy.signal
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from sympy import init_printing
plt.rcParams.update({
    "text.usetex": False
})
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
        temp_function = simplify(self.fw_function).subs(w, s / I)
        self.gs_function = temp_function

    def generate_final_tf(self):
        num, den = [[float(i) for i in Poly(i, s).all_coeffs()] for i in self.gs_function.as_numer_denom()]
        zeros, poles, gain = signal.tf2zpk(num, den)

        new_poles = []
        for pole in poles:
            if re(pole) <= 0:
                new_poles.append(pole)

        self.final_function = signal.TransferFunction(signal.ZerosPolesGain(zeros, new_poles, sqrt(abs(gain))))

    approximation_function = ApproximationFunction.Butter
    basic_template = BasicTemplate
    normalized_template = BasicTemplate

    def print_filter(self, normalized: bool = True):
        assert (self.fw_function is not None)
        template = None
        if normalized:
            template = self.normalized_template
        else:
            template = self.basic_template

        num = self.final_function.num
        den = self.final_function.den
        final_function_expr = (Poly(num, s)) / Poly(den, s)
        tf_abs = abs(final_function_expr.subs(s, I * w))
        dp.display(tf_abs, "Transfer Function Absolute Value")
        numerical_function = lambdify(w, tf_abs, 'numpy')
        x_values = linspace(0, 2.5, 400)
        y_values = numerical_function(x_values)

        # Plot the function using Matplotlib
        plt.plot(x_values, y_values, label=str("filter"))

        # Pass Region
        square_x = 0
        square_y = template.Gp
        square_height = 1 - template.Gp
        square_width = 1

        rectangle_pass = Rectangle((square_x, square_y), square_width, square_height, fill=True, color='gray',
                                   alpha=0.3)
        plt.gca().add_patch(rectangle_pass)

        # Attenuation Region
        square_x = template.w_array[1]
        square_y = 0
        square_height = template.Ga
        square_width = 2.5 - square_x

        rectangle_att = Rectangle((square_x, square_y), square_width, square_height, fill=True, color='gray',
                                  alpha=0.3)
        plt.gca().add_patch(rectangle_att)

        plt.xlabel('w')
        plt.ylabel('f(w)')
        plt.legend()
        plt.grid(True)

        plt.show()

    def show_transformations(self):
        # Initial function
        dp.display(
            "The goal is to have a function that is symmetric (either even or odd) and is restricted between -1 and 1 "
            "for x values between -1 and 1")
        dp.display("This is the approximation function Y(w):")

        numerical_function = lambdify(w, self.approximation_function_expr, 'numpy')
        x_values = linspace(-3, 5, 600)
        y_values = numerical_function(x_values)

        plt.plot(x_values, y_values)

        square_x = -1
        square_y = -1
        square_width = 2
        square_height = 2
        rectangle = Rectangle((square_x, square_y), square_width, square_height, fill=True, color='grey', alpha=0.3)

        plt.gca().add_patch(rectangle)
        plt.xlabel(r'$\omega$')
        plt.ylabel(r'$Y(\omega)$')
        plt.xlim(-2, 2)
        plt.ylim(-3, 3)
        plt.grid(True)

        plt.show()

        dp.display("Then, the vertical scale of the function is reduced to accommodate for the height of the pass band."
                   "This is achieved by multiplying the function by the xi factor")
        numerical_function = lambdify(w, self.xi_val * self.approximation_function_expr, 'numpy')
        y_values = numerical_function(x_values)

        plt.plot(x_values, y_values, label=str(r'$Y(\omega)\xi$'))

        square_x = -1
        square_y = -self.xi_val
        square_width = 2
        square_height = 2 * self.xi_val
        rectangle = Rectangle((square_x, square_y), square_width, square_height, fill=True, color='grey', alpha=0.3)

        plt.gca().add_patch(rectangle)
        plt.xlabel(r'$\omega$')
        plt.ylabel(r'$Y(\omega)\xi$')
        plt.xlim(-2, 2)
        plt.ylim(-3, 3)
        plt.grid(True)

        plt.show()

        print("To ensure the function is even and its codomain is the positive real numbers, the entire expression is "
              "squared")
        numerical_function = lambdify(w, (self.xi_val * self.approximation_function_expr) ** 2, 'numpy')
        y_values = numerical_function(x_values)

        plt.plot(x_values, y_values)

        square_x = -1
        square_y = 0
        square_width = 2
        square_height = self.xi_val ** 2
        rectangle = Rectangle((square_x, square_y), square_width, square_height, fill=True, color='grey', alpha=0.3)

        plt.gca().add_patch(rectangle)
        plt.xlabel(r'$\omega$')
        plt.ylabel(r'$(Y(\omega)\xi)^2$')
        plt.xlim(-2, 2)
        plt.ylim(0, 3)
        plt.grid(True)

        plt.show()

        print("Since the normalized low pass filter has a DC gain of 0 (1 in linear scale), a 1 is added to the "
              "expression so that the starting point of the function is 1")

        numerical_function = lambdify(w, (self.xi_val * self.approximation_function_expr) ** 2 + 1, 'numpy')
        y_values = numerical_function(x_values)

        plt.plot(x_values, y_values)

        square_x = -1
        square_y = 1
        square_width = 2
        square_height = self.xi_val ** 2
        rectangle = Rectangle((square_x, square_y), square_width, square_height, fill=True, color='grey', alpha=0.3)

        plt.gca().add_patch(rectangle)
        plt.xlabel(r'$\omega$')
        plt.ylabel(r'$(Y(\omega)\xi)^2 + 1$')
        plt.xlim(-2, 2)
        plt.ylim(0, 3)
        plt.grid(True)

        plt.show()

        print("Lastly, the entire expression is inverted so that the it behaves like a low pass amplitude transfer "
              "function")

        numerical_function = lambdify(w, 1 / ((self.xi_val * self.approximation_function_expr) ** 2 + 1), 'numpy')
        y_values = numerical_function(x_values)

        plt.plot(x_values, y_values)

        square_x = 0
        square_y = self.normalized_template.Gp ** 2
        square_width = 1
        square_height = 1 - self.normalized_template.Gp ** 2
        rectangle = Rectangle((square_x, square_y), square_width, square_height, fill=True, color='grey', alpha=0.3)

        plt.gca().add_patch(rectangle)
        plt.xlabel(r'$\omega$')
        plt.ylabel(r'$\frac{1}{(Y(\omega)\xi)^2 + 1}$')
        plt.xlim(0, 4)
        plt.ylim(0, 1)
        plt.grid(True)

        plt.show()


def normalize_w_array(w_array, template_type):
    if template_type == TemplateType.LP:
        w_array[1] /= w_array[0]
        w_array[0] = 1
    if template_type == TemplateType.HP:
        w_array[0] /= w_array[1]
        w_array[1] = 1
    # TODO: Complete for band pass and band reject
    return w_array
