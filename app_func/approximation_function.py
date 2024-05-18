from sympy import *
from abc import ABC, abstractmethod
from enum import Enum


class ApproximationFunction(Enum):
    Null = 0
    Butter = 1
    Ch1 = 2
    Cauer = 3


class ApproximationFunctionManager:
    @abstractmethod
    def get_approximation_function(self, template):
        pass

    @staticmethod
    def find_xi(template):
        return float(sqrt((1 / template.normalized_template.Gp ** 2) - 1))

    @abstractmethod
    def find_n(self, template):
        pass

    def validate_template(self, template):
        if template.updated == false:
            template.xi_val = self.find_xi(template)
            template.n = self.find_n(template)
            template.updated = true
        if template.xi_val == 0:
            template.xi_val = self.find_xi(template)
        if template.n == 0:
            template.n = self.find_n(template)
        template.updated = true
