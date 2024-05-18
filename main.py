# dependencies

from app_func.full import *
from template import *

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    example_template = Template(0.01, 0.80, [1, 2], TemplateType.LP)

    example_template.set_approximation_function(ApproximationFunction.Butter)
    example_template.compute()
    pprint(example_template.get_correct_approx_function())

    # example_template.set_approximation_function(ApproximationFunction.Ch1)
    # xi_2 = example_template.find_xi()
    # order_2 = example_template.find_n()
    # pprint(example_template.get_correct_approx_function())

    print("This is the last time")
