import scipy.signal
from numpy import *
import matplotlib.pyplot as plt
import scipy.signal as sg
import IPython.display as dp

def log_base(value, base):
    result = log(value) / log(base)
    return result


def print_pzm(transfer_function: sg.TransferFunction):

    poles = transfer_function.poles
    zeros = transfer_function.zeros

    # Plot poles and zeros
    plt.scatter(poles.real, poles.imag, marker='x', color='red', label='Poles')
    plt.scatter(zeros.real, zeros.imag, marker='o', color='blue', label='Zeros')

    # Add labels and legend
    plt.xlabel('Real')
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
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

    pole_count = 1
    # print poles and zeros
    dp.display("Poles:")
    for pole in poles:
        dp.display("{}) Wo:{}, Phase:{}, Q:{}" .format(pole_count, abs(pole), angle(pole), absolute(1/(2*sin(absolute(angle(pole) - pi/2))))))
        pole_count += 1

    dp.display("Zeros:")
    for zero in zeros:
        dp.display("{}) Wo:{}, Phase:{}".format(pole_count, abs(zero), angle(zero)))

