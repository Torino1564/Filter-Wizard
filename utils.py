import numpy as np
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

    plt.figure(figsize=(8, 8))

    # Plot poles and zeros
    plt.scatter(poles.real, poles.imag, marker='x', color='red', label='Poles')
    plt.scatter(zeros.real, zeros.imag, marker='o', color='blue', label='Zeros')

    # Add labels and legend
    plt.xlabel('Real')
    plt.ylabel('Imaginary')
    plt.title('Poles and Zeros Map')
    plt.legend()
    # plt.gca().set_aspect('equal', adjustable='box')

    # Add unit circle (optional)
    unit_circle = plt.Circle((0, 0), 1, fill=False, linestyle='dotted', color='gray')
    plt.gca().add_artist(unit_circle)

    # Set aspect ratio to equal
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.gca().autoscale_view()
    # Show plot
    plt.grid(True)
    plt.show()
    pole_count = 1
    zero_count = 1
    # print poles and zeros
    dp.display("Poles:")
    for pole in poles:
        dp.display("{}) Wo:{}, Phase:{}, Q:{}".format(pole_count, abs(pole), angle(pole),
                                                      absolute(1 / (2 * sin(absolute(angle(pole) - pi / 2))))))
        pole_count += 1

    dp.display("Zeros:")
    for zero in zeros:
        dp.display("{}) Wo:{}, Phase:{}".format(zero_count, abs(zero), angle(zero)))
        zero_count += 1


def print_bode(transfer_function, wi=0, wf=0, points=0, wx_array=None, wy_array=None):
    if wx_array is None:
        wx_array = []
    w_array = None
    if wi != 0 and wf != 0:
        if points != 0:
            w_array = np.logspace(wi, wf, points)
        else:
            w_array = np.logspace(wi, wf)

    numerator = transfer_function.num
    denominator = transfer_function.den

    system = sg.TransferFunction(numerator, denominator)

    if w_array is None:
        freq, amplitude, phase = sg.bode(system)
    else:
        freq, amplitude, phase = sg.bode(system, w_array)

    plt.figure()

    if wx_array is not None:
        for w in wx_array:
            plt.axvline(x=w, color='red', linestyle='--')

    if wy_array is not None:
        for fy in wy_array:
            plt.axhline(y=fy, color='blue', linestyle='--')

    plt.semilogx(freq, amplitude)
    plt.title("Filter Bode Amplitude Diagram")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Amplitude [dB]")
    plt.grid(which='both', axis='both')
    plt.show()


def round_complex(num, decimals):
    return complex(round(num.real, decimals), round(num.imag, decimals))


def remove_duplicates(complex_array, decimals=2):
    seen = set()
    unique_list = []
    for num in complex_array:
        rounded_num = round_complex(num, decimals)
        if rounded_num not in seen:
            seen.add(rounded_num)
            unique_list.append(num)  # Append the original number
    return np.array(unique_list)
