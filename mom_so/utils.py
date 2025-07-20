""" It is a Python script that contains the main function and a class called User. """

import os
import numpy as np


# Clear the screen
def clear_screen():
    """ Clears the console screen. """
    os.system('cls' if os.name == 'nt' else 'clear')


def complex_formatter(x):
    """ Format complex numbers in scientific notation."""
    if x == 0:
        return " 0 "
    else:
        return f"{x.real:.3e}{x.imag:+.3e}j"

# Set the print options to use the custom formatter
def set_numpy_print_options():
    """ Set the numpy print options to use the custom formatter."""
    # Define a custom formatter for the numpy print options
    my_formatter = {
        'float_kind': lambda x: f"{x:.2e}",
        'complex_kind': complex_formatter
    }
    np.set_printoptions(formatter=my_formatter)
