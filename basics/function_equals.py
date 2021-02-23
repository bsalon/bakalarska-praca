#!/usr/bin/env python3.7

from create_all_functions import AllBooleanFunctions
from random_function import RandomFunction

import sys

import matplotlib.pyplot as plt
import numpy as np


def parse_values(input_values):
    result = []
    with open(input_values, "r") as f:
        for value in f:
            result.append(int(value))

    return result


def is_equal_to(all_f, seed, min_c, values, constraints, plot=False):
    function = RandomFunction(all_f.vars, 1, seed)

    function_results = []
    for value in values:
        function_results.append(function.eval_f(value))

    length = len(values)
    it_length = range(length)

    functions_equals = []

    for i_function in all_f.functions:
        current_results = []
        for value in values:
            current_results.append(i_function.eval_f(value))
        
        # compute eq coefficient
        equals = len([1 for i in it_length if current_results[i] == function_results[i]]) / length
        if equals > min_c:
            functions_equals.append((i_function.__str__(), equals))

    functions_equals.sort(key=lambda tup: (tup[1], -len(tup[0])), reverse=True)

    # save the results into plot_file
    if plot_file != None:
        if len(functions_equals) > 40:
            functions_equals = functions_equals[:40]
        functions_equals.sort(key=lambda tup: (tup[1], -len(tup[0])), reverse=False)
        plot_me(all_f.vars, values, function_results, functions_equals, constraints, plot_file)
    # print the results
    else:
        for (function, equals) in functions_equals:
            print("(" + function + ") %.2f" % equals)


def plot_me(in_vars, in_vals, out_vals, results, constraints, out_file):
    figure = plt.figure(figsize=(5, 10))

    # functions strings
    names = tuple([name for name, value in results])
    
    y_pos = np.arange(len(names))

    # functions values
    correct = [value for name, value in results]
    plt.barh(y_pos, correct)

    plt.yticks(y_pos, names)
    
    title_string = "function {0,1}^" + str(in_vars) + " -> {0,1}\n"
    
    # variables
    if in_vars < 11:
        for i in range(in_vars):
            title_string += str(in_vars - i - 1)
        title_string += "  \n"

    # results
    bin_string = "0" + str(in_vars) + "b"
    for i in range(len(in_vals)):
        title_string += format(in_vals[i], bin_string) + " " + str(out_vals[i]) + "\n"

    # max degree, max terms
    title_string += constraints

    plt.title(title_string)
    plt.subplots_adjust(top=0.7)

    figure.savefig(out_file, bbox_inches="tight")


if __name__ == "__main__":
    if len(sys.argv) > 7:
        in_size = int(sys.argv[1])
        max_deg = int(sys.argv[2])
        max_terms = int(sys.argv[3])
        max_vars = int(sys.argv[4])
        seed = int(sys.argv[5])
        min_c = float(sys.argv[6])
        input_values = parse_values(sys.argv[7])
        plot_file = sys.argv[8] if len(sys.argv) > 8 else None

        all_functions = AllBooleanFunctions(in_size, max_deg, max_terms, max_vars)

        constraints = "max_degree=" + str(max_deg) + " max_terms=" + str(max_terms)
        is_equal_to(all_functions, seed, min_c, input_values, constraints, plot_file)

    else:
        print("Usage: ./function_equals.py input_size max_degree max_terms max_vars seed min_coefficient input_values_file [plot_file]", file=sys.stderr)
        exit(1)
