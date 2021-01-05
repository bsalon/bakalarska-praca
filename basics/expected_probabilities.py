#!/usr/bin/env python3.7

from arybo.lib import MBA

import sys

import itertools


# returns True/False wrapped in pytanque.Expr
def pytanque_expr_const(boolean):
    mba = MBA(1)
    x = mba.var("x")
    x = x | 1 if boolean else x & 0
    return x[0]


class BooleanFunction:

    def __init__(self, monomials, degree):
        mba = MBA(1)

        self.vars = [mba.var("x0")]    # TODO: naming may not be the best
        for i in range(1, degree):
            self.vars.append(mba.var("x" + str(i)))

        self._create_function(monomials)


    # create monomial from chosen indices
    def _create_expr(self, indices):
        if len(indices) == 0:
            return pytanque_expr_const(False)

        expression = pytanque_expr_const(True)
        for index in indices:
            expression *= self.vars[index][0]

        return expression


    # create function in ANF form from monomials
    def _create_function(self, monomials):
        mba = MBA(1)
        polynom = mba.var("")

        polynom_expr = pytanque_expr_const([] in monomials)

        for monomial in monomials:
            monomial_expr = self._create_expr(monomial)

            polynom_expr = polynom_expr + monomial_expr

        self.function = (polynom | 0x1) & polynom_expr


    # evaluate the function on value input
    def eval_f(self, value):
        # we need to set each bit to appropriate value
        dictionary = {self.vars[var] : (value & (0x1 << var)) >> var for var in range(len(self.vars))}
        return self.function.eval(dictionary)


def create_key_string(polynom):
    if len(polynom) == 0:
        return ""

    key_string = ""

    for monomial in polynom[:-1]:
        if len(monomial) == 0:
            key_string += "1 + "
            continue
        for var in monomial:
            key_string += "x" + str(var)
        key_string += " + "

    for var in polynom[-1]:
        key_string += "x" + str(var)

    return key_string


def compute_possibilities(degree):
    if degree > 4:
        degree = 4

    monomials = [[]]

    variables = list(range(degree))
    for i in range(degree):
        monomials.extend([list(monomial) for monomial in itertools.combinations(variables, i + 1)])

    polynoms = []
    for i in range(len(monomials)):
        polynoms.extend([list(polynom) for polynom in itertools.combinations(monomials, i + 1)])

    possibilities = {}
    for polynom in polynoms:
        bool_function = BooleanFunction(polynom, degree)

        one = 0
        for val in range(2 ** degree):
            one += bool_function.eval_f(val)

        possibilities[create_key_string(polynom)] = one / (2.0 ** degree)

    return possibilities


# mandatory command line argument
#  sys.argv[1]: number of variables = maximal degree of monomial
if __name__ == "__main__":
    degree = int(sys.argv[1])
    expectations = compute_possibilities(degree)

    for polynom, probability in expectations.items():
        print("{:>125}   {}".format(polynom, probability))
