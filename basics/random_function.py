#!/usr/bin/env python3.7

from boolean_function import BooleanFunction

from arybo.lib import MBA

import random
import itertools

import sys


class RandomFunction(BooleanFunction):

    # {0,1}^in_size -> {0,1}^out_size
    def __init__(self, in_size, out_size, seed):
        if in_size > 0:
            mba = MBA(1)

            self.v = [mba.var("x0")]
            for i in range(1, in_size):
                self.v.append(mba.var("x" + str(i)))

            random.seed(seed)

            self.x = []
            for i in range(out_size):
                self.x.append(self._random_poly())


    # create random polynomial in ANF from self.v variables
    def _random_poly(self):
        mba = MBA(1)
        polynom = mba.var("")

        polynom_expr = self.pytanque_expr_const(random.random() < (1 / 2.0))

        indices = list(range(len(self.v)))

        for i in indices:
            possibilities = [self._create_expr(list(poss)) for poss in itertools.combinations(indices, i + 1)]

            for possibility in possibilities:
                pick = random.random() < (1 / 2.0)
                if pick:
                    polynom_expr = polynom_expr + possibility

        polynom = (polynom | 0x1) & polynom_expr

        return polynom




# mandatory command line arguments:
#  sys.argv[1]: number of input bits of the function
#  sys.argv[2]: number of output bits of the function
#  sys.argv[3]: seed
#  sys.argv[4]: file of values to be evaluated by the created function (one value per line)
# optional command line arguments
#  sys.argv[5]: file for output values

# TODO: better manipulation with command line arguments
# TODO: error control?
# TODO: optional logging

if __name__ == "__main__":
    if len(sys.argv) > 4:
        in_size = int(sys.argv[1])
        out_size = int(sys.argv[2])
        seed = int(sys.argv[3])
        input_values = sys.argv[4]
        output_file = sys.argv[5] if len(sys.argv) > 5 else sys.stdout

        function = RandomFunction(in_size, out_size, seed)

        with open(input_values, "r") as f:
            for value in f:
                print(function.eval_f(int(value)), file=output_file)
        print(function.x)
        print(function)
    else:
        print("Usage: ./random_function.py input_size output_size seed input_values_file", file=sys.stderr)
        exit(1)
