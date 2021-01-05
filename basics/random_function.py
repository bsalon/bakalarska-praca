#!/usr/bin/env python3.7

from arybo.lib import MBA

import random
import itertools

import sys
    

# returns True/False wrapped in pytanque.Expr
def pytanque_expr_const(boolean):
    mba = MBA(1)
    x = mba.var("x")
    x = x | 1 if boolean else x & 0
    return x[0]


class RandomFunction:

    # {0,1}^size -> {0,1}^size
    def __init__(self, size, seed):
        mba = MBA(1)

        self.v = [mba.var("v0")]    # TODO: naming may not be the best
        for i in range(1, size):
            self.v.append(mba.var("v" + str(i)))
        
        random.seed(seed)

        self.x = []
        for i in range(size):    
            self.x.append(self._random_polynom())


    # create monomial from chosen indices
    def _create_expr(self, indices):
        if len(indices) == 0:
            return pytanque_expr_const(False)

        expression = pytanque_expr_const(True)
        for index in indices:
            expression *= self.v[index][0]

        return expression


    # create random polynom in ANF form from self.v variables
    def _random_polynom(self):
        mba = MBA(1)
        polynom = mba.var("")

        polynom_expr = pytanque_expr_const(random.random() < (1 / 2.0))
        
        indices = list(range(len(self.v)))

        for i in indices:
            possibilities = [self._create_expr(list(poss)) for poss in itertools.combinations(indices, i + 1)]

            for possibility in possibilities:
                pick = random.random() < (1 / 2.0)
                if pick:
                    polynom_expr = polynom_expr + possibility

        polynom = (polynom | 0x1) & polynom_expr

        return polynom
        


    # evaluate the function on value input
    def eval_f(self, value):
        # we need to set each bit to appropriate value
        dictionary = {self.v[var] : (value & (0x1 << var)) >> var for var in range(len(self.v))}
        result = 0
        
        # evaluate each bit of the function
        for bit in range(len(self.x)):
            thisone = self.x[bit].eval(dictionary)
            result += thisone << bit

        return result


# mandatory command line arguments:
#  sys.argv[1]: number of input-output bits of the function
#  sys.argv[2]: seed
#  sys.argv[3]: file of values to be evaluated by the created function (one value per line)
# optional command line arguments
#  sys.argv[4]: file for output values

# TODO: better manipulation with command line arguments
# TODO: error control?
# TODO: optional logging

if __name__ == "__main__":
    if len(sys.argv) > 3:
        size = int(sys.argv[1])
        seed = int(sys.argv[2])
        input_values = sys.argv[3]
        output_file = sys.argv[4] if len(sys.argv) > 4 else sys.stdout

        function = RandomFunction(size, seed)

        with open(input_values, "r") as f:
            for value in f:
                print(function.eval_f(int(value)), file=output_file)
    else:
        print("Usage: ./random_function.py (input/output)_size seed input_values_file", file=sys.stderr)
        exit(1)
