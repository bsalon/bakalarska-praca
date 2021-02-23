#!/usr/bin/env python3.7 

from boolean_function import BooleanFunction

import itertools

import sys


class AllBooleanFunctions:
    
    # TODO: max vars
    def __init__(self, in_size, max_deg=-1, max_terms=-1, max_vars=-1):
        max_deg = in_size if max_deg == -1 else max_deg
        max_terms = 2**(in_size + 1) if max_terms == -1 else max_terms

        indices = list(range(in_size))
        all_monomials = [[]]

        # create all monomials with regards to the maximal degree
        for deg in range(1, max_deg + 1):
            all_monomials.extend([list(monomials) for monomials in itertools.combinations(indices, deg)])

        all_options = []
        max_terms = min(max_terms, len(all_monomials))

        # create all combinations of monomials with regards to the maximal number of terms/monomials
        for terms in range(1, max_terms + 1):
            all_options.extend([list(option) for option in itertools.combinations(all_monomials, terms)])
        
        self.vars = in_size
        self.functions = []
        for option in all_options:
            self.functions.append(BooleanFunction(in_size, 1, [option]))


if __name__ == "__main__":
    in_size = int(sys.argv[1])
    max_deg = int(sys.argv[2]) if len(sys.argv) > 2 else -1
    max_terms = int(sys.argv[3]) if len(sys.argv) > 3 else -1
    max_vars = int(sys.argv[4]) if len(sys.argv) > 4 else -1
    all_functions = AllBooleanFunctions(in_size, max_deg, max_terms, max_vars)
