#!/usr/bin/env python3.7

import bitarray
import itertools
import sys


class IterativeMethod:
    def __init__(self, blocklen=128, max_monomials=5):
        self.blocklen      = blocklen
        self.max_monomials = max_monomials
        self.terms         = None
        self.output        = None
        self.saved_combs   = []
        self.saved_arrays  = []
        self.monomials     = []
        self.approximated  = False
        self.rate          = 0.0
        self.abs_rate      = 0


    # TODO different blocklen, maybe we do not need integers
    # biggest results - [[13, 77], [14, 88], [46, 116], [46, 48], [54, 109]], 503144
    def load_data(self, filename):
        # Add also the true vector
        self.terms  = [bitarray.bitarray(0) for i in range(self.blocklen + 1)]
        self.output = bitarray.bitarray(0)

        with open(filename, "rb") as f:
            data = f.read(16)

            in_number = bitarray.bitarray()
            in_number.frombytes(data)

            while True:
                data = f.read(16)
                if data == b'':
                    break

                next_number = bitarray.bitarray()
                next_number.frombytes(data)
                out_number = next_number[0]

                self.terms[0].append(True)
                for term_index in range(self.blocklen):
                    self.terms[1 + term_index].append(in_number[term_index])
                self.output.append(out_number)

                in_number = bitarray.bitarray(next_number)


    # find index of combination in self.saved_combs
    def __find_saved(self, combination):
        for i in range(len(self.saved_combs)):
            if self.saved_combs[i] == combination:
                return i
        return -1


    # add monomial to the function and xor the output with the monomial
    def __apply_monomial(self, closest_monomial):
        # we can not enhance the approximation more
        if closest_monomial[1] < self.abs_rate:
            self.approximated = True
            return

        # add monomial to the function
        if closest_monomial[0] in self.monomials:
            self.monomials.remove(closest_monomial[0])
        else:
            self.monomials.append(closest_monomial[0])

        # get bits of closest monomial
        if len(closest_monomial[0]) == 0:
            index = closest_monomial[0][0]
            closest_monomial_bits = self.terms[index]
        else:
            # it must be here
            index = self.__find_saved(closest_monomial[0])
            closest_monomial_bits = self.saved_arrays[index]

        # modify the output
        self.output ^= closest_monomial_bits

        # set new approximation rate
        self.abs_rate = closest_monomial[1]

        # the function has maximal number of monomials
        if len(self.monomials) >= self.max_monomials:
            self.approximated = True


    # pick n terms with the smallest hamming distance from the output
    def __pick_n_closest(self, n=10):
        term_counts = []

        # count hamming distance from output
        for i in range(len(self.terms)):
            term_xor_output = self.terms[i] ^ self.output
            term_counts.append(([i], term_xor_output.count(False)))

        # sort terms by smallest hamming distance from output
        term_counts.sort(key = lambda term : term[1], reverse=True)

        return term_counts[:n]


    # pick monomial with the smallest hamming distance from the output
    # chooses only from the n best terms picked by pick_n_closest function
    def pick_closest_combined(self):
        closest_monomials = self.__pick_n_closest()
        indices = [index[0] for index, count in closest_monomials]

        # count hamming distance of degree=2 combinations of closest terms
        combinations = [list(comb) for comb in itertools.combinations(indices, 2)]
        for combination in combinations:
            # combinations can be unordered
            combination.sort()

            # get the monomial from the class
            combination_index = self.__find_saved(combination)
            if combination_index != -1:
                terms_and = self.saved_arrays[combination_index]
            else:
                # FIXME: this currently works only for degree=2
                terms_and = self.terms[combination[0]] & self.terms[combination[1]]

            # save monomial
            self.saved_combs.append(combination)
            self.saved_arrays.append(terms_and)

            terms_and_xor_output = terms_and ^ self.output                
            closest_monomials.append((combination, terms_and_xor_output.count(False)))

        # sort monomials by smallest hamming distance from output
        closest_monomials.sort(key = lambda terms : terms[1], reverse=True)

        # find the closest monomial which enhances the approximation and change the output array
        self.__apply_monomial(closest_monomials[0])
