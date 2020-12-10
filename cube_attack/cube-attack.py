#!/usr/bin/env python3.7

import sys
import itertools
import subprocess
import numpy as np


def pick_indices(whole_list, indices):
    return [whole_list[index] for index in indices]


def pick_keys(dictionary, values):
    return [key for key, value in dictionary.items() if value in values]


def find_independent_matrix(vectors):
    degree = len(vectors[0])
    possibilities = [list(combination) for combination in itertools.combinations(vectors, degree)]
    independent = None

    for possibility in possibilities:
        # count determinant of the matrix
        if abs(np.linalg.det(np.array(possibility)) - 0.0) > 1e-09:
            independent = possibility
            break

    if independent == None:
        return ([], [])

    # return also indices on which were the chosen vectors
    indices = []
    for index in range(len(vectors)):
        if vectors[index] in independent:
            indices += [index]

    return (independent, indices)


def bin_to_dec(bin_list):
    count = 0
    for index in range(len(bin_list)):
        if bin_list[index]:
            count += int(bin_list[index]) << index

    return count


class CubeAttack:
    
    def __init__(self, cipher_filename, public_degree):
        self.cipher = cipher_filename
        self.count_subterms(public_degree)
        self.superpolys = {}
        self.free_terms = []
        self.cipher_runs = {}


    def count_subterms(self, public_degree):
        self.subterms = []
        indices = range(public_degree)
        for i in range(1, public_degree + 1):
            self.subterms += [list(subterm) for subterm in itertools.combinations(indices, i)]


    def preprocessing_phase(self, secret_degree, limit=None):
        # maximal number of unique superpolys
        if limit == None:
            limit = len(self.subterms) + 1
        
        public_input = 0
        for subterm in self.subterms:
            bins = [list(p) for p in itertools.product([0, 1], repeat=len(subterm))]
            
            # summation cube values
            subterm_inputs = []
            for binary_option in bins:
                subterm_inputs.append(sum([b * (2 ** i) for (i, b) in zip(subterm, binary_option)]))
            
            print("subterm: " + str(subterm), file=sys.stderr)
         
            cipher_evals = self.constant_test(subterm_inputs, secret_degree)
            if cipher_evals == []:
                print("X: superpoly is constant", file=sys.stderr, end="\n\n")
                continue
            
            if not self.linearity_test(cipher_evals, secret_degree):
                print("X: superpoly is not linear", file=sys.stderr, end="\n\n")
                # self.delete_subterms(subterm)
                continue

            print("O: superpoly is linear", file=sys.stderr)

            if self.superpoly_coefficients(cipher_evals, subterm_inputs):
                print("x" + str(self.superpolys[tuple(subterm_inputs)]), file=sys.stderr)
                limit -= 1
            print(file=sys.stderr)

            if limit <= 0:
                return



    def preprocessing_cipher_run(self, subterm_inputs, secret_key):
        output = False

        for subterm_input in subterm_inputs:

            # public input, secret key
            options = (str(subterm_input), str(secret_key))
            
            # we already have the result
            if options in self.cipher_runs:
                cipher_bit = self.cipher_runs[options]
            # evaluate cipher on options
            else:
                cipher_process = subprocess.run(["./" + str(self.cipher), options[0], options[1]])
                cipher_bit = cipher_process.returncode
                self.cipher_runs[options] = cipher_bit

            output ^= True if cipher_bit == 1 else False

        return output



    def constant_test(self, subterm_inputs, secret_degree):
        is_constant = True
        zero_key = self.preprocessing_cipher_run(subterm_inputs, 0)
        
        # list of cipher evaluations on all possible secret keys
        cipher_evals = [zero_key]
        
        for secret_key in range(1, 2 ** secret_degree):
            current_key = self.preprocessing_cipher_run(subterm_inputs, secret_key)
            if zero_key != current_key:
                is_constant = False
            cipher_evals.append(current_key)

        return [] if is_constant else cipher_evals



    # test cipher(0) + cipher(x) + cipher(y) == cipher(x + y), for all possible x,y
    def linearity_test(self, cipher_evals, secret_degree):
        is_linear = True
        for x in range(2 ** secret_degree):
            if not is_linear:
                break
            for y in range(x + 1, 2 ** secret_degree):
                if (cipher_evals[0] ^ cipher_evals[x] ^ cipher_evals[y]) != cipher_evals[x ^ y]:
                    is_linear = False
                    break
        
        return is_linear



    def superpoly_coefficients(self, cipher_evals, subterm):
        free_term = int(cipher_evals[0])

        # get superpoly coefficients
        superpoly = []
        secret_key = 1
        while secret_key < len(cipher_evals):
            superpoly += [1] if free_term != int(cipher_evals[secret_key]) else [0]
            secret_key *= 2
        
        # add new superpoly to dictionary
        if superpoly not in self.superpolys.values():
            self.free_terms += [free_term]
            self.superpolys[tuple(subterm)] = superpoly
            return True

        # superpoly has already been found
        return False



    def online_phase(self):
        public_runs = {}
        print("free terms: " + str(self.free_terms))
        print("superpolys: " + str(list(self.superpolys.values())), end="\n\n")
        print("===  looking for independent superpolys  ===")
        
        matrix, indices = find_independent_matrix(list(self.superpolys.values()))
        if matrix == []:
            print("Not enough lineary independent superpolys to attack the cipher")
            return
        
        terms = pick_indices(self.free_terms, indices)
        superpolys_public_inputs = pick_keys(self.superpolys, matrix)
        print("free terms: " + str(terms))
        print("superpolys: " + str(matrix), end="\n\n")
        
        # online phase of the attack
        index = 0
        for public_inputs in superpolys_public_inputs:
            output = False
            for public_input in public_inputs:
                
                if public_input in public_runs:
                    cipher_bit = public_runs[public_input]
                else:
                    cipher_process = subprocess.run(["./" + str(self.cipher), str(public_input)])
                    cipher_bit = cipher_process.returncode
                    public_runs[public_input] = cipher_bit

                output ^= True if cipher_bit == 1 else False

            terms[index] += output
            index += 1

        result = np.linalg.solve(np.array(matrix), np.array(terms))
        secret_key = [val % 2 for val in result]
        print("THE KEY: " + str(secret_key))
        print("number:  " + str(bin_to_dec(secret_key)))


# For now, accept only public_degree == secret_degree
if __name__ == "__main__":
    cipher = sys.argv[1]
    degree = int(sys.argv[2])
    cube_attack = CubeAttack(cipher, degree)
    cube_attack.preprocessing_phase(degree)
    cube_attack.online_phase()
