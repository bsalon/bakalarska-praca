#!/usr/bin/env python3.7

# Create random boolean function, with chosen public and secret degree
# The created function is in C language form:
#   bool bool_function = ...;

import random
import itertools
import sys


def create_rand_polynomial(public_degree, secret_degree):
    public_indices = range(public_degree)
    secret_indices = range(secret_degree)
    polynomial = ""

    for i in range(public_degree + 1):
        # all possible public polynomials of degree i
        i_public_deg_polys = [list(comb) for comb in itertools.combinations(public_indices, i)]
        
        for j in range(secret_degree + 1):
            # all possible secret polynomials of degree j
            j_secret_deg_polys = [list(comb) for comb in itertools.combinations(secret_indices, j)]

            for public_poly in i_public_deg_polys:
                for secret_poly in j_secret_deg_polys:
                    # with (1.0 / (2 ** ((2 * j * j) + 1)))
                    # polynomials with lower secret degree have bigger probability to be chosen
                    # pick_me = random.random() < (1.0 / (2 ** ((2 * j * j) + 1)))
                    pick_me = random.random() < (1 / 2.0)
                    if pick_me:
                        polynomial += " ^ " + poly_str_representation(public_poly, secret_poly)

    return polynomial.replace(" ^ ", "", 1)  + ";"



def poly_str_representation(public_poly, secret_poly):
    if public_poly == [] and secret_poly == []:
        return "(true)"
    polynomial = "("
    for term in public_poly:
        polynomial += " && v[" + str(term) + "]"
    for term in secret_poly:
        polynomial += " && x[" + str(term) + "]"
    polynomial += ")"
    return polynomial.replace(" && ", "", 1)
    

# For now, accept only public_degree == secret_degree
if __name__ == "__main__":
    degree = int(sys.argv[1])
    print("\tbool bool_function = " + create_rand_polynomial(degree, degree))
