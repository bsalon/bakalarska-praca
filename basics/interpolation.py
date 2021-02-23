#!/usr/bin/env python3.7

import random_function

import sys


# source: http://www.selmer.uib.no/odbf/help/ttanf.pdf

# Let C = [c0 c1 ... c2ⁿ-1] be the coefficient vector of
# the polynomial representing the Boolean function f.
# If C = [0 0 0 0 1 0 0 1] then the ANF is x₀ + x₀x₁x₂.

#   THEOREM
# Let f be the truth table of an n-variable Boolean function.
# Let be as defined above. Then
#                   C = fA_n
# where
#                   A_n = [ 1 1 ] (^ xor n)
#                         [ 0 1 ]
# That is A_n is the nth tensor power mod 2 of the matrix
# [ 1 1 ]
# [ 0 1 ], or in other notation,
#                   A_n = [ A_n-1 A_n-1 ]
#                         [ 0     A_n-1 ]   and A₀ = [ 1 ]

# The above theorem helps us to convert from truth table to
# AND and vice versa in almost 2²^ⁿ binary operations. The
# following algorithm reduces the conversion to only O(n2²)
# operations:

# ALGORITHM: TT -> ANF
# Input: TT of a Boolean function f
# Output: The coefficient vector of the ANF of f
# For 0 <= k <= n, define f_(k,a) ∈ F₂k,
#     where 0 <= a <= 2^(n-k) - 1
# 1. Set f_0,a = f(a) for 0 <= a <= 2ⁿ - 1
# 2. for k = 0 to n - 1 do
#        for b = 0 to 2^(n - k - 1) - 1 do
#            f_(k+1, b) = [f_(k, 2b) f_(k, 2b + 1) + f_(k, 2b)]
# 3. C = f_(n, 0)

# The following example illustrates what the above algorithm does.
# Let us find the ANF of the Boolean function represented by the
# truth table in Table 2. We have f = [0 1 1 0 0 1 0 1]. Looking
# at the 1's positions in C we see that the ANF of f
# is x₂ + x₁ + x₀x₁.
# |    f     | 0 | 1 | 1 | 0 | 0 | 1 | 0 | 1 |
# |  k = 0   | 0 | 1 | 1 | 1 | 0 | 1 | 0 | 1 |
# |  k = 1   | 0 | 1 | 1 | 0 | 0 | 1 | 0 | 0 |
# |  k = 2   | 0 | 1 | 1 | 0 | 0 | 0 | 1 | 0 |
# | C = f₃,₀ | 0 | 1 | 1 | 0 | 0 | 0 | 1 | 0 |
#   Table 2: Converting TT to ANF algorithm

# The above algorithm can also be used to convert from ANF to TT
# by changing the input from TT to C (the coefficients) of the
# ANF of f.


def get_truth_table(function, degree):
    return [function.eval_f(value) for value in range(2 ** degree)]


def update_values(the_list, k, b):
    power = 2 ** k
    index = power + (b * (2 << k))
    for i in range(power):
        the_list[index + i] ^= the_list[index + i - power]


def get_anf_from_tt(results, degree):
    coeficients = results

    for k in range(degree):
        for b in range(2 ** (degree - k - 1)):
            update_values(coeficients, k, b)
    
    return coeficients


def get_one_bit_list(the_list, bit):
    return [i & (2 ** bit) for i in the_list]



# mandatory command line arguments:
#  sys.argv[1]: number of input-output bits of the function
#  sys.argv[2]: seed

# TODO: work with more general functions, not only what I create
# TODO: better manipulation with command line arguments
# TODO: create the same function

if __name__ == "__main__":
    if len(sys.argv) > 2:
        size = int(sys.argv[1])
        seed = int(sys.argv[2])

        function = random_function.RandomFunction(size, size, seed)
        print(function.x, end="\n\n")

        indices = list(range(2 ** size))
        print("The indices: ", end="")
        print(*indices)

        tt = get_truth_table(function, size)
        print("Truth table: ", end="")
        print(*tt)

        anf = get_anf_from_tt(tt, size)
        print("ANF of TT:   ", end="")
        print(*anf)

        print()
        bits = []
        for bit in range(size):
            bits.append(get_one_bit_list(anf, bit))
            print(str(bit) + "-bit ANF:   ", end="")
            print(*bits[bit])
    pass
