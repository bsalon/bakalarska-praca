#!/usr/bin/env python3.7

import fileinput
import sys


def char_is_digit(char):
    return ord(char) >= ord('0') and ord(char) <= ord('9')


def integer_substring(string, index):
    number = ""
    while char_is_digit(string[index]):
        number += string[index]
        index += 1
    return number


def transform_to_latex(string):
    transformed = "$ "
    length = len(string)
    index = 0

    while index < length:
        if string[index] in ('(', ')', ' '):
            index += 1
            continue
        
        elif string[index] == '^':
            transformed += " + "
        
        elif string[index] == '&' and string[index + 1] == '&':
            index += 1
            continue
        
        elif string[index] == 'v' and string[index + 1] == '[':
            number = integer_substring(string, index + 2)
            index += 2 + len(number)
            transformed += "v_{" + number + "}"
        
        elif string[index] == 'x' and string[index + 1] == '[':
            number = integer_substring(string, index + 2)
            index += 2 + len(number)
            transformed += "x_{" + number + "}" 

        index += 1

    return transformed + " $"


# we use stdin instead of arguments because the boolean function can be long
if __name__ == "__main__":
    to_transform = ""
    for line in fileinput.input():
        to_transform += line
    print(transform_to_latex(to_transform))
