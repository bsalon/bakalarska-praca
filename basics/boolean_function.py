from arybo.lib import MBA


class BooleanFunction:
    # {0,1}^in_size -> {0,1}^out_size
    def __init__(self, in_size, out_size, monomials_lists):
        mba = MBA(1)

        self.v = [mba.var("x0")]
        for i in range(1, in_size):
            self.v.append(mba.var("x" + str(i)))

        self.x = []
        for i in range(out_size):
            self.x.append(self._create_poly(monomials_lists[i]))


    # viac to generalizuj cele
    # bins = [list(p) for p in itertools.product([0, 1], repeat=in_size)]

    # returns True/False wrapped in pytanque.Expr
    def pytanque_expr_const(self, boolean):
        mba = MBA(1)
        x = mba.var("x")
        x = x | 1 if boolean else x & 0
        return x[0]


    # create polynomial in ANF from values monomials
    def _create_poly(self, monomials):
        mba = MBA(1)
        polynom = mba.var("")

        polynom_expr = self.pytanque_expr_const(len(monomials) == 0 or monomials[0] == [-1])

        for monomial in monomials:
            polynom_expr = polynom_expr + self._create_expr(monomial)

        polynom = (polynom | 0x1) & polynom_expr

        return polynom


    # create monomial in ANF from chosen indices
    def _create_expr(self, indices):
        if len(indices) == 0:
            return self.pytanque_expr_const(True)

        expression = self.pytanque_expr_const(True)
        for index in indices:
            expression *= self.v[index][0]

        return expression


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


    # transform the function string into more readable one
    def _create_str_function(self, string):
        i = 0
        while string[i] != '\n':
            i += 1
        i += 1

        inside = False
        new_string = ""

        while string[i] != '\n':
            if string[i] == 'x':
                inside = True

            elif inside and string[i] == '0' and (string[i + 1] in (' ',')','\n')):
                inside = False
                i += 1
                continue

            elif string[i] in (' ', "*", ")", "("):
                if (string[i - 1] == '+' or string[i + 1] == '+'):
                    new_string += ' '
                i += 1
                continue

            new_string += string[i]
            i += 1

        return new_string


    def __str__(self):
        string = ""
        for bit in range(len(self.x)):
            string += self._create_str_function(self.x[bit].__str__()) + '\n'

        return string[:-1]
