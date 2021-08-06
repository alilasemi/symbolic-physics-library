import sympy as sp

# This is a class to be inherited from for any mathematical expression that
# needs to be created. Each term in an expression must be defined as a
# subexpression, which is eventually plugged in.
class Expression:

    # Symbol representing this expression
    symbol = None

    # Sympy expression
    expression = None

    def __init__(self, expression):
        self.expression = expression

    # TODO: Robust error handling and messages. Right now, it just assumes you
    # pass in either a dict or a scalar expression. It also fails silently,
    # which is bad in this context. Not user friendly.
    def plug_in(self, symbol, expression):

        # For dictionaries, plug in one by one
        if isinstance(expression, dict):
            # Loop over dictionary
            for i, expr in expression.items():
                self.expression = self.expression.subs(symbol[i], expr)

        # Else, just assume it's a scalar (for now)
        else:
            self.expression = self.expression.subs(symbol, expression)

        return self

    def doit(self):
        self.expression = self.expression.doit()
        return self
