import sympy as sp

# This is a class to be inherited from for any mathematical expression that
# needs to be created. Each term in an expression must be defined as a
# subexpression, which is eventually plugged in.
class Expression:

    # Symbol representing this expression
    symbol = None

    # Sympy expression
    expression = None

    def __init__(self, expression, *args):
        # If the expression given is actually a function, then use the function
        # to initialize, with any given arguments
        if callable(expression):
            args_exprs = [arg.expression for arg in args]
            self.expression = expression(*args_exprs)
        # Otherwise, it is expected that the given expression is a Sympy
        # expression
        else:
            self.expression = expression

    def __repr__(self):
        return self.expression.__repr__()

    # TODO: Robust error handling and messages. Right now, it just assumes you
    # pass in either a dict or a scalar expression. It also fails silently,
    # which is bad in this context. Not user friendly.
    def plug_in(self, symbol, expression):

        # For dictionaries, plug in one by one
        if isinstance(expression, dict):
            # Loop over dictionary
            for i, expr in expression.items():
                self.expression = self.expression.subs(symbol[i], expr)

        # For other expressions, plug in the expression
        if isinstance(expression, Expression):
            self.expression = self.expression.subs(symbol,
                    expression.expression)

        # Else, just assume it's a scalar (for now)
        else:
            self.expression = self.expression.subs(symbol, expression)

        return self

    def doit(self):
        self.expression = self.expression.doit()
        return self

    def simplify(self):
        self.expression = self.expression.simplify()
        return self
