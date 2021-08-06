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
        if expression.isinstance(dict):
            # Loop over dictionary
            for i, expr in expression.items():
                self.expression = self.expression.subs(symbol[i], expr)

        # Else, just assume it's a scalar (for now)
        else:
            self.expression = self.expression.subs(symbol, expression)

        # TODO: Remove the graveyard below
        ## Plug in sizes first. This is so that sums and products are expanded.
        #for symbol, size in self.sizes.items():
        #    # Plug in size and expand
        #    self.expression = self.expression.subs(symbol, size).doit()
        #    # Plug in to the symbol (useful for indexed symbols)
        #    self.symbol = self.symbol.subs(symbol, size)

        ## Next, data is plugged in. This way, terms that become zero are
        ## cancelled out as soon as possible.
        #for symbol, value in self.data.items():
        #    # Plug in arrays elementwise
        #    if isinstance(value, np.ndarray):
        #        for i in range(value.size):
        #            # Get the index of this element
        #            index = np.unravel_index(i, value.shape)
        #            # Plug in
        #            self.expression = self.expression.subs(symbol[index],
        #                    value[index])
        #    # Plug in everything else directly
        #    else:
        #        self.expression = self.expression.subs(symbol, value)

        ## Loop over subexpressions
        #for sub_expression in self.sub_expressions:
        #    # If it's a list of expressions, plug in each one
        #    if isinstance(sub_expression, list):
        #        self.expression = self.expression.subs([(
        #            sub_expression[i].symbol[i], sub_expression[i].expression)
        #            for i in range(len(sub_expression))])
        #    # Otherwise, just plug it in
        #    else:
        #        self.expression = self.expression.subs(sub_expression.symbol,
        #                sub_expression.expression)
