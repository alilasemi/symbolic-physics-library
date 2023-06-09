import sympy as sp
from sympy.codegen.ast import Assignment

import symbols as syms

class SourceCode:

    tab = '    '

    def __init__(self, file_name, marker, inputs_line, expressions, pointer=False):
        self.file_name = file_name + '.c'
        self.marker = marker

        # If there are a list of names, connect them all for a name
        if isinstance(self.marker, list):
            self.name = ''
            for marker in self.marker:
                self.name += marker + '_and_'
            self.name = self.name[:-5]
        else:
            self.name = self.marker

        self.inputs_line = inputs_line
        self.expressions = expressions
        self.pointer = pointer
        # Read in base file
        #self.read_base_file()
        self.create_base_file()
        # Generate expressions
        self.generate_expressions()
        # Write to file
        self.write()

    def read_base_file(self):
        self.text = []
        # Get name of base file
        name, extension = self.file_name.split('.')
        base_file_name = name + '_base.' + extension
        with open(base_file_name, 'r') as base_file:
            for line in base_file:
                self.text.append(line)

    def create_base_file(self):
        self.text = []
        self.text.append('#include <math.h>\n')
        self.text.append('\n')
        self.text.append(f'void compute_{self.name}({self.inputs_line}) {{\n')
        self.text.append('\n')
        self.text.append(f'#pragma {self.name}')
        self.text.append('\n')
        self.text.append('}')

    def generate_expressions(self):
        # Loop over lines in source code
        for i, line in enumerate(self.text):
            # If this is the marker for wdot
            if line.startswith(f'#pragma {self.name}'):
                # Remove marker
                self.text.pop(i)
                # Loop over expressions
                text = []
                # Try iterating over elements
                try:
                    for j in range(len(self.expressions)):
                        # Generate code
                        code = self.code(self.expressions[j], j)
                        text.append(self.format(code))
                        text.append('\n')
                # If it doesn't have elements, then just generate the one
                # expression
                except TypeError:
                    # Generate code
                    code = self.code(self.expressions)
                    text.append(self.format(code))
                # Insert text
                self.text[i:i] = text
                break

    def code(self, expr, idx=None):
        prefix = '*' if self.pointer else ''
        suffix = ';'
        # If there are multiple markers, use those as the variable name, instead
        # of treating it as an array
        if isinstance(self.marker, list):
            variable_name = self.marker[idx]
        else:
            indexing = f'[{idx}]' if idx is not None else ''
            variable_name = self.marker + indexing
        return (prefix + variable_name + ' = ' +
                sp.ccode(expr.expression, contract=False) + suffix)

    def format(self, text, indent=1):
        white_space = indent * self.tab
        return white_space + text + '\n'

    def write(self):
        with open(self.file_name, 'w') as out_file:
            for line in self.text:
                out_file.write(line)
