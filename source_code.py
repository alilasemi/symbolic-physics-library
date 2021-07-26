import sympy as sp
from sympy.codegen.ast import Assignment

import symbols as syms

class SourceCode:

    tab = '    '

    def __init__(self, file_name, marker, equations, pointer=False):
        self.file_name = file_name
        self.marker = marker
        self.equations = equations
        self.pointer = pointer
        # Read in base file
        self.read_base_file()
        # Generate equations
        self.generate_equations()
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

    def generate_equations(self):
        # Loop over lines in source code
        for i, line in enumerate(self.text):
            # If this is the marker for wdot
            if line.startswith(f'#pragma {self.marker}'):
                # Remove marker
                self.text.pop(i)
                # Loop over equations
                text = []
                for j in range(len(self.equations)):
                    # Generate code
                    code = self.code(self.equations[j])
                    text.append(self.format(code))
                # Insert text
                self.text[i:i] = text
                break

    def code(self, expr):
        prefix = '*' if self.pointer else ''
        return prefix + sp.ccode(expr, contract=False)

    def format(self, text, indent=1):
        white_space = indent * self.tab
        return white_space + text + '\n'

    def write(self):
        with open(self.file_name, 'w') as out_file:
            for line in self.text:
                out_file.write(line)
