import numpy as np


class KineticsData():

    def __init__(self, kinetics_file_path):
        self.kinetics_file_path = kinetics_file_path

        # Read in kinetics read
        self.read()

    def read(self):

        # Initialize
        self.three_body = []
        self.epsilon = {}

        # Open kinetics  file
        with open(self.kinetics_file_path, 'r') as kinetics_file:
            # Keep reading until there are no more lines
            line = kinetics_file.readline()
            while line:
                # Get the species
                if line.startswith('Species Order:'):
                    species_string = line.replace('Species Order:','')
                    species_string = species_string.replace(' ', '').strip()
                    self.species = species_string.split(',')
                    self.ns = len(self.species)

                # Get the number of reactions
                if line.startswith('nr'):
                    line = kinetics_file.readline()
                    self.nr = int(line.strip().split(',')[0])
                    # Allocate arrays for stoichiometric coefficients
                    self.alpha = np.empty((self.ns, self.nr), dtype=int)
                    self.beta = np.empty_like(self.alpha)
                    # Allocate arrays for Arrhenius parameters
                    self.C = np.empty(self.nr)
                    self.eta = np.empty(self.nr)
                    self.theta = np.empty(self.nr)

                # Get the reactions
                if line.startswith('Reactant/product species'):
                    line = kinetics_file.readline().strip()
                    # Skip the line breaks
                    while not line:
                        line = kinetics_file.readline().strip()
                    # The reactions use p for positive ions, not +
                    species_with_p = [s.replace('+', 'p') for s in
                            self.species]
                    # For each reaction
                    for i in range(self.nr):
                        # Clean beginning
                        line = line.split(f'#{i+1}')[-1].strip()
                        # Reactants
                        reactants = line.split('<')[0].strip().replace(' ', '')
                        reactants = reactants.split('+')
                        self.alpha[:, i] = [reactants.count(s) for s in
                                species_with_p]
                        # Products
                        products = line.split('>')[-1].strip().replace(' ', '')
                        products = products.split('+')
                        self.beta[:, i] = [products.count(s) for s in
                                species_with_p]
                        # Mark whether or not this reaction has a third body
                        self.three_body.append('M' in reactants and 'M' in
                                products)
                        # Read the next line
                        line = kinetics_file.readline()
                    # Net stoichiometric coefficient
                    self.nu = np.sum(self.beta - self.alpha, axis=0)
                    # Reaction indices of third-body reactions
                    three_body_indices = np.where(self.three_body)[0]

                # Get the collision constants
                if line.startswith('Collision constants'):
                    line = kinetics_file.readline().strip()
                    # Skip the line breaks
                    while not line:
                        line = kinetics_file.readline().strip()
                    # For each third-body reaction
                    for three_body_index in three_body_indices:
                        # Get collision constants
                        line = line.strip().replace(' ', '').replace('d', 'e').split(',')
                        self.epsilon[three_body_index] = np.array([float(string)
                            for string in line])
                        # Read the next line
                        line = kinetics_file.readline()

                # Get the Arrhenius parameters
                if line.startswith('Arrhenius parameters'):
                    line = kinetics_file.readline().strip()
                    # Skip the line breaks
                    while not line:
                        line = kinetics_file.readline().strip()
                    # For each reaction
                    for i in range(self.nr):
                        line = line.strip().split('#')[0].replace(' ', '')
                        line = line.replace('d', 'e').split(',')
                        self.C[i] = float(line[0])
                        self.eta[i] = float(line[1])
                        self.theta[i] = float(line[2])
                        # Read the next line
                        line = kinetics_file.readline()

                    # Normalize the collision constants using the Arrhenius C
                    for i in self.epsilon.keys():
                        self.epsilon[i] /= self.C[i]

                    # The units on C need to be in SI units:
                    #     m^(3 a) * mol^(1 - a) * s,
                    # where a is the sum of reactant stoichimetric coeffs (including M).
                    # The DPLR input file assumes kilomoles, so this needs to be converted.
                    kmol_to_mol = 1e3
                    a = np.sum(self.alpha, axis=0) + self.three_body
                    conversion = kmol_to_mol ** (1 - a)
                    self.C /= conversion

                # Read the next line
                line = kinetics_file.readline()

            # Change electrons from 'e' to 'e-'
            self.species = ['e-' if s == 'e' else s for s in self.species]
