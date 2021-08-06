import numpy as np
import os
import sys
import urllib.request
import physics.constants as constants


# Professor Joseph Shephard at Caltech provides the NASA Glenn thermodynamic
# data files on his website, so these are downloaded from there
nasa9_url = 'https://shepherd.caltech.edu/EDL/PublicResources/sdt/SDToolbox/cti/NASA9/nasa9.dat'
nasa7_url = 'https://shepherd.caltech.edu/EDL/PublicResources/sdt/SDToolbox/cti/NASA7/nasa7.dat'

# Location to download data files
root_dir = sys.path[0]
nasa9_file_name = root_dir + '/physics/thermodynamics/data/nasa9.dat'
nasa7_file_name = root_dir + '/physics/thermodynamics/data/nasa7.dat'

class ThermochemicalData:

    def __init__(self, model):
        self.model = model
        # Set files according to the model chosen
        if model == 'NASA9':
            self.file_name = nasa9_file_name
            self.data_url = nasa9_url
        elif model == 'NASA7':
            self.file_name = nasa7_file_name
            self.data_url = nasa7_url
        else: print('Invalid model chosen')

        # Download the data, if needed
        self.get_data()
        # Read in the data file
        self.read_data()

    def __getitem__(self, key):
        return self.data[key]

    def get_data(self):
        # If the data file doesn't exist yet
        if not os.path.isfile(self.file_name):
            urllib.request.urlretrieve(self.data_url, self.file_name)

    def read_data(self):
        '''Load thermo file into NASA9 objects'''
        self.data = {}
        # Open the thermo file
        with open(self.file_name, 'r') as thermo_file:
            # Keep reading, until it hits the word "thermo" (case
            # insensitive)
            for line in thermo_file:
                thermo_file.readline()
                if line.lower().startswith('thermo'): break

            # Keep reading until there is no more data left
            line = thermo_file.readline()
            while not line.lower().startswith('end'):
                # Start recording the text for this species
                text = []
                text.append(line)
                # Keep going until the next species
                line = thermo_file.readline()
                while line.startswith((' ', '-')):
                    text.append(line)
                    line = thermo_file.readline()
                # Get species name from the first line
                species = text[0].split()[0]
                # Create NASA9 object and store
                self.data[species] = NASA9(text)


class NASA9:
    '''
    species_name
    M
    dHf
    n_ranges
    n_coeffs
    temperatures
    temperature_ranges
    a
    text
    '''

    n_coeffs = 9;

    def __init__(self, text):
        # Store text from NASA9 file for this species
        self.text = text
        # Parse text
        self.parse_text()

    # Parse text from a NASA9 file into data structures.
    def parse_text(self):
        # The first line contains the species name
        self.species_name = self.text[0].split()[0]
        # The second line contains the number of temperature ranges, the
        # molar mass, and the heat of formation
        self.n_ranges = int(self.text[1].split()[0])
        self.M        = float(self.text[1].split()[-2]) / 1000. # Convert to kg
        self.dHf      = float(self.text[1].split()[-1])
        self.temperatures = np.empty(self.n_ranges + 1)
        self.a = np.empty((self.n_ranges, self.n_coeffs))
        # Loop over temperature ranges
        for i in range(self.n_ranges):
            # Get text for this range
            range_text = self.text[2 + i*3 : 5 + i*3]
            # The first line in the range contains the temperatures
            self.temperatures[i] = range_text[0][:11]
            self.temperatures[i + 1] = range_text[0][11:22]
            # The second line in the range contains the first five coefficients.
            # Must convert D to e for exponentials
            line = range_text[1].replace('D', 'e')
            for num in range(5):
                self.a[i, num] = line[num*16 : (num + 1)*16]
            # The third line in the range contains the last four coefficients.
            # Must convert D to e for exponentials
            line = range_text[2].replace('D', 'e')
            for num in range(2):
                self.a[i, 5 + num] = line[num*16 : (num + 1)*16]
            for num in range(2):
                self.a[i, 7 + num] = line[(num + 3)*16 : (num + 4)*16]
        # Species mass
        self.m = self.M / constants.N_A
