import numpy as np


# This data comes from the IUPAC Commission on Isotopic Abundances and Atomic
# Weights. The data below was taken from the table at:
#     https://iupac.qmul.ac.uk/AtWt/

atomic_molar_mass = np.array([
[1,   'H',   'Hydrogen',    1.008,],
[2,   'He',  'Helium',  4.002602,],
[3,   'Li',  'Lithium', 6.94,],
[4,   'Be',  'Beryllium',   9.0121831,],
[5,   'B',   'Boron',   10.81,],
[6,   'C',   'Carbon',  12.011,],
[7,   'N',   'Nitrogen',    14.007,],
[8,   'O',   'Oxygen',  15.999,],
[9,   'F',   'Fluorine',    18.998403163,],
[10,  'Ne',  'Neon',    20.1797,],
[11,  'Na',  'Sodium',  22.98976928,],
[12,  'Mg',  'Magnesium',   24.305,],
[13,  'Al',  'Aluminium',   26.9815384,],
[14,  'Si',  'Silicon', 28.085,],
[15,  'P',   'Phosphorus',  30.973761998,],
[16,  'S',   'Sulfur',  32.06,],
[17,  'Cl',  'Chlorine',    35.45,],
[18,  'Ar',  'Argon',   39.948,],
[19,  'K',   'Potassium',   39.0983,],
[20,  'Ca',  'Calcium', 40.078,],
[21,  'Sc',  'Scandium',    44.955908,],
[22,  'Ti',  'Titanium',    47.867,],
[23,  'V',   'Vanadium',    50.9415,],
[24,  'Cr',  'Chromium',    51.9961,],
[25,  'Mn',  'Manganese',   54.938043,],
[26,  'Fe',  'Iron',    55.845,],
[27,  'Co',  'Cobalt',  58.933194,],
[28,  'Ni',  'Nickel',  58.6934,],
[29,  'Cu',  'Copper',  63.546,],
[30,  'Zn',  'Zinc',    65.38,],
[31,  'Ga',  'Gallium', 69.723,],
[32,  'Ge',  'Germanium',   72.630,],
[33,  'As',  'Arsenic', 74.921595,],
[34,  'Se',  'Selenium',    78.971,],
[35,  'Br',  'Bromine', 79.904,],
[36,  'Kr',  'Krypton', 83.798,],
[37,  'Rb',  'Rubidium',    85.4678,],
[38,  'Sr',  'Strontium',   87.62,],
[39,  'Y',   'Yttrium', 88.90584,],
[40,  'Zr',  'Zirconium',   91.224,],
[41,  'Nb',  'Niobium', 92.90637,],
[42,  'Mo',  'Molybdenum',  95.95,],
[43,  'Tc',  'Technetium',  97,],
[44,  'Ru',  'Ruthenium',   101.07,],
[45,  'Rh',  'Rhodium', 102.90549,],
[46,  'Pd',  'Palladium',   106.42,],
[47,  'Ag',  'Silver',  107.8682,],
[48,  'Cd',  'Cadmium', 112.414,],
[49,  'In',  'Indium',  114.818,],
[50,  'Sn',  'Tin', 118.710,],
[51,  'Sb',  'Antimony',    121.760,],
[52,  'Te',  'Tellurium',   127.60,],
[53,  'I',   'Iodine',  126.90447,],
[54,  'Xe',  'Xenon',   131.293,],
[55,  'Cs',  'Caesium', 132.90545196,],
[56,  'Ba',  'Barium',  137.327,],
[57,  'La',  'Lanthanum',   138.90547,],
[58,  'Ce',  'Cerium',  140.116,],
[59,  'Pr',  'Praseodymium',        140.90766,],
[60,  'Nd',  'Neodymium',   144.242,],
[61,  'Pm',  'Promethium',  145,],
[62,  'Sm',  'Samarium',    150.36,],
[63,  'Eu',  'Europium',    151.964,],
[64,  'Gd',  'Gadolinium',  157.25,],
[65,  'Tb',  'Terbium', 158.925354,],
[66,  'Dy',  'Dysprosium',  162.500,],
[67,  'Ho',  'Holmium', 164.930328,],
[68,  'Er',  'Erbium',  167.259,],
[69,  'Tm',  'Thulium', 168.934218,],
[70,  'Yb',  'Ytterbium',   173.045,],
[71,  'Lu',  'Lutetium',    174.9668,],
[72,  'Hf',  'Hafnium', 178.486,],
[73,  'Ta',  'Tantalum',    180.94788,],
[74,  'W',   'Tungsten',    183.84,],
[75,  'Re',  'Rhenium', 186.207,],
[76,  'Os',  'Osmium',  190.23,],
[77,  'Ir',  'Iridium', 192.217,],
[78,  'Pt',  'Platinum',    195.084,],
[79,  'Au',  'Gold',    196.966570,],
[80,  'Hg',  'Mercury', 200.592,],
[81,  'Tl',  'Thallium',    204.38,],
[82,  'Pb',  'Lead',    207.2,],
[83,  'Bi',  'Bismuth', 208.98040,],
[84,  'Po',  'Polonium',    209,],
[85,  'At',  'Astatine',    210,],
[86,  'Rn',  'Radon',   222,],
[87,  'Fr',  'Francium',    223,],
[88,  'Ra',  'Radium',  226,],
[89,  'Ac',  'Actinium',    227,],
[90,  'Th',  'Thorium', 232.0377,],
[91,  'Pa',  'Protactinium',    231.03588,],
[92,  'U',   'Uranium', 238.02891,],
[93,  'Np',  'Neptunium',   237,],
[94,  'Pu',  'Plutonium',   244,],
[95,  'Am',  'Americium',   243,],
[96,  'Cm',  'Curium',  247,],
[97,  'Bk',  'Berkelium',   247,],
[98,  'Cf',  'Californium', 251,],
[99,  'Es',  'Einsteinium', 252,],
[100, 'Fm',  'Fermium', 257,],
[101, 'Md',  'Mendelevium', 258,],
[102, 'No',  'Nobelium',    259,],
[103, 'Lr',  'Lawrencium',  262,],
[104, 'Rf',  'Rutherfordium',   267,],
[105, 'Db',  'Dubnium', 270,],
[106, 'Sg',  'Seaborgium',  269,],
[107, 'Bh',  'Bohrium', 270,],
[108, 'Hs',  'Hassium', 270,],
[109, 'Mt',  'Meitnerium',  278,],
[110, 'Ds',  'Darmstadtium',    281,],
[111, 'Rg',  'Roentgenium', 281,],
[112, 'Cn',  'Copernicium', 285,],
[113, 'Nh',  'Nihonium',    286,],
[114, 'Fl',  'Flerovium',   289,],
[115, 'Mc',  'Moscovium',   289,],
[116, 'Lv',  'Livermorium', 293,],
[117, 'Ts',  'Tennessine',  293,],
[118, 'Og',  'Oganesson',   294,],
], dtype=object)
