# Information about the IUPAC alphabets

protein_letters = "ACDEFGHIKLMNPQRSTVWY"
extended_protein_letters = "ACDEFGHIKLMNPQRSTVWYBXZJUO"
#   B = "Asx";  aspartic acid or asparagine (D or N)
#   X = "Xxx";  unknown or 'other' amino acid
#   Z = "Glx";  glutamic acid or glutamine (E or Q)
#   J = "Xle";  leucine or isoleucine (L or I, used in mass-spec)
#   U = "Sec";  selenocysteine
#   O = "Pyl";  pyrrolysine
ambiguous_dna_letters = "GATCRYWSMKHBVDN"
unambiguous_dna_letters = "GATC"
ambiguous_rna_letters = "GAUCRYWSMKHBVDN"
unambiguous_rna_letters = "GAUC"

#   B == 5-bromouridine
#   D == 5,6-dihydrouridine
#   S == thiouridine
#   W == wyosine
extended_dna_letters = "GATCBDSW"

# are there extended forms?
#extended_rna_letters = "GAUCBDSW"

ambiguous_dna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "GATC",
    }
ambiguous_rna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "U": "U",
    "M": "AC",
    "R": "AG",
    "W": "AU",
    "S": "CG",
    "Y": "CU",
    "K": "GU",
    "V": "ACG",
    "H": "ACU",
    "D": "AGU",
    "B": "CGU",
    "X": "GAUC",
    "N": "GAUC",
    }

ambiguous_dna_complement = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
    }

ambiguous_rna_complement = {
    "A": "U",
    "C": "G",
    "G": "C",
    "U": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
    }


def _make_ranges(dict):
    d = {}
    for key, value in dict.items():
        d[key] = (value, value)
    return d

# From bioperl's SeqStats.pm
unambiguous_dna_weights = {
    "A": 347.,
    "C": 323.,
    "G": 363.,
    "T": 322.,
    }
unambiguous_dna_weight_ranges = _make_ranges(unambiguous_dna_weights)

unambiguous_rna_weights = {
    "A": unambiguous_dna_weights["A"] + 16.,  # 16 for the oxygen
    "C": unambiguous_dna_weights["C"] + 16.,
    "G": unambiguous_dna_weights["G"] + 16.,
    "U": 340.,
}
unambiguous_rna_weight_ranges = _make_ranges(unambiguous_rna_weights)

def _make_ambiguous_ranges(dict, weight_table):
    range_d = {}
    avg_d = {}
    for letter, values in dict.items():
        #Following line is a quick hack to skip undefined weights for U and O
        if len(values)==1 and values[0] not in weight_table : continue
        weights = map(weight_table.get, values)
        range_d[letter] = (min(weights), max(weights))
        total_w = 0.0
        for w in weights:
            total_w = total_w + w
        avg_d[letter] = total_w / len(weights)
    return range_d, avg_d

ambiguous_dna_weight_ranges, avg_ambiguous_dna_weights = \
               _make_ambiguous_ranges(ambiguous_dna_values,
                                      unambiguous_dna_weights)

ambiguous_rna_weight_ranges, avg_ambiguous_rna_weights = \
               _make_ambiguous_ranges(ambiguous_rna_values,
                                      unambiguous_rna_weights)

protein_weights = {
    "A": 89.09,
    "C": 121.16,
    "D": 133.10,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.18,
    "K": 146.19,
    "L": 131.18,
    "M": 149.21,
    "N": 132.12,
    #"O": 0.0, # Needs to be recorded!
    "P": 115.13,
    "Q": 146.15,
    "R": 174.20,
    "S": 105.09,
    "T": 119.12,
    #"U": 168.05, # To be confirmed
    "V": 117.15,
    "W": 204.23,
    "Y": 181.19
    }

extended_protein_values = {
    "A": "A",
    "B": "ND",
    "C": "C",
    "D": "D",
    "E": "E",
    "F": "F",
    "G": "G",
    "H": "H",
    "I": "I",
    "J": "IL",
    "K": "K",
    "L": "L",
    "M": "M",
    "N": "N",
    "O": "O",
    "P": "P",
    "Q": "Q",
    "R": "R",
    "S": "S",
    "T": "T",
    "U": "U",
    "V": "V",
    "W": "W",
    "X": "ACDEFGHIKLMNPQRSTVWY",
    #TODO - Include U and O in the possible values of X?
    #This could alter the extended_protein_weight_ranges ...
    "Y": "Y",
    "Z": "QE",
}
    
protein_weight_ranges = _make_ranges(protein_weights)

extended_protein_weight_ranges, avg_extended_protein_weights = \
               _make_ambiguous_ranges(extended_protein_values,
                                      protein_weights)


# Joao Rodrigues @ 2010
# For Center of Mass Calculation.
# Taken from http://www.chem.qmul.ac.uk/iupac/AtWt/ & PyMol
atom_weigths = {
    "?"  :   0.0,
    'H'  :   1.00794,
    'HE' :   4.002602,
    'LI' :   6.941,
    'BE' :   9.012182,
    'B'  :  10.811,
    'C'  :  12.0107,
    'N'  :  14.0067,
    'O'  :  15.9994,
    'F'  :  18.9984032,
    'NE' :  20.1797,
    'NA' :  22.989770,
    'MG' :  24.3050,
    'AL' :  26.981538,
    'SI' :  28.0855,
    'P'  :  30.973761,
    'S'  :  32.065,
    'CL' :  35.453,
    'AR' :  39.948,
    'K'  :  39.0983,
    'CA' :  40.078,
    'SC' :  44.955910,
    'TI' :  47.867,
    'V'  :  50.9415,
    'CR' :  51.9961,
    'MN' :  54.938049,
    'FE' :  55.845,
    'CO' :  58.933200,
    'NI' :  58.6934,
    'CU' :  63.546,
    'ZN' :  65.39,
    'GA' :  69.723,
    'GE' :  72.64,
    'AS' :  74.92160,
    'SE' :  78.96,
    'BR' :  79.904,   
    'KR' :  83.80,
    'RB' :  85.4678,
    'SR' :  87.62,
    'Y'  :  88.90585,
    'ZR' :  91.224,
    'NB' :  92.90638,
    'MO' :  95.94,
    'TC' :  98,
    'RU' : 101.07,
    'RH' : 102.90550,
    'PD' : 106.42,
    'AG' : 107.8682,
    'CD' : 112.411,
    'IN' : 114.818,
    'SN' : 118.710,
    'SB' : 121.760,
    'TE' : 127.60,
    'I'  : 126.90447,
    'XE' : 131.293,
    'CS' : 132.90545,
    'BA' : 137.327,
    'LA' : 138.9055,
    'CE' : 140.116,
    'PR' : 140.90765,
    'ND' : 144.24,
    'PM' : 145,
    'SM' : 150.36,
    'EU' : 151.964,
    'GD' : 157.25,
    'TB' : 158.92534,
    'DY' : 162.50,
    'HO' : 164.93032,
    'ER' : 167.259,
    'TM' : 168.93421,
    'YB' : 173.04,
    'LU' : 174.967,
    'HF' : 178.49,
    'TA' : 180.9479,
    'W'  : 183.84,
    'RE' : 186.207,
    'OS' : 190.23,
    'IR' : 192.217,
    'PT' : 195.078,
    'AU' : 196.96655,
    'HG' : 200.59,
    'TL' : 204.3833,
    'PB' : 207.2,
    'BI' : 208.98038,
    'PO' : 208.98,
    'AT' : 209.99,
    'RN' : 222.02,
    'FR' : 223.02,
    'RA' : 226.03,
    'AC' : 227.03,
    'TH' : 232.0381,
    'PA' : 231.03588,
    'U'  : 238.02891,
    'NP' : 237.05,
    'PU' : 244.06,
    'AM' : 243.06,
    'CM' : 247.07,
    'BK' : 247.07,
    'CF' : 251.08,
    'ES' : 252.08,
    'FM' : 257.10,
    'MD' : 258.10,
    'NO' : 259.10,
    'LR' : 262.11,
    'RF' : 261.11,
    'DB' : 262.11,
    'SG' : 266.12,
    'BH' : 264.12,
    'HS' : 269.13,
    'MT' : 268.14,    
}
