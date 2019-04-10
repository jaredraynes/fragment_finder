import argparse
import numpy
import matplotlib
import lxml
import pandas
import pyteomics
from pyteomics import mass


#A program to identify peptide fragments and their location within a protein
#using mass spectrometry data

# https://pyteomics.readthedocs.io/en/latest/index.html

# potential code to use:

#>>> from pyteomics import parser
#>>> parser.cleave('AKAKBK', parser.expasy_rules['trypsin'], 0)
#{'AK', 'BK'}

#from pyteomics import parser
#forms = parser.isoforms('PEPTIDE', variable_mods={'p': ['T'], 'ox': ['P']})
#for seq in forms: print seq

#>>> mass.calculate_mass(sequence='PEPTIDE')
#799.359964027207

#>>> from pyteomics import parser
#>>> ps = parser.parse('PEPTIDE', show_unmodified_termini=True)
#>>> mass.calculate_mass(parsed_sequence=ps)
#799.359964027207

#>>> from pyteomics import mass
#>>> mass.calculate_mass(sequence='PEPTIDE', ion_type='M', charge=2)
#400.6872584803735

#>>> mass.calculate_mass(sequence='PEP', ion_type='b', charge=1)
#324.15539725264904

#>>> mass.calculate_mass(sequence='TIDE', ion_type='y', charge=1)
#477.219119708098

#from https://pyteomics.readthedocs.io/en/latest/examples/example_msms.html

def fragments(peptide, types=('b', 'y'), maxcharge=1):
    """
    The function generates all possible m/z for fragments of types
    `types` and of charges from 1 to `maxharge`.
    """
    for i in range(1, len(peptide)-1):
        for ion_type in types:
            for charge in range(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    return mass.fast_mass(
                            peptide[:i], ion_type=ion_type, charge=charge)
                else:
                    return   mass.fast_mass(
                            peptide[i:], ion_type=ion_type, charge=charge)



######starting

#def fragments(peptide_seq):

def mass_cal(peptide_seq):
    return(round(mass.calculate_mass(peptide_seq), 2))


if __name__ == "__main__":
    main()