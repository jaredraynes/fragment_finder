import argparse
import numpy
import matplotlib
import lxml
import pandas
import pyteomics
import csv
import math
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




######starting

#def fragments(peptide_seq):

def mass_cal(peptide_seq):
    return(round(mass.calculate_mass(peptide_seq, average = True), 1))

def import_obs_masses(file_location):
    with open(file_location, "r") as csv_file:
        obs_masses = []
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for lines in csv_reader:
            obs_masses.append(float(lines[1]))
    return(obs_masses)

def mass_diff(prot_mass, obs_masses):
    
    mass_diffs = [prot_mass - masses for masses in obs_masses]
    return(int(min(mass_diffs) // 100))

def fragments(prot_seq, obs_masses, mass_diffs, tolerance):

    found = []
    start = 0
    s = int(min(obs_masses)//140)
    e = len(prot_seq)
    for frag in prot_seq:
        for i in range(s, e):
            for num in obs_masses:
                if math.isclose(round(mass.calculate_mass(prot_seq[start:i], average = True), 1), num, abs_tol = tolerance):
                    if prot_seq[start:i] not in found:
                        found.append(prot_seq[start:i]) 
                        found.append(round(mass.calculate_mass(prot_seq[start:i], average = True), 1))
        s += 1
        e += 1
        start += 1
    print(found)

def main():
    input = argparse.ArgumentParser()
    input.add_argument("protein_sequence", help = 'Input your protein sequence e.g. "PEPTIDE".', type = str)
    input.add_argument("obs_mass_input_file", help = 'Input the exact location of your .csv file, e.g. "C:/Users/ray07c/Documents/Parkville_data/fragment_finder/files/VCLH_T-145-DSP-04_input.csv"', type = str)
    input.add_argument("-t", "--mass_tolerance", help = 'kDa tolerance for observed masses vs calculated masses', type = float, default = 0.5)
    input.add_argument("-s", "--save_output_file", help = 'Input the exact location where you would like your file saved, e.g. "~/XXX/XXXX/masses.xlsx"', type = str)
    args = input.parse_args()
    whole_prot_mass = mass_cal(args.protein_sequence)
    observed_masses = import_obs_masses(args.obs_mass_input_file)
    mass_differences = mass_diff(whole_prot_mass, observed_masses)
    fragments(args.protein_sequence, observed_masses, mass_differences, args.mass_tolerance)

if __name__ == "__main__":
    main()