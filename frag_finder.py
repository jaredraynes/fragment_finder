import argparse
import numpy
import matplotlib
import lxml
import pandas
import pyteomics
import csv
import math
import multiprocessing
from pyteomics import mass


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
    return(mass_diffs)

def fragments(prot_seq, obs_masses, tolerance):

    found = []
    start = 0
    s = int(min(obs_masses)//105)
    e = int(min(obs_masses)//90)
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

def fragments_multi(prot_seq, obs_mass, tolerance):

    found = []
    start = 0
    s = int(obs_mass)//107
    e = int(obs_mass)//95
    for frag in prot_seq:
        for i in range(s, e):
            if math.isclose(round(mass.calculate_mass(prot_seq[start:i], average = True), 1), obs_mass, abs_tol = tolerance):
                if prot_seq[start:i] not in found:
                    found.append(prot_seq[start:i])
                    found.append(round(mass.calculate_mass(prot_seq[start:i], average = True), 1))
        s += 1
        e += 1
        start += 1
    if len(found) != 0:
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
    #mass_differences = mass_diff(whole_prot_mass, observed_masses)
    multi = [(args.protein_sequence, mass, args.mass_tolerance) for mass in observed_masses]
    #fragments(args.protein_sequence, observed_masses, args.mass_tolerance)

    if __name__ == '__main__':
        with multiprocessing.Pool(processes=2) as pool:
            results = pool.starmap(fragments_multi, multi)
        print(results)

if __name__ == "__main__":
    main()
