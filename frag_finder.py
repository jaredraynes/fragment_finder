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


def import_dataframe(file_location):

    if 'csv' in file_location:
        dataframe = pandas.read_csv(file_location)
        dataframe.rename(columns={'m/z':'M(obs)'}, inplace=True)
    else:
        dataframe = pandas.read_excel(file_location)
        dataframe.rename(columns={'m/z':'M(obs)'}, inplace=True)
    return(dataframe)

def mass_cal(peptide_seq):

    return(round(mass.calculate_mass(peptide_seq, average = True), 1))

def import_obs_masses(dataframe):

    return(list(dataframe['M(obs)']))

def fragments(prot_seq, obs_masses, dataframe, tolerance):

    single_cut = []
    double_cut = []
    start = 0
    s = int(min(obs_masses)//105)
    e = int(max(obs_masses)//90)
    for frag in prot_seq:
        for i in range(s, e):
            if i > len(prot_seq):
                break
            for num in obs_masses:
                if math.isclose(round(mass.calculate_mass(prot_seq[start:i], average = True), 1), num, abs_tol = tolerance):
                    if i == len(prot_seq):
                        find = [prot_seq[start] + str(start + 1),
                                str(i),
                                num,
                                round(mass.calculate_mass(prot_seq[start:i], average = True), 1),
                                round(num - round(mass.calculate_mass(prot_seq[start:i], average = True), 1), 1)]
                        single_cut.append(find)
                    else:
                        find = [prot_seq[start] + str(start + 1),
                                str(i),
                                num,
                                round(mass.calculate_mass(prot_seq[start:i], average = True), 1),
                                round(num - round(mass.calculate_mass(prot_seq[start:i], average = True), 1), 1)]
                        double_cut.append(find)
        s += 1
        e += 1
        start += 1

    df1 = pandas.DataFrame(single_cut, columns = ['Cutsite (Nterm)', 'Cterm', 'M(obs)', 'M(calc)', 'deltaM'])
    df1.sort_values('M(obs)', inplace=True)
    df2 = pandas.DataFrame(double_cut, columns = ['Cutsite (Nterm)', 'Cutsite (Cterm)', 'M(obs)', 'M(calc)', 'deltaM'])
    df2.sort_values('M(obs)', inplace=True)
    df_i = dataframe[['M(obs)', 'I']]
    df1_i = pandas.merge(df1, df_i, on= 'M(obs)', how='right')
    df1_i.dropna(how = 'any', inplace = True)
    percent_i = [round(((num / max(df1_i['I'])) * 100), 2) for num in df1_i['I']]
    df1_i['I'] = percent_i
    df1_i.rename(columns={'I':'% Intensity'}, inplace=True)

    print(df1_i.to_string(index=False))
    print(df2.to_string(index=False))

def fragments_multi(prot_seq, obs_mass, dataframe, tolerance):

    found = []
    start = 0
    s = int(obs_mass)//107
    e = int(obs_mass)//95
    for frag in prot_seq:
        for i in range(s, e):
            if i > len(prot_seq):
                break
            if math.isclose(round(mass.calculate_mass(prot_seq[start:i], average = True), 1), obs_mass, abs_tol = tolerance):
                if i == len(prot_seq):
                    find = ['single', 
                            prot_seq[start] + str(start + 1),
                            str(i),
                            obs_mass,
                            round(mass.calculate_mass(prot_seq[start:i], average = True), 1),
                            round(obs_mass - round(mass.calculate_mass(prot_seq[start:i], average = True), 1), 1)]
                    found.append(find)
                else:
                    find = ['double', 
                            prot_seq[start] + str(start + 1),
                            str(i),
                            obs_mass,
                            round(mass.calculate_mass(prot_seq[start:i], average = True), 1),
                            round(obs_mass - round(mass.calculate_mass(prot_seq[start:i], average = True), 1), 1)]
                    found.append(find)
        s += 1
        e += 1
        start += 1

    return(found)

def main():
    input = argparse.ArgumentParser()
    input.add_argument("protein_sequence", help = 'Input your protein sequence e.g. "PEPTIDE".', type = str)
    input.add_argument("obs_mass_input_file", help = 'Input the exact location of your .csv file, e.g. "C:/Users/ray07c/Documents/Parkville_data/fragment_finder/files/VCLH_T-145-DSP-04_input.csv"', type = str)
    input.add_argument("-t", "--mass_tolerance", help = 'kDa tolerance for observed masses vs calculated masses', type = float, default = 0.5)
    input.add_argument("-s", "--save_output_file", help = 'Input the exact location where you would like your file saved, e.g. "~/XXX/XXXX/masses.xlsx"', type = str)
    input.add_argument("-c", "--number_of_cores", help = 'Input the number of processing cores your computer has. The default is 2', type = int, default = 2)
    
    args = input.parse_args()
    dataframe = import_dataframe(args.obs_mass_input_file)
    observed_masses = import_obs_masses(dataframe)
    multi = [(args.protein_sequence, mass, dataframe, args.mass_tolerance) for mass in observed_masses]

    if __name__ == '__main__':
        with multiprocessing.Pool(processes = args.number_of_cores) as pool:
            results = pool.starmap(fragments_multi, multi)
            combined = [index for line in results for index in line]
            df1 = pandas.DataFrame(combined , columns = ['# Cuts', 'Cutsite (Nterm)', 'Cterm', 'M(obs)', 'M(calc)', 'deltaM'])
            df_i = dataframe[['M(obs)', 'I']]
            df1_i = pandas.merge(df1, df_i, on= 'M(obs)', how='right')
            df1_i.dropna(how = 'any', inplace = True)
            percent_i = [round(((num / max(df1_i['I'])) * 100), 2) for num in df1_i['I']]
            df1_i['I'] = percent_i
            df1_i.rename(columns={'I':'% Intensity'}, inplace=True)
            df1_i.sort_values(['# Cuts', 'M(obs)'], ascending = [False, True], inplace=True)

        print(df1_i.to_string(index=False))

if __name__ == "__main__":
    main()
