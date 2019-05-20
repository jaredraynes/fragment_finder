import argparse
import numpy
import matplotlib
import lxml
import pandas
import pyteomics
import csv
import math
import multiprocessing
import itertools as it
import xlsxwriter
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
                    find = ['Single',
                            prot_seq[start],
                            int(start + 1),
                            prot_seq[i-1],
                            int(i),
                            obs_mass,
                            round(mass.calculate_mass(prot_seq[start:i], average = True), 1),
                            round(obs_mass - round(mass.calculate_mass(prot_seq[start:i], average = True), 1), 1)]
                    found.append(find)
                else:
                    find = ['Double',
                            prot_seq[start],
                            int(start + 1),
                            prot_seq[i-1],
                            int(i),
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
    input.add_argument("-s", "--save_output_file", help = 'Input the exact location where you would like your file saved and its name, e.g. "~/XXX/XXXX/masses"', type = str)
    input.add_argument("-c", "--number_of_cores", help = 'Input the number of processing cores your computer has. The default is 2', type = int, default = 2)

    args = input.parse_args()
    dataframe = import_dataframe(args.obs_mass_input_file)
    observed_masses = import_obs_masses(dataframe)
    #making a list of argument inputs for the multiprocessing
    multi = [(args.protein_sequence, mass, dataframe, args.mass_tolerance) for mass in observed_masses]
    #converting the protein amino acid sequence into a list of individual amino acids so that the cut
    #location can be simply inserted by using indexing
    amino_list_single_print = [a for a in args.protein_sequence]
    amino_list_single_save = [a for a in args.protein_sequence]
    amino_list_double_print = [a for a in args.protein_sequence]
    amino_list_double_save = [a for a in args.protein_sequence]
    #creating empty string to be able to convert the above lists back into single strings
    rejoined_single_print = ''
    rejoined_single_save = ''
    rejoined_double_print = ''
    rejoined_double_save = ''

    #building the multiprocessing capability of the program
    if __name__ == '__main__':
        with multiprocessing.Pool(processes = args.number_of_cores) as pool:
            results = pool.starmap(fragments_multi, multi)
            #combining the multiprocessing results to put into a pandas dataframe
            combined = [index for line in results for index in line]
            df1 = pandas.DataFrame(combined, columns = ['# Cuts', 'Nterm AA', 'Nterm Num', 'Cterm AA', 'Cterm Num', 'M(obs)', 'M(calc)', 'deltaM'])
            #Selecting the intensity 'I' parameter and 'M(obs)' so that can merge it to the dataframe
            df_i = dataframe[['M(obs)', 'I']]
            df1_i = pandas.merge(df1, df_i, on= 'M(obs)', how='right')
            df1_i.dropna(how = 'any', inplace = True)
            #converting the intensity to a percentage
            percent_i = [round(((num / max(df1_i['I'])) * 100), 2) for num in df1_i['I']]
            df1_i['I'] = percent_i
            df1_i.rename(columns={'I':'% Intensity'}, inplace=True)
            df1_i.sort_values(['# Cuts', 'M(obs)'], ascending = [False, True], inplace=True)
            #changing the data types so the amino acid number is an integer
            df1_i['Nterm Num'] = df1_i['Nterm Num'].astype(dtype='int64')
            df1_i['Cterm Num'] = df1_i['Cterm Num'].astype(dtype='int64')
            df1_i = df1_i.reset_index(drop = True)
            df1_i.index += 1
            #splitting the datafram into two dataframes, one for single cuts and one for double cuts
            mask = df1_i['# Cuts'] == 'Single'
            df2_i = df1_i[mask]
            df3_i = df1_i[~mask]

            #iterating over the protein amino acid sequence list and inserting colour formatting and the
            #index identification number from the dataframe so you know which mass has come about from
            #which cut location
            for a, b, in it.zip_longest(df2_i['Nterm Num'], df2_i.index):
                amino_list_single_print[a-2] = amino_list_single_print[a-2] + '\x1b[0;33;40m' + str(b) + '\x1b[0m'
                amino_list_single_save[a-2] = amino_list_single_save[a-2] + str(b)

            for a, b, in it.zip_longest(df3_i['Nterm Num'], df3_i.index):
                amino_list_double_print[a-2] = amino_list_double_print[a-2] + '\x1b[6;36;40m' + str(b) + '\x1b[0m'
                amino_list_double_save[a-2] = amino_list_double_save[a-2] + str(b)

            for a, b, in it.zip_longest(df3_i['Cterm Num'], df3_i.index):
                amino_list_double_print[a-2] = amino_list_double_print[a-2] + '\x1b[0;35;40m' + str(b) + '\x1b[0m'
                amino_list_double_save[a-2] = amino_list_double_save[a-2] + str(b)

        #rejoining the amino acid sequence lists back into a single string
        rejoined_single_save = rejoined_single_save.join(amino_list_single_save)
        rejoined_double_save = rejoined_double_save.join(amino_list_double_save)

        #using xlsxwriter to write the dataframes and amino acid sequences to excel spreadsheets
        writer = pandas.ExcelWriter(args.save_output_file + '.xlsx', engine='xlsxwriter')
        df2_i.to_excel(writer, sheet_name='Single_Cut')
        df3_i.to_excel(writer, sheet_name='Double_Cut')

        workbook = writer.book

        # Set up some formats to use
        red = workbook.add_format({'color': 'red', 'font_name':'Courier New', 'bold': True})
        font = workbook.add_format({'font_name':'Courier New'})
        text_wrap = workbook.add_format({'text_wrap': True})
        alignment = workbook.add_format({'align': 'center'})

        #converting the dataframe index to a list of strings so can search for them below
        index = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']

        #adding red colour and bold to the cut locations in the amino acid string sequence
        worksheet_single = writer.sheets['Single_Cut']
        format_seq_single = []
        for base in rejoined_single_save:
            if base in index:
                format_seq_single.extend((red, base))
            else:
                format_seq_single.extend((font, base))
        worksheet_single.write_rich_string(
        'A' + str(len(df2_i.index) + 2),
        *format_seq_single)
        worksheet_single.set_column('A:A', 25, text_wrap)
        worksheet_single.set_column('B:J', 10, alignment)

        #adding red colour and bold to the cut locations in the amino acid string sequence
        worksheet_double = writer.sheets['Double_Cut']
        format_seq_double = []
        print(rejoined_double_save)
        for base in rejoined_double_save:
            if base in index:
                format_seq_double.extend((red, base))
            else:
                format_seq_double.extend((font, base))
        worksheet_double.write_rich_string(
        'A' + str(len(df3_i.index) + 2),
        *format_seq_double)
        worksheet_double.set_column('A:A', 25, text_wrap)
        worksheet_double.set_column('B:J', 10, alignment)

        workbook.close()

        print(df2_i.to_string())
        print(rejoined_single_print.join(amino_list_single_print))
        print(df3_i.to_string())
        print(rejoined_double_print.join(amino_list_double_print))

if __name__ == "__main__":
    main()
