import pandas as pd
import numpy as np
from ast import literal_eval
import os
from decimal import *
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("Directory", type=str, help="Directory to where raw sumstat file is located and where formatted "
                                                "file will appear. Don't forget a slash at the end")
parser.add_argument("File_in", type=str, help="Name of raw sumstat file")
parser.add_argument("File_out", type=str, help="Name of formatted sumstat file")
parser.add_argument("--N_case", type=str, help="Set the number of cases")
parser.add_argument("--N_control", type=str, help="Set the number of controls")


class Format:

    def __init__(self, in_file, out_file, directory, N_case, N_control):

        if N_case != None and N_control != None:
            self.N_flag = True
            self.N_case = N_case
            self.N_control = N_control
        elif np.logical_xor(N_case == None, N_control == None):
            raise ValueError("You must include both N_case and N_control or neither")
        else:
            self.N_flag = False

        self.in_file = in_file
        self.out_file = out_file
        self.directory = directory
        self.ext_strip = in_file.replace(".tsv", "")
        self.clean_file = out_file.replace(".tsv", "") + "_temp.tsv"
        self.v_dir = os.getcwd() + "/Variant/"
        self.i_dir = os.getcwd() + "/Info/"
        self.t_dir = os.getcwd() + "/Temp/"
        getcontext().prec = 13

        self.flow_control()

    def EAF_converter(self, raw_EAF):
        calc_EAF = Decimal(1) - Decimal(raw_EAF)
        clean_EAF = str(calc_EAF).rstrip("0")
        return clean_EAF


    def clean_raw_file(self):

        o = open(self.directory + self.in_file, "r")
        n = open(self.directory + self.clean_file, "w+")

        n.write("variant\tCHR\tPOS\tA1\tA2\tEAF\tN\tBeta\tse\tP\tminor_allele\n")

        count = 1
        for line in o:
            if count > 1:
                strip = line.strip("\n")
                split = strip.split()
                if len(split) == 12:
                    v_data = split[0].split(":")
                    pos_key = split[0]
                    chrm = v_data[0]
                    bp = v_data[1]
                    maf = split[2]
                    N = split[5]
                    beta = split[8]
                    se = split[9]
                    p = split[11]
                    A1 = v_data[3]
                    A2 = v_data[2]
                    if v_data[3] == split[1]:
                        minor_allele = "A1"
                    elif v_data[2] == split[1]:
                        minor_allele = "A2"
                    else:
                        raise ValueError("No MA")

                    record = pos_key + "\t" + chrm + "\t" + bp + "\t" + A1 + "\t" + A2 + "\t" + maf + "\t" + N + "\t" \
                             + beta + "\t" + se + "\t" + p + "\t" + minor_allele

                    n.write(record + "\n")

                elif len(split) == 11:
                    v_data = split[0].split(":")
                    pos_key = split[0]
                    chrm = v_data[0]
                    bp = v_data[1]
                    maf = split[2]
                    N = split[4]
                    beta = split[7]
                    se = split[8]
                    p = split[10]
                    A1 = v_data[3]
                    A2 = v_data[2]
                    if v_data[3] == split[1]:
                        minor_allele = "A1"
                    elif v_data[2] == split[1]:
                        minor_allele = "A2"
                    else:
                        raise ValueError("No MA")

                    record = pos_key + "\t" + chrm + "\t" + bp + "\t" + A1 + "\t" + A2 + "\t" + maf + "\t" + N + "\t" \
                             + beta + "\t" + se + "\t" + p + "\t" + minor_allele

                    n.write(record + "\n")

                else:
                    raise IndexError("Unexpected Number of Columns")

            count += 1


    def split_file(self):

        o = open(self.directory + self.clean_file, "r")
        n = open(self.t_dir + "temp.txt", "w+")

        current = 0
        line_count = 0
        for line in o:
            if line_count > 0:
                s1 = line.split()
                s2 = s1[0].split(":")
                try:
                    chr = literal_eval(s2[0])

                except ValueError:
                    chr = 23

                if chr == current:
                    n.write(line)
                else:
                    if s2[0] != "X":
                        current = chr
                    else:
                        current = 23
                    n.close()
                    n = open(self.t_dir + self.ext_strip + "_" + str(current) + ".tsv", "w+")
                    n.write("variant\tCHR\tPOS\tA1\tA2\tEAF\tN\tBeta\tse\tP\tminor_allele\n")
                    n.write(line)
            line_count += 1


    def attach_info(self, chr_num):

        raw = pd.read_csv(self.t_dir + self.ext_strip + "_" + str(chr_num) + ".tsv", delim_whitespace=True,
                          index_col=0, dtype=str)

        if self.N_flag:
            raw.loc[:, "N_case"]  = self.N_case
            raw.loc[:, "N_control"] = self.N_control
            raw.rename(columns={"N": "N_total"}, inplace=True)

        self.with_info = pd.DataFrame()

        for chunk in pd.read_csv(self.i_dir + "Info_chr" + str(chr_num) + ".tsv", delim_whitespace=True, dtype=str,
                                  index_col=0, chunksize=300000):

            chunk = chunk[["INFO"]]
            portion = raw.join(chunk, how="inner")
            portion.dropna(subset=["Beta", "se", "P"], inplace=True)
            self.with_info = self.with_info.append(portion)

        del raw


    def attach_EAF_rsid(self, chr_num):

        variant_files = ["variants_studyEAF", "variants_popEAF"]

        if self.N_flag:
            col_order = ["SNP", "CHR", "POS", "A1", "A2", "minor_allele", "EAF", "Beta", "se", "P", "N_total", "N_case",
                         "N_control", "INFO", "EAF_source"]
        else:
            col_order = ["SNP", "CHR", "POS", "A1", "A2", "minor_allele", "EAF", "Beta", "se", "P", "N", "INFO",
                         "EAF_source"]

        self.EAF_no_calc = pd.DataFrame()
        self.EAF_calc = pd.DataFrame()
        self.EAF_pop = pd.DataFrame()

        for file in variant_files:

            for chunk in pd.read_csv(self.v_dir + file + "_" + str(chr_num) + ".tsv", delim_whitespace=True, dtype=str,
                                     index_col=0, chunksize=300000):

                chunk = chunk[["rsid", "AF"]]
                portion = self.with_info.join(chunk, how="inner")

                if file == variant_files[0]:

                    portion.drop(columns=["AF"], inplace=True)
                    portion.rename(columns={"rsid": "SNP"}, inplace=True)

                    portion_no_calc = portion.loc[portion["minor_allele"] == "A1", :].copy()
                    portion_no_calc.loc[:, "EAF_source"] = "STUDY"
                    portion_no_calc = portion_no_calc[col_order]
                    self.EAF_no_calc = self.EAF_no_calc.append(portion_no_calc)

                    portion_calc = portion.loc[portion["minor_allele"] == "A2", :].copy()
                    EAF = pd.DataFrame(portion_calc.loc[:, "EAF"])
                    EAF = EAF.applymap(lambda x: self.EAF_converter(x))
                    portion_calc.loc[:, "EAF"] = EAF.loc[:, "EAF"]
                    portion_calc.loc[:, "EAF_source"] = "CALC"
                    portion_calc = portion_calc[col_order]
                    self.EAF_calc = self.EAF_calc.append(portion_calc)

                elif file == variant_files[1]:

                    portion.drop(columns=["EAF"], inplace=True)
                    portion.rename(columns={"AF": "EAF", "rsid": "SNP"}, inplace=True)
                    portion.loc[:, "EAF_source"] = "POP"
                    portion = portion[col_order]
                    self.EAF_pop = self.EAF_pop.append(portion)

        del self.with_info


    def final_append(self):

        self.final = pd.DataFrame()
        self.final = self.final.append(self.EAF_no_calc)
        del self.EAF_no_calc
        self.final = self.final.append(self.EAF_calc)
        del self.EAF_calc
        self.final = self.final.append(self.EAF_pop)
        del self.EAF_pop
        self.final.drop(columns=["minor_allele"], inplace=True)
        self.final.sort_values("POS", inplace=True)


    def flow_control(self):

        self.clean_raw_file()
        self.split_file()

        for chr_num in range(1, 24):
            self.attach_info(chr_num)
            self.attach_EAF_rsid(chr_num)
            self.final_append()
            if chr_num == 1:
                self.final.to_csv(self.directory + self.out_file, sep="\t", index=False)
            else:
                self.final.to_csv(self.directory + self.out_file, sep="\t", index=False, header=False, mode="a")

        for file in os.listdir(self.t_dir):
            os.remove(self.t_dir + file)


if __name__ == "__main__":

    args = parser.parse_args()
    f = Format(args.Directory, args.File_in, args.File_out, args.N_case, args.N_control)