#!/usr/bin/env python3

""" The aim of this script is for cleanning up the PDB files.
Most of time, the PDB files are complicated, which have lots of redundant information.
The detailed descriptions about them are shown in the README file.

##################################################################################
author: Dr. Huan Wang
email: huan.wang@mail.huji.ac.il
copyright: The Hebrew University of Jerusalem, and The Open University of Israel.
##################################################################################

# How to run this script:
python pdb_cleaner.py

Then, the program will ask you to specified the directory that the PDB files located, 
and how to deal with multiple chains (keep all the chains or just one of them).
If you choose "one", the program will choose the longest chain in the PDB file.
"""

from datetime import datetime
from numpy import char
import numpy as np
import pandas as pd
import os, re, sys, time


pathstr = "\nPlease type the directory contains PDB files: \n"
keepstr = ("\nIf you want to retain all chains, please type: all\n"
           "If you want to keep only one chain, please type: one\n")
rmh_str = "\nDo you want to remove all hydrogen atoms? (y or n)\n"
rep_str = ("\nIf you only want to get a report.txt file, please type: r\n"
           "If you want to clean PDB files and get a report, please type: c\n")
           
filepath = input(pathstr)
k_option = input(keepstr).lower()
remove_H = input(rmh_str).lower()
r_option = input(rep_str).lower()


def find_PDB_files(path):
    suffix = ".pdb"
    return (f for f in os.listdir(path) if f.endswith(suffix))


def pdb_reader(filename):
    method = ""
    data = []

    items = ["Records",
             "AtomSeq",
             "AtomTyp",
             "Alt_Loc",
             "ResName",
             "ChainID",
             "Seq_Num",
             "InsCode",
             "Coord_X",
             "Coord_Y",
             "Coord_Z",
             "SD_Occp",
             "SD_Temp",
             "Element",
             "Charges"]

    with open(filename, "r") as fo:
        temp = []
        mods = []
        for line in fo:
            if line.startswith("EXPDTA"):
                method = re.search(r"EXPDTA\s+(.+\w+)", line).group(1)
                print("\nExperiment Method:\t{:}\n".format(method))
                break

        if "NMR" in method:
            pdb_df = np.array([])
            
            for line in fo:
                #### read PDB into a multi-index DataFrame
                if line.startswith("MODEL"):
                    model = line.strip()
                    mods.append(model)
                    
                elif line.startswith(("ATOM", "HETATM", "TER")):
                    info = pdb_structure(line)
                    temp.append(info)
                              
                elif line.startswith("ENDMDL"):
                    data = np.asarray(temp).reshape(-1, len(items))
                    temp = []

                    data = np.column_stack(([model] * len(data), data))
                    
                    if pdb_df.size:
                        pdb_df = np.concatenate((pdb_df, data))
                    else:
                        pdb_df = data
                        
            cols = ["Model"] + items
            pdb_df = pd.DataFrame(pdb_df, columns=cols)
            
        #### X-ray or other non-NMR experiments
        elif "NMR" not in method or method == "":
            for line in fo:
                if line.startswith(("ATOM", "HETATM", "TER")):
                    info = pdb_structure(line)
                    data.append(info)
                        
            pdb_df = pd.DataFrame(data, columns=items)
            """
            if line.startswith("REMARK 465"):
                pass
            
            elif line.startswith(("SIGATM", "ANISOU", "SIGUIJ")):
                pass
            """
    pdb_df = pd.DataFrame(pdb_df, columns=items)
    
    #### locates the last "TER" in the sequence.
    terminal_id = pdb_df[pdb_df["Records"] == "TER"].index.tolist()

    #### return the experiment method, pdb information, and ligands/solvents.
    if pdb_df.shape[0] > terminal_id[-1]:
        pdb = pdb_df.iloc[:terminal_id[-1]+1, :]
        ligand = pdb_df.iloc[terminal_id[-1]+1:, :]
        return method, pdb, ligand

    elif pdb_df.shape[0] == terminal_id[-1]:
        return method, pdb_df, pd.DataFrame()

    else: # pdb_df.shape[0] < terminal_id[-1], 
        print("Error!")


def pdb_structure(string):
    pdb = [string[0:6].strip(),      # 0.  record name
           string[6:12].strip(),     # 1.  atom serial number 
           string[12:16].strip(),    # 2.  atom name with type
           string[16],               # 3.  alternate locatin indicator
           string[17:20].strip(),    # 4.  residue name
           string[21],               # 5.  chain identifier
           string[22:26].strip(),    # 6.  residue sequence number
           string[26],               # 7.  insertion code
           string[30:38].strip(),    # 8.  coordinates X
           string[38:46].strip(),    # 9.  coordinates Y
           string[46:54].strip(),    # 10. coordinates Z
           string[54:60].strip(),    # 11. standard deviation of occupancy
           string[60:66].strip(),    # 12. standard deviation of temperature
           string[76:78].strip(),    # 13. element symbol
           string[98:80].strip()]    # 14. charge on the atom
    return pdb


AMINO_ACIDS = char.asarray(["Ala", "Arg", "Asn", "Asp",
                            "Cys", "Gln", "Glu", "Gly",
                            "His", "Ile", "Leu", "Lys",
                            "Met", "Phe", "Pro", "Ser", 
                            "Thr", "Trp", "Tyr", "Val"]).upper()


def check_ligand(filename, ligand):
    """ This function check if solvent exist in the PDB file."""
    if ligand.empty:
        print("No ligands in {:}".format(filename))
    else:
        return (filename, ligand.ResName.unique())


def check_altloc(filename, pdb_df):
    """ This function deals with the alternate locations."""
    altloc = pdb_df.Alt_Loc.unique()
    if len(altloc) > 1:
        return (filename, np.sort(altloc))


def non_std_residues(filename, pdb_df):
    """ This function checks the non-standard amino acid residues."""
    nonstdRes = pdb_df.ResName[~pdb_df.ResName.isin(AMINO_ACIDS)].unique()
    if nonstdRes.any():
        return (filename, nonstdRes)


def check_negative_seqnum(filename, pdb_df):
    """ This function reports the negative sequence numbers,
    although it may not very important."""
    seq_num = pdb_df.Seq_Num.values.astype(int)
    if (seq_num < 0).any():
        return filename


def check_sequence_gaps(filename, pdb_df):
    """ The aim of this function is for find the sequence gaps by checking
    the differences of the sequence numbers."""
    seq_num = pdb_df.Seq_Num.values.astype(int)
    seq_diff = np.abs(np.diff(seq_num))

    if np.any(seq_diff > 1):
        gap_id = np.where(seq_diff > 1)[0]
        gap_head = map(":".join,
                       pdb_df[["ChainID", "Seq_Num"]].values[gap_id])
        gap_tail = map(":".join,
                       pdb_df[["ChainID", "Seq_Num"]].values[gap_id + 1])
        gap = list(zip(gap_head, gap_tail))
        return (filename, gap)


def check_insertion_code(filename, pdb_df):
    """ This function deals with the insertion code"""
    insert = pdb_df.InsCode.unique()
    if len(insert) > 1:
        return (filename, np.sort(insert))


def check_multiple_chains(filename, pdb_df):
    """ This function checks if the multiple chains exist or not 
    in the sequence."""
    chains = pdb_df.ChainID.unique()
    if len(chains) > 1:
        return (filename, chains)


def check_hydrogen(filename, pdb_df):
    """ This function checks if the hydrogen atoms exist or not 
    in the sequence."""
    if pdb_df.AtomTyp.str.startswith("H").any():
        return filename



def save_report(path, ligand_info, altloc_info, non_std_Res, hydrogens,
                seqGap_info, insert_info, multiChains, negativeSeq, drawline):
    report = "special_PDB_files.txt"
    string = "The files below have"

    with open(os.path.join(path, report), "w") as fw:
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        head = "".join(("Summary\t(", now, ")\n"))
        fw.write(head)
        
        title_lig = " ".join((drawline, string, "ligands\n"))
        fw.write(title_lig)
        head_lig = "\t".join(("{:<}", "ligands: "))
        for l in ligand_info:
            fmt_lig = "".join((head_lig, "{:<6s}" * len(l[1]), "\n"))
            fw.write(fmt_lig.format(l[0], *l[1]))
            
        title_alt = " ".join((drawline, string, "alternate locations\n"))
        fw.write(title_alt)
        head_alt = "\t".join(("{:<}", "alternate locations: "))
        for a in altloc_info:
            fmt_alt = "".join((head_alt, "{:<4}" * len(a[1][1:]), "\n"))
            fw.write(fmt_alt.format(a[0], *a[1][1:]))
        
        title_nsr = " ".join((drawline, string, "non-standard residues\n"))
        fw.write(title_nsr)
        head_nsr = "\t".join(("{:<}", "non-standard residues:"))
        for n in non_std_Res:
            fmt_nsr = "".join((head_nsr, "{:<4}" * len(n[1]), "\n"))
            fw.write(fmt_nsr.format(n[0], *n[1]))
        
        title_h = " ".join((drawline, string, "hydrogen atoms\n"))
        fw.write(title_h)
        fmt_h = "{:<}\n"
        for h in hydrogens:
            fw.write(fmt_h.format(h))

        title_gap = " ".join((drawline, string, "sequence gaps\n"))
        fw.write(title_gap)
        head_gap = "\t".join(("{:<}", "sequence gaps:"))
        for g in seqGap_info:
            fmt_gap = "".join((head_gap, "{:}" * len(g[1]), "\n"))
            fw.write(fmt_gap.format(g[0], *g[1]))

        title_ins = " ".join((drawline, string, "insertion code\n"))
        fw.write(title_ins)
        head_ins = "\t".join(("{:<}", "insertion code:"))
        for i in insert_info:
            fmt_ins = "".join((head_ins, "{:<4}" * len(i[1][1:]), "\n"))
            fw.write(fmt_ins.format(i[0], *i[1][1:]))
        
        title_mlc = " ".join((drawline, string, "multiple chains\n"))
        fw.write(title_mlc)
        head_mlc = "\t".join(("{:<}", "multiple chains: "))
        for m in multiChains:
            fmt_mlc = "".join((head_mlc, "{:<4}" * len(m[1]), "\n"))
            fw.write(fmt_mlc.format(m[0], *m[1]))
            
        title_ngs = " ".join((drawline, string, "negative sequence number\n"))
        fw.write(title_ngs)
        fmt_ngs = "{:<}\n"
        for s in negativeSeq:
            fw.write(fmt_ngs.format(s))


def save_cleaned_PDB(outputf, pdb_df):
    space = " "
    normal = "{:<6}{:>5}{:<2}{:<3}{:<}{:>3}{:>2}{:>4}{:<}{:>11}{:>8}{:>8}{:>6}{:>6}{:>12}{:>2}\n"
    special = "{:<6}{:>5}{:<1}{:<3}{:<}{:>3}{:>2}{:>4}{:<}{:>11}{:>8}{:>8}{:>6}{:>6}{:>12}{:>2}\n"
    fmt = ""
    with open(outputf, "w") as fw:
        for line in pdb_df.values:
            if len(line[2]) <= 3:
                fmt = normal
            elif len(line[2]) == 4:
                fmt = special
            fw.write(fmt.format(line[0],   # 0. reconrd name "ATOM"
                                line[1],   # 1. atom serial number
                                space,     # 2. blank spaces
                                line[2],   # 3. atom name with type
                                line[3],   # 4. alternate locatin indicator
                                line[4],   # 5. residue name
                                line[5],   # 6. chain ID
                                line[6],   # 7. residue sequence number
                                line[7],   # 8. insertion code
                                line[8],   # 9. coordinates X
                                line[9],   #10. coordinates Y
                                line[10],  #11. coordinates Z
                                line[11],  #12. standard deviation of occupancy
                                line[12],  #13. standard deviation of temperature
                                line[13],  #14. element symbol
                                line[14])) #15. charge on the atom


def main(path, keep, hydrogen, report):
    """ Workflow:
    (1) Collect all the PDB files in the given directory;
    
    (2) In each PDB file, check the following items:
        (2.1) alternate locations;
        (2.2) non-standard amino acid residues;
        (2.3) negative sequence numbers (less important);
        (2.4) sequence gaps;
        (2.5) insertion code;
        (2.6) multiple chains;
        (2.7) hydrogen atoms;
        (2.8) ** to do: missing atoms **
        
    (3) Clean the PDB files if the aforementioned items exist,
        and user choose "c" for clean, with following options:
        
        (3.1) keep all chains if the user specified "all";
              keep the longest chain (or the 1st chain, when each
              chains have the same length), if the user specified
              "one";
               
        (3.2) remove hydrogein, if the user specified "y";
        
        (3.3) remove insert codes

        (3.4) keep the first appeared alternate location

        (3.5) remove non-standard residue(s)

        (3.6) report minus sequence number(s) and sequence gap(s)
        
    (4) Save the cleaned PDB files one by one, if user specified "c";
    
    (5) Save the summary report.
    """
    ligand_info = []
    altloc_info = []
    non_std_Res = []
    Hatoms_info = []
    negativeSeq = []
    seqGap_info = []
    insert_info = []
    multiChains = []
    
    drawline = "".join(("\n", "-" * 79, "\n"))
    time_fmt = "The Used Time in this step is {:.4f} Seconds"

    str_clean = "The cleaning work on {:} file has completed.\n"
    str_saved = "The cleaned PDB file has been saved as:\n{:}\n"
    print_fmt = "".join((str_clean, str_saved))
    
    fmt_alt = "".join(("~~~~ {:} has alternative location ~~~~\n",
                       "The alternate locations are: {:}\n"))

    fmt_nonstd = "".join(("**** {:} has non-standard Residue(s) ****\n",
                          "The non-standard Residue(s): {:}\n"))
    
    fmt_insert = "".join(("!!!! {:} has insertion code !!!!\n",
                          "The insertion codes are: {:}\n"))

    fmt_chains = "".join(("=-=-= {:} has multiple chains =-=-=\n",
                          "The chains are: {:}\n"))
    
    initial_time = time.time()

    pdbfiles = find_PDB_files(path)

    for i, f in enumerate(pdbfiles):
        start_time = time.time()
        
        fmt0 = "".join((drawline, "Check Point:{:>6},\tPDB IDS:\t{:s}"))
        print(fmt0.format(i + 1, f))

        filename = os.path.join(path, f)
        method, pdb_df, ligand = pdb_reader(filename)
        
        #### treatment for different experiment methods
        if "NMR" in method:
            pdb_df = pdb_df.set_index("Model")
            
        else: # "X-ray, or other experiment methods"
            chains = check_multiple_chains(f, pdb_df)
            if chains:
                print(fmt_chains.format(f, chains[1]))
                multiChains.append(chains)
                
            #### PDB file has multiple chains and choose the longest one.
            if (chains and keep == "one"):
                pdb_df = pdb_df[pdb_df.ChainID == pdb_df.ChainID.mode()[0]]
                fname = "".join((f[:4], "_cleaned_keep_one_chain.pdb"))
               
            #### PDB file has multiple chains and choose all chains.
            elif (chains and keep == "all"):
                fname = "".join((f[:4], "_cleaned_keep_multichains.pdb"))
                
            #### PDB file has only one chain.
            else:
                fname = "".join((f[:4], "_cleaned_original_single_chains.pdb"))

            #### parse and clean the PDB file item by item.
            lig = check_ligand(filename, ligand)    # return a tuple
            if lig:
                fmt_lig = "".join(("Oo.. {:} has ligand(s) ..oO\n",
                                   "The ligand(s) information:\t",
                                   "{:<6s}" * len(lig[1]), "\n"))
                print(fmt_lig.format(filename, *lig[1]))
                ligand_info.append(lig)

            h_atoms = check_hydrogen(filename, pdb_df)
            if h_atoms:
                print("#### {:} has hydrogen atoms ####\n".format(filename))
                Hatoms_info.append(h_atoms)
                #### remove hydrogen atoms
                if report == "c" and remove_H in ["y", "yes"]:
                    pdb_df = pdb_df[pdb_df.Element != "H"]

            insert = check_insertion_code(filename, pdb_df)    # return a tuple
            if insert:
                print(fmt_insert.format(filename, insert[1][1:]))
                insert_info.append(insert)
                #### delete the insertion residues lines
                if report == "c":
                    pdb_df = pdb_df[pdb_df.InsCode == " "]

            altloc = check_altloc(filename, pdb_df)    # return a tuple
            if altloc:
                print(fmt_alt.format(filename, altloc[1][1:]))
                altloc_info.append(altloc)
                #### delete the redundant alternate locations,
                #### only keep the first apperance
                if report == "c":
                    groups = pdb_df.groupby(["Seq_Num", "ChainID"], sort=False)
                    pdb_df = groups.apply(lambda x:
                                          x.drop_duplicates(subset=["AtomTyp"],
                                                            keep="first")
                                          if len(groups["Alt_Loc"]) >= 2 else x)
            
            nonstdRes = non_std_residues(filename, pdb_df)    # return a tuple
            if nonstdRes:
                print(fmt_nonstd.format(filename, nonstdRes[1]))
                non_std_Res.append(nonstdRes)
                #### delete the non-standard amino acid residues, DNA and RNA
                if report == "c":
                    pdb_df = pdb_df[~pdb_df.ResName.isin(nonstdRes[1])]

            minusSeq = check_negative_seqnum(filename, pdb_df)    # return the file name
            if minusSeq:
                print("---- {:} has negative sequence number ----\n".format(f))
                negativeSeq.append(minusSeq)

            gaps = check_sequence_gaps(filename, pdb_df)    # return a tuple
            if gaps:
                fmt_gap = "".join(("\___/ {:} has sequence gap(s) \___/\n",
                                   "The sequence gaps are:",
                                   " {:}" * len(gaps[1]), "\n"))
                print(fmt_gap.format(filename, *gaps[1]))
                seqGap_info.append(gaps)
                
            if report == "c":
                outputf = os.path.join(path, fname)
                save_cleaned_PDB(outputf, pdb_df)
                print(print_fmt.format(f, outputf))

        steptime = time.time() - start_time
        print(time_fmt.format(steptime))
        
    save_report(path, ligand_info, altloc_info, non_std_Res, Hatoms_info,
                 seqGap_info, insert_info, multiChains, negativeSeq, drawline)
    
    total_time = time.time() - initial_time
    end_fmt = "{:}Works Completed! Total Time: {:.4f} Seconds.\n"
    print(end_fmt.format(drawline, total_time))
    return pdb_df, ligand


################################## main #######################################
if __name__ == "__main__":
    pdb_df, ligand = main(filepath, k_option, remove_H, r_option)
