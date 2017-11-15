#!/usr/bin/env python3

''' The aim of this script is for cleanning up the PDB files.
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
'''

from numpy import char
import numpy as np
import pandas as pd
import os, sys, time


pathstr = '\nPlease type the directory contains PDB files: \n'
options = '\nIf you want to retain all chains, please type: all\n' \
              'If you want to keep only one chain, please type: one \n'

filepath = input(pathstr)
k_option = input(options)


def find_PDB_files(path):
    suffix = ".pdb"
    return (f for f in os.listdir(path) if f.endswith(suffix))


def pdb_reader(filename):
    pdb_info = [] 

    items = ['Records',
             'AtomSeq',
             'AtomTyp',
             'Alt_Loc',
             'ResName',
             'ChainID',
             'Seq_Num',
             'InsCode',
             'Coord_X',
             'Coord_Y',
             'Coord_Z',
             'SD_Occp',
             'SD_Temp',
             'Element',
             'Charges']


    with open(filename, 'r') as fo:
        for line in fo:
            if line.startswith("REMARK 465"):
                pass

            elif line.startswith("SEQRES"):
                pass

            elif line.startswith("MODEL"):
                pass
            
            elif line.startswith("ENDMDL"):
                pass
            
            elif line.startswith("ATOM") or \
                 line.startswith("HETATM") or \
                 line.startswith("TER"):

                info = pdb_structure(line)
                pdb_info.append(info)
                
            elif line.startswith("SIGATM"):
                pass
            
            elif line.startswith("ANISOU"):
                pass
            
            elif line.startswith("SIGUIJ"):
                pass
        
    pdb_info = pd.DataFrame(pdb_info, columns=items)
    pdb_info['Seq_Num'] = pdb_info['Seq_Num'].astype(int)

    #####################################################################
    #### locates the last 'TER' in the sequence.
    terminal_id = pdb_info[pdb_info["Records"] == "TER"].index.tolist()

    #####################################################################    
    #### return the pdb information without solvent.
    if pdb_info.shape[0] > terminal_id[-1]:
        return pdb_info.drop(pdb_info.index[terminal_id[-1]+1:])

    elif pdb_info.shape[0] == terminal_id[-1]:
        return pdb_info

    else: # pdb_info.shape[0] < terminal_id[-1], 
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


AMINO_ACIDS = char.asarray(['Ala', 'Arg', 'Asn', 'Asp',
                            'Cys', 'Gln', 'Glu', 'Gly',
                            'His', 'Ile', 'Leu', 'Lys',
                            'Met', 'Phe', 'Pro', 'Ser', 
                            'Thr', 'Trp', 'Tyr', 'Val']).upper()


def check_altloc(f, pdb_info):
    """ This function deals with the alternate locations."""
    altloc = pdb_info.Alt_Loc.unique()
    if len(altloc) > 1:
        return (f, np.sort(altloc))


def non_std_residues(f, pdb_info):
    """ This function checks the non-standard amino acid residues."""
    nonstdRes = pdb_info.ResName[~pdb_info.ResName.isin(AMINO_ACIDS)].unique()
    if nonstdRes.any():
        return (f, nonstdRes)


def check_negative_seqnum(f, pdb_info):
    """ This function reports the negative sequence numbers, although it may not very
    important."""
    seq_num = pdb_info.Seq_Num
    if (seq_num < 0).any():
        return f


def check_sequence_gaps(f, pdb_info):
    """ The aim of this function is for find the sequence gaps by checking
    the differences of the sequence numbers."""
    seq_num = pdb_info.Seq_Num
    seq_diff = np.abs(np.diff(seq_num))

    if np.any(seq_diff > 1):
        gap_id = np.where(seq_diff > 1)[0]
        gap_head = seq_num[gap_id]
        gap_tail = seq_num[gap_id + 1]
        gap = list(zip(gap_head, gap_tail))
        return (f, gap)


def check_insertion_code(f, pdb_info):
    """ This function deals with the insertion code"""
    insert = pdb_info.InsCode.unique()
    if len(insert) > 1:
        return (f, np.sort(insert))


def check_multiple_chains(f, pdb_info):
    """ This function checks the multiple chains exist or not 
    in the sequence."""
    chains = pdb_info.ChainID.unique()
    if len(chains) > 1:
        return (f, chains)


def save_report(path, number, altloc_info, non_std_Res, negativeSeq,
                seqGap_info, insert_info, multiChains, drawline):
    report = ''.join(("special_PDB_in_", str(number), "_PDB_files.txt"))
    string = 'The files below have'

    with open(os.path.join(path, report), 'w') as fw:
        fw.write(''.join(("Summary", drawline)))
        title1 = ' '.join((string, 'alternate locations', '\n'))
        fw.write(title1)
        head = ''.join(('{:<}\t', 'alternate locations: '))
        for a in altloc_info:
            fmt1 = ''.join((head, "{:<4}" * len(a[1][1:]), "\n"))
            fw.write(fmt1.format(a[0], *a[1][1:]))
        
        fw.write(drawline)
        title2 = ' '.join((string, 'non-standard residues', '\n'))
        fw.write(title2)
        head = ''.join(('{:<}\t', 'non-standard residues: '))
        for n in non_std_Res:
            fmt2 = ''.join((head, "{:<4}" * len(n[1]), "\n"))
            fw.write(fmt2.format(n[0], *n[1]))

        fw.write(drawline)
        title3 = ' '.join((string, 'negative sequence number', '\n'))
        fw.write(title3)
        fmt3 = '{:<}\n'
        for s in negativeSeq:
            fw.write(fmt3.format(s))

        fw.write(drawline)
        title4 = ' '.join((string, 'sequence gaps','\n'))
        fw.write(title4)
        head = ''.join(('{:<}\t', 'sequence gaps: '))
        for g in seqGap_info:
            fmt4 = ''.join((head, "{:}" * len(g[1]), "\n"))
            fw.write(fmt4.format(g[0], *g[1]))

        fw.write(drawline)
        title5 = ' '.join((string, 'insertion code', '\n'))
        fw.write(title5)
        head = ''.join(('{:<}\t', 'insertion code: '))
        for i in insert_info:
            fmt5 = ''.join((head, "{:<4}" * len(i[1][1:]), "\n"))
            fw.write(fmt5.format(i[0], *i[1][1:]))
        
        fw.write(drawline)        
        title6 = ' '.join((string, 'multiple chains','\n'))
        fw.write(title6)
        head = ''.join(('{:<}\t', 'multiple chains: '))
        for m in multiChains:
            fmt6 = ''.join((head, "{:<4}" * len(m[1]), "\n"))
            fw.write(fmt6.format(m[0], *m[1]))


def save_cleaned_PDB(path, f, pdb_info, altloc, nonstdRes, keep, printfmt):
    ''' This function first clean up the PDB file, only remain one chain,
    the non-labeled and A-alternate location, throw away the non-standard
    amino acid residues, delete the insertion code lines and also change 
    the HETATM into ATOM.
        Finally, save the cleaned result with the name of XXXX_cleaned.pdb, 
    here XXXX is the corresponding original PDB file name/code.
    '''
    #### delete the insertion residues lines
    pdb_info = pdb_info[pdb_info.InsCode == ' ']
    
    #### delete the redundant alternate locations, only keep the first apperance
    if altloc:
        groups = pdb_info.groupby(['Seq_Num', 'ChainID'], sort=False)
        pdb_info = groups.apply(lambda x:
                                x.drop_duplicates(subset=["AtomTyp"],
                                                  keep='first')
                                if len(groups['Alt_Loc']) >= 2 else x)

    #### delete the non-standard amino acid residues, DNA and RNA
    if nonstdRes:
        pdb_info = pdb_info[~pdb_info.ResName.isin(nonstdRes[1])]

    #### replaceing the 'HETATOM' by 'ATOM'    
    pdb_info['Records'].replace("HETATM", "ATOM", inplace=True)
    
    #### After deleting non-std residues, choose only one chain.
    chains = check_multiple_chains(f, pdb_info)

    #### PDB file has multiple chains and choose the longest one.
    if (chains and keep == 'one'):
        pdb_info = pdb_info[pdb_info.ChainID == pdb_info.ChainID.mode()[0]]
        fname = ''.join((f[:4], "_cleaned_keep_single_chain.pdb"))
        outputf = os.path.join(path, fname)
        output_format(pdb_info, outputf)
    #### PDB file has multiple chains and choose all chains.
    elif (chains and keep == 'all'):
        fname = ''.join((f[:4], "_cleaned_keep_multichains.pdb"))
        outputf = os.path.join(path, fname)
        output_format(pdb_info, outputf)
    #### PDB file has only one chain.
    else:
        fname = ''.join((f[:4], "_cleaned_one_chain.pdb"))
        outputf = os.path.join(path, fname)
        output_format(pdb_info, outputf)


def output_format(pdb_info, outputf):
    space = ' '
    normal = "{:<6}{:>5}{:<2}{:<3}{:<}{:>3}{:>2}{:>4}{:<}{:>11}{:>8}{:>8}{:>6}{:>6}{:>12}{:>2}\n"
    special = "{:<6}{:>5}{:<1}{:<3}{:<}{:>3}{:>2}{:>4}{:<}{:>11}{:>8}{:>8}{:>6}{:>6}{:>12}{:>2}\n"
    
    with open(outputf, 'w') as fw:
        for line in pdb_info.values:
            fmt = ''
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


def main(path, keep):
    ''' Workflow:
    (1) Collect all the PDB files in the given directory;
    
    (2) In each PDB file, check the following items:
        (2.1) alternate locations;
        (2.2) non-standard amino acid residues;
        (2.3) negative sequence numbers (less important);
        (2.4) sequence gaps;
        (2.5) insertion code;
        (2.6) multiple chains;
        (2.7) *** to do: missing atoms ***
        
    (3) Clean the PDB files if the aforementioned items exist,
        with following options if protein has multiple chains;
        (3.1) keep all chains;
        (3.2) keep the longest chain (or the 1st chain, if all
        chains have the same length).
        
    (4) Save the cleaned PDB files one by one;
    
    (5) Save the summary report.
    '''
    altloc_info = []
    non_std_Res = []
    negativeSeq = []
    seqGap_info = []
    insert_info = []
    multiChains = []
    
    drawline = ''.join(("\n", "-" * 50, "\n"))
    time_fmt = "The Used Time in this step is {:.4f} Seconds"

    str_clean = "The cleaning work on {:} file has completed.\n"
    str_saved = "The cleaned PDB file has been saved as:\n{:}\n"
    print_fmt = ''.join((str_clean, str_saved))
    
    count = 0
    initial_time = time.time()

    pdbfiles = find_PDB_files(path)

    for i, f in enumerate(pdbfiles):
        start_time = time.time()
        
        fmt0 = ''.join((drawline, "Check point: {:>5}\t, PDB IDS:\t {:s}"))
        print(fmt0.format(i + 1, f))

        filename = os.path.join(path, f)
        pdb_info = pdb_reader(filename)

        altloc = check_altloc(f, pdb_info)    # return a tuple
        if altloc:
            fmt_alt = ''.join(('~~~~ {:} has alternative location ~~~~\n',
                               'The alternate locations are: {:}\n'))
            print(fmt_alt.format(f, altloc[1][1:]))
            #print('The alternate locations are: {:}\n'.format(altloc[1][1:]))
            altloc_info.append(altloc)

        nonstdRes = non_std_residues(f, pdb_info)    # return a tuple
        if nonstdRes:
            fmt_nonstd = ''.join(('**** {:} has special Residue ****',
                                  'They are: {:}\n'))
            print(fmt_nonstd.format(f, nonstdRes[1]))
            non_std_Res.append(nonstdRes)

        minusSeq = check_negative_seqnum(f, pdb_info)    # return the file name
        if minusSeq:
            print('---- {:} has negative sequence number ----\n'.format(f))
            negativeSeq.append(minusSeq)

        gaps = check_sequence_gaps(f, pdb_info)    # return a tuple
        if gaps:
            fmt_gap = ''.join(("\___/ {:} has sequence gap(s) \___/\n",
                               "The sequence gaps are: ",
                               " {:}" * len(gaps[1]), "\n"))
            print(fmt_gap.format(f, *gaps[1]))
            seqGap_info.append(gaps)

        insert = check_insertion_code(f, pdb_info)    # return a tuple
        if insert:
            fmt_insert = ''.join(("!!!! {:} has insertion code !!!!\n",
                                  "The insertion_code are: {:}\n"))
            #print('!!!!{:} has insertion code!!!!'.format(f))
            print(fmt_insert.format(f, insert[1][1:]))
            insert_info.append(insert)

        chains = check_multiple_chains(f, pdb_info)    # return a tuple
        if chains:
            fmt_chains = ''.join(("=-=-= {:} has multiple chains =-=-=\n",
                                  "The chains are: {:}\n"))
            #print('=-=-= {:} has multiple chains! =-=-='.format(f))
            print(fmt_chains.format(f, chains[1]))
            multiChains.append(chains)
        
        save_cleaned_PDB(path, f, pdb_info, altloc, nonstdRes, keep, print_fmt)
        steptime = time.time() - start_time
        print(time_fmt.format(steptime))
        count += i

    save_report(path, count, altloc_info, non_std_Res, negativeSeq,
                seqGap_info, insert_info, multiChains, drawline)
    
    total_time = time.time() - initial_time
    end_fmt = "{:}Works Completed! Total Time: {:.4f} Seconds.\n"
    print(end_fmt.format(drawline, total_time))


################################## main #######################################
if __name__ == "__main__":
    main(filepath, k_option)
