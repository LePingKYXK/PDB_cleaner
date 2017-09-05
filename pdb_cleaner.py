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
python PATH/pdb_cleaner.py

here, PATH is the directory contains this script, 
Then, the program will ask you to specified the directory that the PDB files located, 
and how to deal with multiple chains (keep all the chains or just one of them).
If you choose "one", the program will choose the longest chain in the PDB file.
'''

from numpy import char
import numpy as np
import pandas as pd
import os, sys, time



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
        print('==== {:} has sequence gaps! ===='.format(f))
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


def save_report(path, number, altloc_info, non_std_Res,
                negativeSeq, seqGap_info, multiChains):
    report = ''.join(("special_PDB_in_", str(number + 1), "_database.txt"))
    line = ''.join(("\n", "-" * 50, "\n"))
    string = 'The files below have'

    with open(os.path.join(path, report), 'w') as fw:
        fw.write(line)
        title1 = ' '.join((string, 'alternate location', '\n'))
        fw.write(title1)
        head = ''.join(('{:<}\t', 'alternate locations: '))
        for a in altloc_info:
            fmt1 = ''.join((head, "{:<4}" * len(a[1][1:]), "\n"))
            fw.write(fmt1.format(a[0], *a[1][1:]))
        
        fw.write(line)
        title2 = ' '.join((string, 'non-standard residues', '\n'))
        fw.write(title2)
        head = ''.join(('{:<}\t', 'non-standard residues: '))
        for n in non_std_Res:
            fmt2 = ''.join((head, "{:<4}" * len(n[1]), "\n"))
            fw.write(fmt2.format(n[0], *n[1]))

        fw.write(line)
        title3 = ' '.join((string, 'negative sequence number', '\n'))
        fw.write(title3)
        fmt3 = '{:<}\n'
        for s in negativeSeq:
            fw.write(fmt3.format(s))

        fw.write(line)
        title4 = ' '.join((string, 'sequence gaps','\n'))
        fw.write(title4)
        head = ''.join(('{:<}\t', 'sequence gaps: '))
        for g in seqGap_info:
            fmt4 = ''.join((head, "{:}" * len(g[1]), "\n"))
            fw.write(fmt4.format(g[0], *g[1]))

        fw.write(line)
        title5 = ' '.join((string, 'insertion code', '\n'))
        fw.write(title5)
        head = ''.join(('{:<}\t', 'insertion code: '))
        for i in insert_info:
            fmt5 = ''.join((head, "{:<4}" * len(i[1][1:]), "\n"))
            fw.write(fmt5.format(i[0], *i[1][1:]))
        
        fw.write(line)        
        title6 = ' '.join((string, 'multiple chains','\n'))
        fw.write(title6)
        head = ''.join(('{:<}\t', 'multiple chains: '))
        for m in multiChains:
            fmt6 = ''.join((head, "{:<4}" * len(m[1]), "\n"))
            fw.write(fmt6.format(m[0], *m[1]))


def save_cleaned_PDB(path, f, pdb_info, nonstdRes, poly):
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
    groups = pdb_info.groupby(['Seq_Num', 'ChainID'], sort=False)
    pdb_info = groups.apply(lambda x:
                            x.drop_duplicates(subset=["AtomTyp"], keep='first')
                            if len(groups['Alt_Loc']) >= 2 else x)

    #### delete the non-standard amino acid residues, DNA and RNA
    if nonstdRes:
        pdb_info = pdb_info[~pdb_info.ResName.isin(nonstdRes[1])]

    #### replaceing the 'HETATOM' by 'ATOM'    
    pdb_info['Records'].replace("HETATM", "ATOM", inplace=True)
    
    #### After deleting non-std residues, choose only one chain.
    chains = check_multiple_chains(f, pdb_info)

    #### PDB file has multiple chains and choose the longest one.
    if (chains and poly == 'one'):
        pdb_info = pdb_info[pdb_info.ChainID == pdb_info.ChainID.mode()[0]]
        outputf = ''.join((path, f[:4], "_cleaned_single_chain.pdb"))
        output_format(pdb_info, outputf)
    #### PDB file has multiple chains and choose all chains.
    elif (chains and poly == 'all'):
        outputf = ''.join((path, f[:4], "_cleaned_multichains.pdb"))
        output_format(pdb_info, outputf)
    #### PDB file has only one chain.
    else:
        outputf = ''.join((path, f[:4], "_cleaned_one_chain.pdb"))
        output_format(pdb_info, outputf)


def output_format(pdb_info, outputf):
    space = ' '
    with open(outputf, 'w') as fw:
        for line in pdb_info.values:
            fmt = ''
            if len(line[2]) <= 3:
                fmt = "{:<6}{:>5}{:<2}{:<3}{:<}{:>3}{:>2}{:>4}{:<}{:>11}{:>8}{:>8}{:>6}{:>6}{:>12}{:>2}\n"
            elif len(line[2]) == 4:
                fmt = "{:<6}{:>5}{:<1}{:<3}{:<}{:>3}{:>2}{:>4}{:<}{:>11}{:>8}{:>8}{:>6}{:>6}{:>12}{:>2}\n"
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
        str1 = "The cleaning work on {:} file has completed.\n"
        str2 = "The cleaned PDB file is saved as {:}.\n"
        fmtprint = ''.join((str1, str2))
        print(fmtprint.format(f, outputf))


#############################################################################
if __name__ == "__main__":
    pathstr = '\nPlease type the directory contains PDB files: \n'
    typestr = '\nIf you want to retain all chains, please type: all\n' \
              'If you want to keep only one chain, please type: one \n'

    path = input(pathstr)
    poly = input(typestr)
    pdbfiles = find_PDB_files(path)
    
    altloc_info = []
    non_std_Res = []
    negativeSeq = []
    seqGap_info = []
    insert_info = []
    multiChains = []
    
    timefmt = "The Used Time in this step is {:.4f} Seconds"
    initial_time = time.time()

    for i, f in enumerate(pdbfiles):
        start_time = time.time()

        line = ''.join(("\n", "-" * 50, "\n"))
        fmt0 = ''.join((line, "Check point: {:>5}\t, PDB IDS:\t {:s}"))
        print(fmt0.format(i + 1, f))

        filename = ''.join((path, f))
        pdb_info = pdb_reader(filename)

        altloc = check_altloc(f, pdb_info)    # return a tuple
        if altloc:
            print('!!!!{:} has alternative location!!!!'.format(f))
            print('The alternate locations are: {:}\n'.format(altloc[1][1:]))
            altloc_info.append(altloc)

        nonstdRes = non_std_residues(f, pdb_info)    # return a tuple
        if nonstdRes:
            print('!!!!{:} has special Residue {:}\n'.format(f, nonstdRes[1]))
            non_std_Res.append(nonstdRes)

        minusSeq = check_negative_seqnum(f, pdb_info)    # return the file name
        if minusSeq:
            print('==== {:} has negative sequence number! ====\n'.format(f))
            negativeSeq.append(minusSeq)

        gaps = check_sequence_gaps(f, pdb_info)    # return a tuple
        if gaps:
            fmt = ''.join(("sequence gaps are: ", "{:}" * len(gaps[1]), "\n"))
            print(fmt.format(*gaps[1]))
            seqGap_info.append(gaps)

        insert = check_insertion_code(f, pdb_info)    # return a tuple
        if insert:
            print('!!!!{:} has insertion code!!!!'.format(f))
            print('The insertion_code are: {:}\n'.format(insert[1][1:]))
            insert_info.append(insert)

        chains = check_multiple_chains(f, pdb_info)    # return a tuple
        if chains:
            print('=-=-= {:} has multiple chains! =-=-='.format(f))
            print('The chains are {:}\n'.format(chains[1]))
            multiChains.append(chains)
        
        save_cleaned_PDB(path, f, pdb_info, nonstdRes, poly)
        steptime = time.time() - start_time
        print(timefmt.format(steptime))

    save_report(path, i, altloc_info, non_std_Res,
                negativeSeq, seqGap_info, multiChains)
    total_time = time.time() - initial_time
    fmtend = ''.join((line, "Works Completed! Total Time: {:.4f} Seconds.\n"))
    print(fmtend.format(total_time))
