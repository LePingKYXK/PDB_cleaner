# PDB_cleaner
A Python3 script to clean up the PDB file

Most of time, the PDB files are complicated, which have lots of redundant information as shown below.

* ANISOU (data copied from 1lk2.pdb) 
<pre>
ATOM      1  N   GLY A   1      66.440  45.780   5.177  1.00 14.10           N  
<b>ANISOU</b>    1  N   GLY A   1     1908   1789   1659     99    -37     -3       N  
ATOM      2  CA  GLY A   1      65.947  45.284   3.863  1.00 12.08           C  
<b>ANISOU</b>    2  CA  GLY A   1     1484   1486   1620     47    -39     75       C  
ATOM      3  C   GLY A   1      64.961  46.275   3.303  1.00 10.99           C  
<b>ANISOU</b>    3  C   GLY A   1     1471   1204   1500     36     50    108       C  
ATOM      4  O   GLY A   1      64.683  47.291   3.943  1.00 11.91           O  
<b>ANISOU</b>    4  O   GLY A   1     1390   1542   1593    -61     88    -19       O  
</pre>

The simplest way is to delete the ANISOU lines.

* HETATM

* non-standard amino acid residues (data copied from 2o2x.pdb)
<pre>
ATOM    821  OD1 ASP A 112      25.580  11.019  35.906  1.00 12.28           O  
ATOM    822  OD2 ASP A 112      24.586   9.016  35.848  1.00 11.81           O  
<b>HETATM</b>  823  N   <b>MSE</b> A 113      25.018  10.050  30.641  1.00  9.26           N  
<b>HETATM</b>  824  CA  <b>MSE</b> A 113      25.494  10.026  29.262  1.00  9.59           C  
<b>HETATM</b>  825  C   <b>MSE</b> A 113      24.291   9.758  28.359  1.00  8.63           C  
<b>HETATM</b>  826  O   <b>MSE</b> A 113      23.362   9.026  28.750  1.00  9.51           O  
<b>HETATM</b>  827  CB  <b>MSE</b> A 113      26.563   8.959  29.078  1.00  8.81           C  
<b>HETATM</b>  828  CG  <b>MSE</b> A 113      27.157   8.896  27.700  1.00  8.23           C  
<b>HETATM</b>  829  SE  <b>MSE</b> A 113      28.681   7.732  27.499  0.75 12.65          SE  
<b>HETATM</b>  830  CE  <b>MSE</b> A 113      30.013   8.895  28.258  1.00 18.85           C  
ATOM    831  N   VAL A 114      24.306  10.362  27.178  1.00  7.56           N  
ATOM    832  CA  VAL A 114      23.308  10.072  26.129  1.00  7.93           C  
</pre>

Since PDBSlicer could not deal with non-standard residues, the simplest way is to delete them.

* Missing Residue(s) or so-called Sequence Gap(s) (data copied from 1nzj.pdb)
<pre>
ATOM   1756  O   LEU A 222      48.274   3.534  34.949  1.00 27.98           O  
ANISOU 1756  O   LEU A 222     3531   3513   3584     47    -26      6       O  
ATOM   1757  CB  LEU A 222      45.906   2.133  33.476  1.00 26.02           C  
ANISOU 1757  CB  LEU A 222     3274   3295   3315     29     58     -5       C  
ATOM   1758  N   ASN A 223      47.050   5.216  34.027  1.00 28.96           N  
ANISOU 1758  N   ASN A 223     3698   3595   3707     53      7     21       N  
ATOM   1759  CA  ASN A 223      47.326   6.262  35.028  1.00 29.62           C  
ANISOU 1759  CA  ASN A 223     3782   3735   3737     11     12    -25       C  
ATOM   1760  C   ASN A 223      48.230   5.851  36.192  1.00 30.19           C  
ANISOU 1760  C   ASN A 223     3871   3833   3767     45     -3     -6       C  
ATOM   1761  O   ASN A 223      47.951   6.165  37.354  1.00 31.23           O  
ANISOU 1761  O   ASN A 223     4074   3965   3824     76     26    -70       O  
<b>ATOM   1762  CB  ASN A 223      46.003   6.831  35.561  1.00 29.98           C</b>  
<b>ANISOU 1762  CB  ASN A 223     3798   3785   3804     36     15      8       C</b>  
<b>ATOM   1763  N   ALA A 237      50.141  13.856  28.172  1.00 30.51           N</b>  
<b>ANISOU 1763  N   ALA A 237     3895   3875   3821     32     28    -17       N</b>  
ATOM   1764  CA  ALA A 237      50.857  13.904  26.900  1.00 30.22           C  
ANISOU 1764  CA  ALA A 237     3816   3837   3827      7      2      7       C  
ATOM   1765  C   ALA A 237      52.347  13.656  27.124  1.00 30.06           C  
ANISOU 1765  C   ALA A 237     3809   3808   3803     26     11    -16       C  
ATOM   1766  O   ALA A 237      52.869  13.962  28.189  1.00 30.54           O  
ANISOU 1766  O   ALA A 237     3901   3866   3834     34    -52     -8       O  
ATOM   1767  CB  ALA A 237      50.648  15.254  26.254  1.00 30.32           C  
ANISOU 1767  CB  ALA A 237     3814   3832   3871     20      0      0       C  
ATOM   1768  N   LEU A 238      53.035  13.117  26.121  1.00 29.79           N  
ANISOU 1768  N   LEU A 238     3760   3773   3785      9     -4    -12       N  
ATOM   1769  CA  LEU A 238      54.470  12.845  26.250  1.00 29.52           C  
ANISOU 1769  CA  LEU A 238     3743   3740   3733      9      9    -18       C  
</pre>

The bold font lines indicate the **discontinuous** sequence numbers (223 ...empty... 237) due to the missing residues. We called this case as sequence gap. It is a very serious problem because **the Ramachandran subunit is defined by three adjacent residues.** It is immpossible to directly choose residue number series (222, 223, 237) and (2233, 237, 238) as the members of the Ramachandran subunit. The solution is that treat the peptide as segments, e.g. from beginning to residue number 223, then from residue number 237 to the end. If the PDB file has more than one gap, we divide it into several segments based on the locations of the gaps. Note: the discontinuous sequence number between different chains also treated as 'gap', just because it is easy for programming.

**Improvement (Nov. 16, 2017)** In the printing and report format, the chain ID was added aside to the sequence number, e.g. ('A:223', 'A:237'). Previously, only the sequence numbers between gap(s) were showed. 


* alternate locations (data copied from 3ife.pdb)
<pre>
ATOM     21  CE1 PHE A  -4      40.991  47.856  19.364  1.00 27.65           C  
ATOM     22  CE2 PHE A  -4      41.948  49.936  20.068  1.00 28.56           C  
ATOM     23  CZ  PHE A  -4      40.841  49.190  19.686  1.00 28.08           C  
ATOM     24  N  <b>A</b>GLN A  -3      46.967  45.549  21.004  0.50 23.13           N  
ATOM     25  N  <b>B</b>GLN A  -3      46.998  45.555  20.982  0.50 22.90           N  
ATOM     26  CA <b>A</b>GLN A  -3      48.373  45.164  21.046  0.50 23.16           C  
ATOM     27  CA <b>B</b>GLN A  -3      48.400  45.139  20.949  0.50 22.72           C  
ATOM     28  C  <b>A</b>GLN A  -3      48.567  43.661  20.812  0.50 22.63           C  
ATOM     29  C  <b>B</b>GLN A  -3      48.554  43.631  20.764  0.50 22.37           C  
ATOM     30  O  <b>A</b>GLN A  -3      49.384  43.259  19.986  0.50 21.47           O  
ATOM     31  O  <b>B</b>GLN A  -3      49.344  43.191  19.930  0.50 21.20           O  
ATOM     32  CB <b>A</b>GLN A  -3      49.002  45.601  22.384  0.50 23.61           C  
ATOM     33  CB <b>B</b>GLN A  -3      49.160  45.591  22.201  0.50 23.04           C  
ATOM     34  CG <b>A</b>GLN A  -3      48.488  44.854  23.614  0.50 25.25           C  
ATOM     35  CG <b>B</b>GLN A  -3      50.631  45.172  22.185  0.50 23.58           C  
ATOM     36  CD <b>A</b>GLN A  -3      46.975  44.863  23.719  0.50 25.86           C  
ATOM     37  CD <b>B</b>GLN A  -3      51.390  45.619  23.424  0.50 24.16           C  
ATOM     38  OE1<b>A</b>GLN A  -3      46.364  43.846  24.056  0.50 23.20           O  
ATOM     39  OE1<b>B</b>GLN A  -3      50.935  46.485  24.167  0.50 26.65           O  
ATOM     40  NE2<b>A</b>GLN A  -3      46.361  45.990  23.375  0.50 25.35           N  
ATOM     41  NE2<b>B</b>GLN A  -3      52.563  45.035  23.640  0.50 27.37           N  
ATOM     42  N   SER A  -2      47.792  42.842  21.521  1.00 21.72           N  
ATOM     43  CA  SER A  -2      47.888  41.386  21.401  1.00 22.23           C  
ATOM     44  C   SER A  -2      47.402  40.921  20.036  1.00 19.65           C  
ATOM     45  O   SER A  -2      48.008  40.034  19.456  1.00 20.72           O  
</pre>

  * special cases in alternate locations (data copied from 5DXX.pdb)
  <pre>
ATOM    448  N   MET A  61      48.127   9.414  21.012  1.00  8.02           N  
ANISOU  448  N   MET A  61      952    878   1219    -50    501     95       N  
ATOM    449  CA <b>A</b>MET A  61      47.494   8.918  22.231  0.58  8.39           C  
ANISOU  449  CA <b>A</b>MET A  61     1091    827   1271     24    428    219       C  
ATOM    450  CA <b>B</b>MET A  61      47.420   8.922  22.202  0.42  8.88           C  
ANISOU  450  CA <b>B</b>MET A  61     1144    895   1334    -61    457    185       C  
ATOM    451  C   MET A  61      47.346   7.404  22.267  1.00  8.78           C  
ANISOU  451  C   MET A  61     1223    782   1330     59    378    169       C  
ATOM    452  O   MET A  61      46.991   6.766  21.272  1.00 10.08           O  
ANISOU  452  O   MET A  61     1398    943   1491     14     75    122       O  
ATOM    453  CB <b>A</b>MET A  61      46.118   9.546  22.410  0.58  8.06           C  
ANISOU  453  CB <b>A</b>MET A  61      903    838   1320    372    380    212       C  
ATOM    454  CB <b>B</b>MET A  61      45.980   9.455  22.241  0.42  8.97           C  
ANISOU  454  CB <b>B</b>MET A  61      930    991   1486      2    458    168       C  
ATOM    455  CG <b>A</b>MET A  61      46.138  11.063  22.501  0.58  8.72           C  
ANISOU  455  CG <b>A</b>MET A  61     1253    809   1251    307    330    165       C  
ATOM    456  CG <b>B</b>MET A  61      45.805  10.973  22.171  0.42  9.66           C  
ANISOU  456  CG <b>B</b>MET A  61     1110   1045   1516     57    274    111       C  
ATOM    457  SD <b>A</b>MET A  61      44.516  11.746  22.852  0.58  9.87           S  
ANISOU  457  SD <b>A</b>MET A  61     1393   1136   1221    329    227    121       S  
ATOM    458  SD <b>B</b>MET A  61      44.071  11.452  21.925  0.42 11.25           S  
ANISOU  458  SD <b>B</b>MET A  61     1357   1344   1573    -97    206    -11       S  
ATOM    459  CE <b>A</b>MET A  61      43.632  11.262  21.374  0.58  8.79           C  
ANISOU  459  CE <b>A</b>MET A  61      912   1125   1304    448    193     62       C  
ATOM    460  CE <b>B</b>MET A  61      43.308  10.818  23.419  0.42 10.70           C  
ANISOU  460  CE <b>B</b>MET A  61     1120   1371   1573     26    296    -15       C  
...
ATOM   2041  N   ARG A 268      68.983  -6.030  20.233  1.00 12.62           N  
ANISOU 2041  N   ARG A 268     1676    819   2299    101    -35    523       N  
ATOM   2042  CA <b>B</b>ARG A 268      68.988  -4.603  20.530  0.60 12.88           C  
ANISOU 2042  CA <b>B</b>ARG A 268     1398    984   2513    141    107    402       C  
ATOM   2043  CA <b>C</b>ARG A 268      68.989  -4.603  20.527  0.40 12.83           C  
ANISOU 2043  CA <b>C</b>ARG A 268     1483    920   2471     82    157    473       C  
ATOM   2044  C   ARG A 268      67.641  -3.953  20.247  1.00 11.56           C  
ANISOU 2044  C   ARG A 268     1170    935   2286    -23     70    342       C  
ATOM   2045  O   ARG A 268      66.930  -4.345  19.316  1.00 12.73           O  
ANISOU 2045  O   ARG A 268     1496   1160   2181    -23     37    354       O  
ATOM   2046  CB <b>B</b>ARG A 268      70.061  -3.890  19.701  0.60 15.09           C  
ANISOU 2046  CB <b>B</b>ARG A 268     1382   1451   2901    308    108    405       C  
ATOM   2047  CB <b>C</b>ARG A 268      70.065  -3.894  19.699  0.40 14.76           C  
ANISOU 2047  CB <b>C</b>ARG A 268     1631   1189   2787    164    317    597       C  
ATOM   2048  CG <b>B</b>ARG A 268      71.428  -4.538  19.755  0.60 20.70           C  
ANISOU 2048  CG <b>B</b>ARG A 268     2380   2183   3300    506    138    227       C  
ATOM   2049  CG <b>C</b>ARG A 268      71.458  -4.466  19.860  0.40 18.82           C  
ANISOU 2049  CG <b>C</b>ARG A 268     2367   1677   3108    350    438    578       C  
ATOM   2050  CD <b>B</b>ARG A 268      72.280  -3.968  20.869  0.60 24.96           C  
ANISOU 2050  CD <b>B</b>ARG A 268     3408   2535   3540    600    301    -84       C  
ATOM   2051  CD <b>C</b>ARG A 268      72.378  -3.477  20.542  0.40 22.04           C  
ANISOU 2051  CD <b>C</b>ARG A 268     3150   1893   3329    467    727    531       C  
ATOM   2052  NE <b>B</b>ARG A 268      73.616  -4.559  20.871  0.60 27.23           N  
ANISOU 2052  NE <b>B</b>ARG A 268     3846   2843   3658    816    402   -233       N  
ATOM   2053  NE <b>C</b>ARG A 268      73.461  -3.031  19.670  0.40 25.28           N  
ANISOU 2053  NE <b>C</b>ARG A 268     3900   2201   3505    545    882    478       N  
ATOM   2054  CZ <b>B</b>ARG A 268      74.606  -4.169  20.074  0.60 29.73           C  
ANISOU 2054  CZ <b>B</b>ARG A 268     4396   3111   3790   1084    535   -418       C  
ATOM   2055  CZ <b>C</b>ARG A 268      74.657  -3.607  19.612  0.40 28.07           C  
ANISOU 2055  CZ <b>C</b>ARG A 268     4528   2513   3625    423    993    395       C  
ATOM   2056  NH1<b>B</b>ARG A 268      74.412  -3.180  19.206  0.60 30.54           N  
ANISOU 2056  NH1<b>B</b>ARG A 268     4601   3217   3787   1268    624   -504       N  
ATOM   2057  NH1<b>C</b>ARG A 268      74.925  -4.665  20.369  0.40 29.81           N  
ANISOU 2057  NH1<b>C</b>ARG A 268     4898   2709   3720    472    964    299       N  
ATOM   2058  NH2<b>B</b>ARG A 268      75.794  -4.766  20.144  0.60 30.36           N  
ANISOU 2058  NH2<b>B</b>ARG A 268     4511   3196   3828   1150    633   -562       N  
ATOM   2059  NH2<b>C</b>ARG A 268      75.586  -3.125  18.795  0.40 28.14           N  
ANISOU 2059  NH2<b>C</b>ARG A 268     4497   2583   3613    248   1151    380       N  
</pre>

In this case (5DXX.pdb), there are three different types of the alternative locations, **`A`**, **`B`**, and **`C`** high-lighted with the bold font. However, they distribute with irregular way. For instance, in **sequence 61**, **`A`** and **`B`** appeared, whereas in **sequence 268**, **`B`** and **`C`** emerged. As a result, it is impossible to simply use the ```pdb_info[(pdb_info.Alt_Loc == ' ') | (pdb_info.Alt_Loc == 'A')]``` because that would delete all **`B`** and **`C`** labeled atoms in **sequence 268**!

**Improvement or Debug (Sep. 04, 2017)** By using pandas df.groupby() on the ['Seq_Num', 'ChainID'] columns, we can focus on each specific residue and keep the first alternative location, no matter the first one is 'A' or 'B' or 'C'. The code is show as following
```python
    #### delete the redundant alternate locations, only keep the first apperance
    if altloc:
        groups = pdb_info.groupby(['Seq_Num', 'ChainID'], sort=False)
        pdb_info = groups.apply(lambda x:
                                x.drop_duplicates(subset=["AtomTyp"],
                                                  keep='first')
                                if len(groups['Alt_Loc']) >= 2 else x)
```

* insertion codes
<pre>
ATOM   1258  CD1 ILE A 185       4.002  11.557  18.921  1.00 19.47           C  
ANISOU 1258  CD1 ILE A 185     2567   2632   2200    -66   -252    125       C  
ATOM   1259  N   PRO A 186       6.584  15.226  16.396  1.00 16.95           N  
ANISOU 1259  N   PRO A 186     2324   2351   1766    -93   -218    271       N  
ATOM   1260  CA  PRO A 186       6.984  16.463  15.718  1.00 17.27           C  
ANISOU 1260  CA  PRO A 186     2382   2394   1786   -103   -219    330       C  
ATOM   1261  C   PRO A 186       6.139  17.642  16.167  1.00 19.26           C  
ANISOU 1261  C   PRO A 186     2626   2598   2094    -86   -245    374       C  
ATOM   1262  O   PRO A 186       4.907  17.532  16.301  1.00 18.40           O  
ANISOU 1262  O   PRO A 186     2500   2480   2011    -67   -280    374       O  
ATOM   1263  CB  PRO A 186       6.742  16.159  14.234  1.00 20.29           C  
ANISOU 1263  CB  PRO A 186     2785   2831   2092   -115   -244    345       C  
ATOM   1264  CG  PRO A 186       6.728  14.695  14.124  1.00 25.31           C  
ANISOU 1264  CG  PRO A 186     3421   3497   2701   -115   -240    282       C  
ATOM   1265  CD  PRO A 186       6.252  14.151  15.432  1.00 19.86           C  
ANISOU 1265  CD  PRO A 186     2702   2765   2078   -100   -238    244       C  
ATOM   1266  N   ASP A 186<b>A</b>      6.812  18.768  16.413  1.00 16.88           N  
ANISOU 1266  N   ASP A 186<b>A</b>    2335   2266   1814    -93   -227    410       N  
ATOM   1267  CA  ASP A 186<b>A</b>      6.193  20.046  16.803  1.00 18.33           C  
ANISOU 1267  CA  ASP A 186<b>A</b>    2517   2396   2051    -76   -248    453       C  
ATOM   1268  C   ASP A 186<b>A</b>      5.389  19.957  18.110  1.00 21.71           C  
ANISOU 1268  C   ASP A 186<b>A</b>   2920   2782   2548    -46   -251    420       C  
ATOM   1269  O   ASP A 186<b>A</b>      4.477  20.754  18.337  1.00 23.99           O  
ANISOU 1269  O   ASP A 186<b>A</b>    3201   3034   2879    -21   -276    447       O  
ATOM   1270  CB  ASP A 186<b>A</b>      5.342  20.626  15.640  1.00 20.86           C  
ANISOU 1270  CB  ASP A 186<b>A</b>    2848   2731   2345    -71   -295    510       C  
ATOM   1271  CG  ASP A 186<b>A</b>      6.138  20.870  14.377  1.00 27.21           C  
ANISOU 1271  CG  ASP A 186<b>A</b>    3681   3578   3078   -102   -290    551       C  
ATOM   1272  OD1 ASP A 186<b>A</b>      7.316  21.272  14.485  1.00 27.14           O  
ANISOU 1272  OD1 ASP A 186<b>A</b>    3686   3563   3064   -125   -254    561       O  
ATOM   1273  OD2 ASP A 186<b>A</b>      5.578  20.677  13.277  1.00 34.63           O  
ANISOU 1273  OD2 ASP A 186<b>A</b>    4630   4560   3967   -104   -324    575       O  
ATOM   1274  N   SER A 186<b>B</b>      5.742  18.999  18.983  1.00 16.28           N  
ANISOU 1274  N   SER A 186<b>B</b>    2218   2098   1871    -47   -223    364       N  
ATOM   1275  CA  SER A 186<b>B</b>      5.050  18.813  20.239  1.00 16.16           C  
ANISOU 1275  CA  SER A 186<b>B</b>    2178   2050   1911    -22   -220    332       C  
ATOM   1276  C   SER A 186<b>B</b>      6.014  18.876  21.407  1.00 16.84           C  
ANISOU 1276  C   SER A 186<b>B</b>    2267   2109   2024    -28   -181    302       C  
ATOM   1277  O   SER A 186<b>B</b>      7.167  18.490  21.277  1.00 17.17           O  
ANISOU 1277  O   SER A 186<b>B</b>    2317   2170   2035    -52   -156    289       O  
ATOM   1278  CB  SER A 186<b>B</b>      4.378  17.452  20.244  1.00 17.47           C  
ANISOU 1278  CB  SER A 186<b>B</b>    2323   2250   2066    -18   -229    294       C  
ATOM   1279  OG  SER A 186<b>B</b>      3.785  17.181  21.503  1.00 16.37           O  
ANISOU 1279  OG  SER A 186<b>B</b>    2158   2085   1978      2   -220    264       O  
ATOM   1280  N   LYS A 187       5.518  19.323  22.546  1.00 14.62           N  
ANISOU 1280  N   LYS A 187     1974   1786   1795     -5   -177    290       N  
</pre>
As shown above, the same sequence number (186) labeled with two insertion codes (A and B), however, there are two kinds of residues ASP and SER! The simplest way is to delete the residues labeled by insertion codes. 

* When I save the cleaned results, I found another alignment issue... (data copied from 1BTY.pdb)
<pre>
ATOM      1  N   ILE A  16      35.700  19.589  20.234  1.00 10.94           N  
ATOM      2  CA  ILE A  16      35.550  20.497  19.066  1.00 10.97           C  
ATOM      3  C   ILE A  16      36.807  20.237  18.234  1.00  9.79           C  
ATOM      4  O   ILE A  16      37.894  20.256  18.772  1.00 10.26           O  
ATOM      5  CB  ILE A  16      35.544  21.989  19.514  1.00 11.47           C  
ATOM      6  CG1 ILE A  16      34.399  22.321  20.484  1.00 12.32           C  
ATOM      7  CG2 ILE A  16      35.560  22.968  18.278  1.00 12.30           C  
ATOM      8  CD1 ILE A  16      33.034  22.335  19.785  1.00 13.18           C  
ATOM      9  HA  ILE A  16      34.673  20.230  18.499  1.00 10.47           H  
ATOM     10  HB  ILE A  16      36.473  22.161  20.042  1.00 11.57           H  
ATOM     11 <b>HG12</b> ILE A  16      34.396  21.655  21.334  1.00 11.90           H  
ATOM     12 <b>HG13</b> ILE A  16      34.579  23.313  20.881  1.00 12.02           H  
ATOM     13 <b>HG21</b> ILE A  16      34.717  22.818  17.621  1.00 11.98           H  
ATOM     14 <b>HG22</b> ILE A  16      35.548  23.994  18.620  1.00 12.00           H  
ATOM     15 <b>HG23</b> ILE A  16      36.462  22.839  17.694  1.00 11.56           H  
ATOM     16 <b>HD11</b> ILE A  16      32.786  21.397  19.326  1.00 12.90           H  
ATOM     17 <b>HD12</b> ILE A  16      32.266  22.577  20.509  1.00 12.70           H  
ATOM     18 <b>HD13</b> ILE A  16      33.010  23.114  19.032  1.00 12.55           H
ATOM     19  N   VAL A  17      36.640  20.021  16.964  1.00 11.69           N  
ATOM     20  CA  VAL A  17      37.785  19.760  16.052  1.00 10.93           C  
ATOM     21  C   VAL A  17      37.896  21.020  15.170  1.00  9.18           C  
ATOM     22  O   VAL A  17      36.905  21.499  14.639  1.00 11.67           O  
ATOM     23  CB  VAL A  17      37.466  18.517  15.170  1.00 12.02           C  
ATOM     24  CG1 VAL A  17      38.603  18.296  14.156  1.00 11.39           C  
ATOM     25  CG2 VAL A  17      37.202  17.225  16.050  1.00 13.94           C  
ATOM     26  H   VAL A  17      35.748  20.036  16.564  1.00 11.21           H  
ATOM     27  HA  VAL A  17      38.694  19.634  16.621  1.00 10.27           H  
ATOM     28  HB  VAL A  17      36.577  18.735  14.593  1.00 11.73           H  
ATOM     29 <b>HG11</b> VAL A  17      39.545  18.156  14.663  1.00 11.37           H  
ATOM     30 <b>HG12</b> VAL A  17      38.402  17.438  13.536  1.00 11.87           H  
ATOM     31 <b>HG13</b> VAL A  17      38.686  19.156  13.503  1.00 11.42           H  
ATOM     32 <b>HG21</b> VAL A  17      38.046  16.986  16.679  1.00 12.88           H  
ATOM     33 <b>HG22</b> VAL A  17      36.338  17.370  16.683  1.00 13.55           H  
ATOM     34 <b>HG23</b> VAL A  17      36.989  16.368  15.427  1.00 13.73           H  
ATOM     35  N   GLY A  18      39.101  21.479  15.085  1.00 10.03           N  
ATOM     36  CA  GLY A  18      39.440  22.677  14.271  1.00 12.83           C  
ATOM     37  C   GLY A  18      38.928  24.015  14.824  1.00 14.65           C  
ATOM     38  O   GLY A  18      38.710  24.947  14.072  1.00 14.74           O  
ATOM     39  H   GLY A  18      39.816  21.025  15.573  1.00 10.61           H  
ATOM     40  HA2 GLY A  18      40.513  22.729  14.176  1.00 11.84           H  
ATOM     41  HA3 GLY A  18      39.023  22.532  13.283  1.00 11.64           H  
</pre>
where "HG12", "HG13", "HG21", "HG22", "HG23", "HD11", "HD12", and "HD13" are one character left-shifted compared with the preceding lines. 
**Improvement or Debug (Jul. 31, 2017)** Implemented two printing formats to deal with this issue.

* Usually, PDB files do not contain hydrogen atoms (might due to the highly dynamic of the motion of hydrogen atoms or the limitation of the X-ray resolution). However, in some PDB file, e.g. 5JRY.pdb, hydrogen atoms were recorded. 
<pre>
ATOM      1  N   MET A   1       3.164  22.103 135.939  1.00 28.43           N  
ANISOU    1  N   MET A   1     3558   4003   3241    245    418   -535       N  
ATOM      2  CA  MET A   1       3.182  20.676 135.533  1.00 27.33           C  
ANISOU    2  CA  MET A   1     3398   3863   3124    254    483   -543       C  
ATOM      3  C   MET A   1       3.889  20.519 134.187  1.00 26.06           C  
ANISOU    3  C   MET A   1     3199   3710   2993    137    508   -477       C  
ATOM      4  O   MET A   1       3.671  21.292 133.254  1.00 26.86           O  
ANISOU    4  O   MET A   1     3353   3787   3064    176    487   -391       O  
ATOM      5  CB  MET A   1       1.755  20.132 135.441  1.00 27.72           C  
ANISOU    5  CB  MET A   1     3472   3915   3143    245    541   -550       C  
ATOM      6  <b>H</b>   MET A   1       2.770  22.181 136.733  1.00 34.12           H  
ATOM      7  <b>HA</b>  MET A   1       3.659  20.162 136.203  1.00 32.80           H  
ATOM      8  N   LEU A   2       4.740  19.510 134.085  1.00 23.71           N  
ANISOU    8  N   LEU A   2     2783   3443   2781    -63    570   -535       N  
ATOM      9  CA  LEU A   2       5.389  19.232 132.817  1.00 21.69           C  
ANISOU    9  CA  LEU A   2     2418   3187   2638   -244    572   -551       C  
ATOM     10  C   LEU A   2       4.358  18.763 131.793  1.00 20.93           C  
ANISOU   10  C   LEU A   2     2144   3167   2643   -261    488   -592       C  
ATOM     11  O   LEU A   2       3.268  18.294 132.137  1.00 21.65           O  
ANISOU   11  O   LEU A   2     2148   3333   2746   -381    548   -733       O  
ATOM     12  CB  LEU A   2       6.449  18.148 133.007  1.00 20.95           C  
ANISOU   12  CB  LEU A   2     2388   3018   2555   -362    577   -509       C  
ATOM     13  CG  LEU A   2       7.526  18.465 134.041  1.00 20.97           C  
ANISOU   13  CG  LEU A   2     2473   2948   2547   -376    504   -471       C  
ATOM     14  CD1 LEU A   2       8.523  17.324 134.149  1.00 21.51           C  
ANISOU   14  CD1 LEU A   2     2624   2963   2586   -440    390   -438       C  
ATOM     15  CD2 LEU A   2       8.236  19.754 133.671  1.00 21.33           C  
ANISOU   15  CD2 LEU A   2     2560   2944   2601   -381    439   -472       C  
ATOM     16  <b>H</b>   LEU A   2       4.955  18.978 134.726  1.00 28.45           H  
ATOM     17  <b>HA</b>  LEU A   2       5.819  20.035 132.484  1.00 26.03           H  
ATOM     18  <b>HB2</b> LEU A   2       6.007  17.331 133.287  1.00 25.14           H  
ATOM     19  <b>HB3</b> LEU A   2       6.894  18.000 132.158  1.00 25.14           H  
ATOM     20  <b>HG</b>  LEU A   2       7.110  18.586 134.909  1.00 25.16           H  
ATOM     21 <b>HD11</b> LEU A   2       9.193  17.553 134.812  1.00 25.81           H  
ATOM     22 <b>HD12</b> LEU A   2       8.053  16.519 134.418  1.00 25.81           H  
ATOM     23 <b>HD13</b> LEU A   2       8.943  17.189 133.285  1.00 25.81           H  
ATOM     24 <b>HD21</b> LEU A   2       8.916  19.941 134.336  1.00 25.60           H  
ATOM     25 <b>HD22</b> LEU A   2       8.646  19.648 132.798  1.00 25.60           H  
ATOM     26 <b>HD23</b> LEU A   2       7.588  20.475 133.648  1.00 25.60           H  
</pre>
**Improvement (Nov. 16, 2017)** In this case, I added a new option in my PDB_cleaner script enabling the users to choose whether remove all hydrogen atoms or not.

## How to run this script
This script can be run in both Linux and Windows system. The command is shown below,

$`python pdb_cleaner.py`

Then, the program will ask you to specified the directory that the PDB files located, and how to deal with multiple chains (keep all the chains or just one of them).

If you choose "one", the program will choose the longest chain in the PDB file (if all chains have the same length, the first chain will be kept).

* Workflow:

1. Collect all the PDB files in the given directory;
2. In each PDB file, check the following items:

    2.1. alternate locations;

    2.2. non-standard amino acid residues;

    2.3. negative sequence numbers (less important);
        
    2.4. sequence gaps;
        
    2.5. insertion code;
        
    2.6. multiple chains;
        
    2.7. hydrogen atoms;
        
    2.8. __** to do: missing atoms **__
        
3. Clean the PDB files if the aforementioned items exist,
   with following options if protein has multiple chains;
    
    3.1. remove hydrogein, if the user specified "y";
    
    3.2. keep all chains if the user specified "all";
        
    3.3. keep the longest chain (or the 1st chain, if all
        chains have the same length), if the user specified "one".
        
4. Save the cleaned PDB files one by one;
    
5. Save the summary report.
