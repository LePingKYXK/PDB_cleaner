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

* HETATM

* non-standard amino acid residues (data copied from 2o2x.pdb)
<pre>
ATOM    821  OD1 ASP A 112      25.580  11.019  35.906  1.00 12.28           O  
ATOM    822  OD2 ASP A 112      24.586   9.016  35.848  1.00 11.81           O  
HETATM  823  N   <b>MSE</b> A 113      25.018  10.050  30.641  1.00  9.26           N  
HETATM  824  CA  <b>MSE</b> A 113      25.494  10.026  29.262  1.00  9.59           C  
HETATM  825  C   <b>MSE</b> A 113      24.291   9.758  28.359  1.00  8.63           C  
HETATM  826  O   <b>MSE</b> A 113      23.362   9.026  28.750  1.00  9.51           O  
HETATM  827  CB  <b>MSE</b> A 113      26.563   8.959  29.078  1.00  8.81           C  
HETATM  828  CG  <b>MSE</b> A 113      27.157   8.896  27.700  1.00  8.23           C  
HETATM  829 SE   <b>MSE</b> A 113      28.681   7.732  27.499  0.75 12.65          SE  
HETATM  830  CE  <b>MSE</b> A 113      30.013   8.895  28.258  1.00 18.85           C  
ATOM    831  N   VAL A 114      24.306  10.362  27.178  1.00  7.56           N  
ATOM    832  CA  VAL A 114      23.308  10.072  26.129  1.00  7.93           C  
</pre>
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
As shown above, the same sequence number (186) labeled with two insertion codes (A and B), however, there are two kinds of residues ASP and SER!

* When I save the cleaned results, I found another alignment issue... (data copied from 1bty.pdb)
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
