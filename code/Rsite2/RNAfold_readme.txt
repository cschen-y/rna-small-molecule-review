RNAfold是预测序列二级结构的软件。

本地版下载地址：
http://www.tbi.univie.ac.at/~ivo/RNA/windoze/

在上面个的网址中还有其他预测结构的软件可供下载~

现在主要介绍一下本地下载版的使用方法：

1.不能够双击应用程序，若双击RNAfold.exe应用,只能手动输入序列，而且只能一条一条输入，极其麻烦

2.应该在Dos环境下，到达应用程序的当前目录。
  
  然后输入命令：
  
  RNAfold.exe  <seq.fasta   >result.txt

  其中“<”代表序列以文件形式作为输入，这就允许输入多个序列进行结构预测，而不单单是一个。若没有这个"<", 则无结果。
      
      “>”代表结果以文件形式输出，最终的预测的结构在result.txt文件中。若没有这个符号，则会输出在dos环境下。

当然还有online服务的网址：http://rna.tbi.univie.ac.at/cgi-bin/RNAfold.cgi

PS: 在使用中发现如果仅用  RNAfold.exe  <seq.fasta   >result.txt 会得到结构文件，但是也会出现大量的.ps文件，因此要使得得到的仅仅是结构文件的话，要使用命令： RNAfold.exe  -noPS <seq.fasta   >result.txt 






RNAFOLD
ViennaRNA (l) Return to Main Contents 
--------------------------------------------------------------------------------

NAME
RNAfold - calculate secondary structures of RNAs 
SYNOPSIS
RNAfold [-p[0|2]] [-C] [-T temp] [-4] [-d[0|1|2|3]] [-noLP] [-noGU] [-noCloseGU] [-e 1|2] [-P paramfile] [-nsp pairs] [-S scale] [-circ] [-MEA [gamma]] 
DESCRIPTION
RNAfold reads RNA sequences from stdin, calculates their minimum free energy (mfe) structure and prints to stdout the mfe structure in bracket notation and its free energy. If the -p option was given it also computes the partition function (pf) and base pairing probability matrix, and prints the free energy of the thermodynamic ensemble, the frequency of the mfe structure in the ensemble, and the ensemble diversity to stdout. 
It also produces PostScript files with plots of the resulting secondary structure graph and a "dot plot" of the base pairing matrix. The dot plot shows a matrix of squares with area proportional to the pairing probability in the upper right half, and one square for each pair in the minimum free energy structure in the lower left half. For each pair i-j with probability p>10E-6 there is a line of the form
i j sqrt(p) ubox
in the PostScript file, so that the pair probabilities can be easily extracted.

Sequences are read in a simple text format where each sequence occupies a single line. Each sequence may be preceded by a line of the form
> name
to assign a name to the sequence. If a name is given in the input
 PostScript files "name_ss.ps" and "name_dp.ps" are produced for the structure and dot plot, respectively. Otherwise the file names default to rna.ps and dot.ps. Existing files of the same name will be overwritten.
The input format is similar to fasta except that even long sequences may not be interrupted by line breaks, and the header lines are optional. The program will continue to read new sequences until a line consisting of the single character @ or an end of file condition is encountered.

OPTIONS
-p 
Calculate the partition function and base pairing probability matrix in addition to the mfe structure. Default is calculation of mfe structure only. In addition to the MFE structure we print a coarse representation of the pair probabilities in form of a pseudo bracket notation, followed by the ensemble free energy, as well as the centroid structure derived from the pair probabilities together with its free energy and distance to the ensemble. Finally it prints the frequency of the mfe structure, and the structural diversity (mean distance between the structures in the ensemble). See the description of pf_fold() and mean_bp_dist() and centroid() in the RNAlib documentation for details.
Note that unless you also specify -d2 or -d0, the partition function and mfe calculations will use a slightly different energy model. See the discussion of dangling end options below. 
-p0 
Calculate the partition function but not the pair probabilities, saving about 50% in runtime. Prints the ensemble free energy -kT ln(Z). 
-p2 
In addition to pair probabilities compute stack probabilities, i.e. the probability that a pair (i,j) and the immediately interior pair (i+1,j-1) are formed simultaneously. A second postscript dot plot called "name_dp2.ps", or "dot2.ps" (if the sequence does not have a name), is produced that contains pair probabilities in the upper right half and stack probabilities in the lower left. 
-MEA [gamma] 
Calculate an MEA (maximum expected accuracy) structure, where the expected accuracy is computed from the pair probabilities: each base pair (i,j) gets a score 2*gamma*p_ij and the score of an unpaired base is given by the probability of not forming a pair. The parameter gamma tunes the importance of correctly predicted pairs versus unpaired bases. Thus, for small values of gamma the MEA structure will contain only pairs with very high probability. The default value is gamma=1. Using -MEA implies -p for computing the pair probabilities. 
-C 
Calculate structures subject to constraints. The program reads first the sequence, then a string containing constraints on the structure encoded with the symbols: | (the corresponding base has to be paired x (the base is unpaired) < (base i is paired with a base j>i) > (base i is paired with a base j<i) and matching brackets ( ) (base i pairs base j) With the exception of "|", constraints will disallow all pairs conflicting with the constraint. This is usually sufficient to enforce the constraint, but occasionally a base may stay unpaired in spite of constraints. PF folding ignores constraints of type "|". 
-T temp 
Rescale energy parameters to a temperature of temp C. Default is 37C. 
-4 
Do not include special stabilizing energies for certain tetra-loops. Mostly for testing. 
-d[0|1|2|3] 
How to treat "dangling end" energies for bases adjacent to helices in free ends and multi-loops: With (-d1) only unpaired bases can participate in at most one dangling end, this is the default for mfe folding but unsupported for the partition function folding. With -d2 this check is ignored, dangling energies will be added for the bases adjacent to a helix on both sides in any case; this is the default for partition function folding (-p). -d or -d0 ignores dangling ends altogether (mostly for debugging).
With -d3 mfe folding will allow coaxial stacking of adjacent helices in multi-loops. At the moment the implementation will not allow coaxial stacking of the two interior pairs in a loop of degree 3 and works only for mfe folding.
Note that by default (as well as with -d1 and -d3) pf and mfe folding treat dangling ends differently. Use -d2 in addition to -p to ensure that both algorithms use the same energy model. 
-noLP 
Produce structures without lonely pairs (helices of length 1). For partition function folding this only disallows pairs that can only occur isolated. Other pairs may still occasionally occur as helices of length 1. 
-noGU 
Do not allow GU pairs. 
-noCloseGU 
Do not allow GU pairs at the end of helices. 
-e 1|2 
Rarely used option to fold sequences from the artificial ABCD... alphabet, where A pairs B, C-D etc. Use the energy parameters for GC (-e 1) or AU (-e 2) pairs. 
-P <paramfile> 
Read energy parameters from paramfile, instead of using the default parameter set. A sample parameter file should accompany your distribution. See the RNAlib documentation for details on the file format. 
-nsp pairs 
Allow other pairs in addition to the usual AU,GC,and GU pairs. pairs is a comma separated list of additionally allowed pairs. If a the first character is a "-" then AB will imply that AB and BA are allowed pairs. e.g. RNAfold -nsp -GA will allow GA and AG pairs. Nonstandard pairs are given 0 stacking energy. 
-S scale 
In the calculation of the pf use scale*mfe as an estimate for the ensemble free energy (used to avoid overflows). The default is 1.07, useful values are 1.0 to 1.2. Occasionally needed for long sequences. 
-circ 
Assume a circular (instead of linear) RNA molecule. 
-noPS 
Do not produce postscript drawing of the mfe structure. 
REFERENCES
The calculation of mfe structures is based on dynamic programming algorithm originally developed by M. Zuker and P. Stiegler. The partition function algorithm is based on work by J.S. McCaskill. The energy parameters are taken from:
D.H. Mathews, J. Sabina, M. Zuker and H. Turner "Expanded Sequence Dependence of Thermodynamic Parameters Provides Robust Prediction of RNA Secondary Structure" JMB, 288, pp 911-940, 1999
A. Walter, D Turner, J Kim, M Lyttle, P M[:u]ller, D Mathews, M Zuker "Coaxial stacking of helices enhances binding of oligoribonucleotides.." PNAS, 91, pp 9218-9222, 1994 
If you use this program in your work you might want to cite:

I.L. Hofacker, W. Fontana, P.F. Stadler, S. Bonhoeffer, M. Tacker, P. Schuster (1994) Fast Folding and Comparison of RNA Secondary Structures. Monatshefte f. Chemie 125: 167-188
M. Zuker, P. Stiegler (1981) Optimal computer folding of large RNA sequences using thermodynamic and auxiliary information, Nucl Acid Res 9: 133-148
J.S. McCaskill (1990) The equilibrium partition function and base pair binding probabilities for RNA secondary structures, Biopolymers 29: 1105-1119
I.L. Hofacker & P.F. Stadler (2006) Memory Efficient Folding Algorithms for Circular RNA Secondary Structures, Bioinformatics
A.F. Bompf??newerer, R. Backofen, S.H. Bernhart, J. Hertel, I.L. Hofacker, P.F. Stadler, S. Will (2007) "Variations on RNA Folding and Alignment: Lessons from Benasque" J. Math. Biol.
D. Adams (1979) The hitchhiker's guide to the galaxy, Pan Books, London

VERSION
This man page documents version 1.8.5 Vienna RNA Package. 
AUTHORS
Ivo L Hofacker, Walter Fontana, Sebastian Bonhoeffer, Peter F Stadler. 
BUGS
If in doubt our program is right, nature is at fault. Comments should be sent to rna@tbi.univie.ac.at. 
--------------------------------------------------------------------------------
This document was created by man2html, using the manual pages.
Time: 07:19:16 GMT, February 23, 2011 