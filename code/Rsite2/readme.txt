Rsite2

STEP 1: getNtDis2.py
	This script is designed to parse PostScript files ('.ps') produced by RNAfold.
	Or be sure that '/coor [...] def' rows are recorded in the PostScript files ('.ps').
	The results are the two metrics calculated and separately written into files for each RNA.
	The results will be saved into the folder 'PS/SS_NDC' and the folder 'PS/SS_NDS' respectively.
	
	Command line usage: The working directory is where the script runs.
	Option 1: e.g. --> 'python getNtDis2.py[ test.fst]' (If you give the ncRNA primary sequences in a FASTA file ('.fst').)
		'RNAfold.exe' must be in the same folder as 'getNtDis2.py'. 
		'test.fst' is provided for testing this option. It contains the primary sequences of two ncRNAs (PDB ID: 1FIR&1YKV).
	Option 2: e.g. --> 'python getNtDis2.py ps' (If you give the ncRNA secondary structure PostScript files ('.ps') in the folder 'PS'.)
		The two PostScript files ('1FIR_ss.ps'&'1YKV_ss.ps') in the folder 'PS' produced by RNAfold are provided for testing this option.

STEP 2: getGaussianDenoisedExtrema.R
	Seeking extrema after Gaussian denoising.
	The results will also be saved into the folder 'PS/SS_NDC' and the folder 'PS/SS_NDS' respectively.
	
	Command line usage: You can give THREE PROPER POSITIVE INTEGERS to serve as the MINIMUM and MAXIMUM window for Gaussian denoising and the MAXIMUM number of spaced nucleotides in merging predicted sites. Default values are 2, 2 and 1.
	e.g. --> 'Rscript getGaussianDenoisedExtrema.R' or 'Rscript getGaussianDenoisedExtrema.R 1 5 2' or 'Rscript getGaussianDenoisedExtrema.R 2 2 1' or ...
