#Seeking extrema after Gaussian denoising.
#Command line usage: You can give THREE PROPER POSITIVE INTEGERS to serve as the MINIMUM and MAXIMUM window for Gaussian denoising
# and the MAXIMUM number of spaced nucleotides in merging predicted sites. Default values are 2, 2 and 1.
#e.g. --> 'Rscript getGaussianDenoisedExtrema.R' or 'Rscript getGaussianDenoisedExtrema.R 1 5 2' or 'Rscript getGaussianDenoisedExtrema.R 2 2 1' or ...
#knownSites = list(1:4, 13:22, 34:36, 49:51, 53:61, 63:65, 72:76)	#tRNA(Lys)(PDB#1FIR)
#knownSites = list(1:4, 23:25, 42:45)	#Diels-Alder Ribozyme(PDB#1YKV)
#knownSites = list(13:28, 31:42, 62:89, 151:152, 155:159, 170, 172, 183:184, 186:189, 194:209, 213:215, 237:239, 253:266, 298:307, 311:313, 323:325, 341:353)	#RNase P
#knownSites = list(1:4, 12:27, 34:45, 64:77, 84:95, 106:111, 301:332)	#7SK
knownSites = list()	#UNKNOWN
knownSitesNum = length(knownSites)

args = commandArgs(TRUE)
if (length(args) == 0){
	minWin = 2	#Output equals input.
	maxWin = 2
	spacedNt = 1
}
if (length(args) > 0){
	minWin = as.integer(args[1])
	maxWin = as.integer(args[2])
	spacedNt = as.integer(args[3])
}

outDir = getwd()
dirs = c('D:/Desktop/rna/code/RLBind/v2/reproduction/Rsite2/PS/SS_NDC', 'D:/Desktop/rna/code/RLBind/v2/reproduction/Rsite2/PS/SS_NDS')
# # 如果目录不存在，可以创建
# dir.create("./PS/SS_NDS", recursive = TRUE)
# dir.create("./PS/SS_NDC", recursive = TRUE)

for (d in dirs){
	setwd(d)
	flist = dir()
	ndlist = grep('((nds)|(ndc))\\.txt', flist, value=T)
	for (i in ndlist){
		name = gsub('\\.', paste('_gd_extrema_', minWin, '-', maxWin, '\\.', sep=''), i)
		if (!(name %in% flist)){
			nd = read.table(i, header=T)
			nts = nd[,1]
			signal = nd[,2]
			ntNum = length(signal)
			Q1 = quantile(signal,0.25)	#The lower quartile.
			Q3 = quantile(signal,0.75)	#The upper quartile.
			
			for (win in minWin:maxWin){
				paddedSig = c(rep(signal[1], ceiling((win-1)/2)), signal, rep(rev(signal)[1], floor((win-1)/2)))	#Pad the signal by corresponding window size at both ends.
				gaussianDenoiser = exp(-0.5*seq(-(win-1)/2, (win-1)/2)^2)
				gaussianDenoiser = gaussianDenoiser/sum(gaussianDenoiser)
				denoisedSig = convolve(paddedSig, gaussianDenoiser, type='o')
				denoisedSig = round(denoisedSig, 2)
				left = win
				right = length(denoisedSig)-(win-1)
				denoisedSig = denoisedSig[left:right]	#Keep the length of signal as before.

				#Extrema seeking.
				extremaTF = (abs(diff(sign(diff(denoisedSig)))) == 2)
				extremaLoc = c()
				for (j in 1:length(extremaTF)){
					if (any(is.na(extremaTF[j]))) {
					  print(paste("NA found in extremaTF at file:", i, "win:", win))
					  print(head(denoisedSig))
					}
					if (extremaTF[j]){
						extremaLoc = c(extremaLoc, j+1)
					}
				}
				#extremaLoc = unique(c(1, extremaLoc, ntNum))
				#extremaLoc = extremaLoc[which(signal[extremaLoc] <= Q1 | signal[extremaLoc] >= Q3)]

				if (!(1 %in% extremaLoc) && ((signal[1] <= Q1) || (signal[1] >= Q3))){
					extremaLoc = c(1, extremaLoc)
				}
				if (!(ntNum %in% extremaLoc) && ((signal[ntNum] <= Q1) || (signal[ntNum] >= Q3))){
					extremaLoc = c(extremaLoc, ntNum)
				}

				extremaNt = nts[extremaLoc]
				extremaVal = denoisedSig[extremaLoc]
				extremaRank = rank(extremaVal)
				extremaNum = length(extremaLoc)
				extrema = data.frame(extremaRank, extremaLoc, extremaNt, extremaVal)
				extrema = extrema[order(extrema$extremaRank),]

				#Merge multiple predicted nucleotides into one if they are close.
				predictedNts = sort(extremaLoc)
				predictedSites = list()
				tmp = predictedNts[1]
				predictedSitesNum = 1
				for (k in 2:extremaNum){
					if (predictedNts[k]-predictedNts[k-1] > spacedNt+1){
						predictedSitesNum = predictedSitesNum+1
						predictedSites = append(predictedSites, list(tmp))
						tmp = predictedNts[k]
					}
					else{
						tmp = c(tmp, predictedNts[k])
					}
				}
				predictedSites = append(predictedSites, list(tmp))
				
				#If there are known functional sites, then evaluate the prediction.
				if (knownSitesNum > 0){
					matchSites = c()
					matchSitesNum = 0
					site = c(0)
					for (i in 1:knownSitesNum){
						for (j in 1:length(knownSites[[i]])){
							if (knownSites[[i]][j] %in% predictedNts){
								if (i != site[length(site)]){
									site = c(site, i)
									matchSitesNum = matchSitesNum+1
									matchSites = c(matchSites, knownSites[[i]][j])
								}
								else{
									site = c(site, i)
									matchSites = c(matchSites, knownSites[[i]][j])
								}
							}
						}
					}
					matchPredictedSitesNum = matchSitesNum
					if (length(matchSites) > 1){
						matchPredictedSitesNum = 1
						for (k in 2:length(matchSites)){
							if (matchSites[k]-matchSites[k-1] > spacedNt+1){
								matchPredictedSitesNum = matchPredictedSitesNum+1
							}
						}
					}
					sensitivity = round(matchSitesNum/knownSitesNum, 2)*100
					PPV = round(matchPredictedSitesNum/predictedSitesNum, 2)*100
					rate = paste(sensitivity, '%(', matchSitesNum, '/', knownSitesNum, '),',
					 PPV, '%(', matchPredictedSitesNum, '/', predictedSitesNum, ')\n')
				}
				else{
					rate=paste('Predicted sites number:', predictedSitesNum, '\nPredicted nucleotides number:', extremaNum, '\n')
				}

				write(paste('Gaussian Denoiser Window Size:',win),
				 file=name, append=T)
				write(rate, file=name, append=T)
				siteID = 1
				for (predSite in predictedSites){
					write(paste('Site#', siteID, ':', sep=''), file=name, append=T)
					siteID = siteID+1
					tmp = list()
					for (predNt in predSite){
						tmp = append(tmp, predNt)
					}
					write.table(tmp, file=name, sep=',', append=T,
					 quote=F, row.names=F, col.names=F)
				}
				write('\nextremaRank\textremaLoc\textremaNt\textremaNDS_gd', file=name, append=T)
				write.table(extrema, file=name, sep='\t', append=T,
				 quote=F, row.names=F, col.names=F)
				write('\nGaussian Denoised Signal:', file=name, append=T)
				write.table(data.frame(seq(1,ntNum),nts,denoisedSig), file=name, sep='\t', append=T,
				 quote=F, row.names=F, col.names=F)
				write('====================\n', file=name, append=T)
			}
		}
	}
	setwd(outDir)
}
