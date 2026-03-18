#This script is designed to parse PostScript files ('.ps') produced by RNAfold.
#Or be sure that '/coor [...] def' rows are recorded in the PostScript files ('.ps').
#The results are the two metrics calculated and separately written into files for each RNA.
#The results are saved into the folder 'PS/SS_NDC' and the folder 'PS/SS_NDS' respectively.
#Command line usage: The working directory is where the script runs.
#e.g. --> 'python getNtDis2.py[ test.fst]' (If you give the RNA primary sequences in a FASTA file ('.fst').)
#e.g. --> 'python getNtDis2.py ps' (If you give the RNA secondary structure PostScript files ('.ps') in the folder 'PS'.)

#!/usr/bin/env python

import os
import sys
import re
import csv
import shutil
import subprocess

class Coor(object):
    def __init__(self, name, coor):
        self.name = name
        self.coor = coor
        self.coorNum = len(self.coor)
        try:
            self.dim = len(self.coor[0])
        except IndexError as diag:
            print(diag)
        tmp = list()
        for i in range(self.dim):
            tmp.append(sum([j[i] for j in self.coor])/self.coorNum)
        self.center = tuple(tmp)

    def __repr__(self):
        return self.name

    def getOneDis(self, point1=0, point2=1):
        try:
            tmp = 0
            for i in range(self.dim):
                tmp += pow((self.coor[point1][i]-self.coor[point2][i]), 2)
            dis = pow(tmp, 0.5)
            return dis
        except IndexError as diag:
            print(diag)

    def getAllDis(self):
        self.allDis = list([('nt1', 'nt2', 'dis')])
        for i in range(self.coorNum):
            for j in range((i+1), self.coorNum):
                tmp = self.getOneDis(i, j)
                self.allDis.append(((i+1), (j+1), tmp))
        return self.allDis

    def getOneNDS(self, point=0):
        try:
            nds = 0
            for i in range(self.coorNum):
                nds += self.getOneDis(point, i)
            return nds
        except IndexError as diag:
            print(diag)

    def getAllNDS(self):
        self.allNDS = list([('nt', 'nds')])
        for i in range(self.coorNum):
            nds = self.getOneNDS(i)
            self.allNDS.append(((i+1), nds))
        return self.allNDS

    def getOneNDC(self, point=0):
        try:
            tmp = 0
            for i in range(self.dim):
                tmp += pow((self.coor[point][i]-self.center[i]), 2)
            ndc = pow(tmp, 0.5)
            return ndc
        except IndexError as diag:
            print(diag)

    def getAllNDC(self):
        self.allNDC = list([('nt', 'ndc')])
        for i in range(self.coorNum):
            tmp = self.getOneNDC(i)
            self.allNDC.append(((i+1), tmp))
        return self.allNDC

class PS_RNAFOLD(object):
    def __init__(self, fileName):
        self.fName = fileName

    def __repr__(self):
        return self.fName

    def getCoor(self):
        pattern = re.compile(r'^\[(?P<x>-?\d+\.\d+) (?P<y>-?\d+\.\d+)\]$')
        self.coor = list()
        flag = 0
        with open(self.fName, 'r') as f:
            for l in f:
                if (l == '/coor [\n') and (not flag):	#Before reaching the 'coor' rows.
                    flag = 1
                if (l == '] def\n') and flag:	#Finished the 'coor' rows.
                    flag = 0
                if (l != '/coor [\n') and flag:
                    match = pattern.search(l)
                    self.coor.append((float(match.group('x')),float(match.group('y'))))
        return self.coor

def getNDCorS2(fName='test.fst', mode='ndc'):
    pattern = re.compile(r'\..+')
    if (not os.path.exists('PS')):
        os.makedirs('PS')
    if ('.fst' in fName):
        os.chdir('PS')
        if (not os.path.exists(fName)):
            os.chdir('..')
            shutil.copy(fName, 'PS')
            os.chdir('PS')
        if (not os.path.exists('RNAfold.exe')):
            os.chdir('..')
            shutil.copy('RNAfold.exe', 'PS')
            os.chdir('PS')
        cmd = 'RNAfold.exe < %s > %s' % (fName, fName.replace('.fst', '_ss.txt'))
        sp = subprocess.call(cmd, shell=True)
    else:
        os.chdir('PS')
    folders = dict(zip(['ndc','nds'], ['SS_NDC/','SS_NDS/']))
    suffixes = dict(zip(['ndc','nds'], ['_ndc.txt','_nds.txt']))
    if os.path.exists(folders[mode]):
        shutil.rmtree(folders[mode])
    os.makedirs(folders[mode])
    fs = os.listdir(os.getcwd())
    psfs = [f for f in fs if '.ps' in f]
    for psf in psfs:
        psid = pattern.sub('', psf)
        ps = PS_RNAFOLD(psf)
        coor = Coor(psf, ps.getCoor())
        rows = dict(zip(['ndc','nds'], [coor.getAllNDC(),coor.getAllNDS()]))
        try:
            with open(folders[mode]+psid+suffixes[mode], 'w') as fo:
                writer = csv.writer(fo, dialect='tab')
                writer.writerows(rows[mode])
        except IOError as diag:
            print(diag)
    os.chdir('..')

def main():
    csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONNUMERIC, lineterminator='\n')
    if (len(sys.argv) == 1):
        getNDCorS2(fName='test.fst', mode='ndc')
        getNDCorS2(fName='test.fst', mode='nds')
    elif ('.fst' in sys.argv[1]):
        getNDCorS2(fName=sys.argv[1], mode='ndc')
        getNDCorS2(fName=sys.argv[1], mode='nds')
    elif (sys.argv[1] == 'ps'):
        getNDCorS2(fName='', mode='ndc')
        getNDCorS2(fName='', mode='nds')

if __name__ == '__main__':
    main()
