#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 15:53:37 2021

@author: cathy
"""
from os import path,popen,mkdir
sCurrentFilePath = path.dirname(path.abspath(__file__)) + '/'

import random
import primer3
from datetime import datetime
from traceback import format_exc
import pandas as pd

import primer3Setting

import logging
log = logging.getLogger('mylog')

dNumBase = {'0': 'A', '1': 'C', '2': 'G', '3': 'T'}

dBaseNum = {'A': '0', 'C': '1', 'G': '2', 'T': '3'}

dMap = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        '0': '3', '1': '2', '2': '1', '3': '0'}

lBase = ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G']

iWordSize = 12

def getTimeStr():
    '''
    get time string for file name
    '''
    s = datetime.now()
    return '{}'.format(s.hour).zfill(2) + '{}'.format(s.minute).zfill(2) + '{}'.format(s.second).zfill(2)

class PrimerDesign:
    def __init__(self, lTemplateFasta, iPrimerLen=20, iPrimer3PairNum=5, iRandTempNumInCycle=5,
                 iRandTempLen=200, iBreakNum=5, iBreakSec=180, sBaseDir = sCurrentFilePath, bTimeTempName = True, bLenStrict=True):
        """
        PrimerDesign Class
        @param lTemplateFasta: path list of the template files
        @param iPrimerLen: optimal primer length
        @param iPrimer3PairNum: The number of primer pairs designed by Primer3 for each template (the maximum value, there may not be so many designed)
        @param iRandTempNumInCycle: The number of random sequence template strips used in each cycle. Because each template will use primer3 to design a set of primers,
                                    the size of the value directly affects the time of each cycle, and it does not need to be set too large
        @param iRandTempLen: Random sequence length used to design primers, the default minimum of primer3 is 100 bases
        @param iBreakNum: Primer pairs that jump out of the loop (and time co-limited)
        @param iBreakSec: Run seconds to jump out of loop (with primer pair limit)
        @param sBaseDir:
        @param bTimeTempName: whether to name temp files with a timestamp
        @param bLenStrict: Whether the primer length is strictly limited
        """
        '''
        PrimerDesign Class
        :param lTemplateFasta: path list of the template files
        :param iPrimerLen: optimal primer length
        :param iPrimer3PairNum: The number of primer pairs designed by Primer3 for each template (the maximum value, there may not be so many designed)
        :param iRandTempNumInCycle: The number of random sequence template strips used in each cycle. Because each template will use primer3 to design a set of primers,
                                    the size of the value directly affects the time of each cycle, and it does not need to be set too large
        :param iRandTempLen: Random sequence length used to design primers, the default minimum of primer3 is 100 bases
        :param iBreakNum: Primer pairs that jump out of the loop (and time co-limited)
        :param iBreakSec: Run seconds to jump out of loop (with primer pair limit)
        '''
        if sBaseDir[-1]!='/' and sBaseDir[-1]!='\\':
            sBaseDir += '/'
        sBaseDir += 'tmp/'
        if not path.isdir(sBaseDir):
            mkdir(sBaseDir)
        if bTimeTempName:
            self._tempFileName = getTimeStr()
        else:
            self._tempFileName = 'c'
            
        self._bLenStrict = bLenStrict
        self._lTemplateFasta = lTemplateFasta
        self._iPrimerLen = iPrimerLen
        self._iPrimer3PairNum = iPrimer3PairNum
        self._iRandTempNumInCycle = iRandTempNumInCycle
        self._iRandTempLen = iRandTempLen
        self._iBreakNum = iBreakNum
        self._iBreakSec = iBreakSec
        self._basePath = sBaseDir

        self._sAllTemplateFastFile = sBaseDir + self._tempFileName + '_allTemplate.fasta'  # Processed template fasta file
        self._sPrimerFastaFile = sBaseDir + self._tempFileName + '_primer.fasta'  # Primer file with 3-end verification
        self._sBlastOutputFile = sBaseDir + self._tempFileName + '_blast_output.txt'
        self._lTemplateAll = list()  # quaternary template sequence
        self._setTemp6 = set()  # Template 5-terminal six bases and 3-terminal reverse complementary six bases all combinations (decimal numbers)
        self._lPrimerPair = list()  # 3-end-tested primer pair designed by primer3
        self._dfOutPut = pd.DataFrame()  # Final sequenced primer pair
        self._succeed = True # False represents actual no primer design results, but uses the last set of results

    def getPrimer(self):
        '''
        main function
        :return:
        '''
        try:
            self._preProcess()
            self._primerDesign(self._iRandTempNumInCycle, self._iRandTempLen, self._iPrimer3PairNum,
                               self._iBreakNum, self._iBreakSec)
            self._getPrimerFasta()
            self._blastCheck()

            if len(self._dfOutPut) != 0:
                sOutPut = self._basePath + self._tempFileName + '_primerOutput_{}.csv'.format(
                    datetime.strftime(datetime.today(), "%Y%m%d%H%M%S"))
                self._dfOutPut.to_csv(sOutPut)
                log.info("output file: {}".format(sOutPut))
                log.info('A total of {} pairs of primers, the best recommended primer pairs are: {} {}'.format(len(self._dfOutPut.index), self._dfOutPut.iloc[0].left,
                                                          self._dfOutPut.iloc[0].right))
            else:
                log.warning("Eligible primers have not yet been designed to increase cycle time for redesign")
            return self._dfOutPut
        except Exception as e:
            log.error("getPrimer failed : {}".format(e) + str(format_exc()))
            raise

    def _preProcess(self):
        '''
        Working with template files
        :return:
        '''
        try:
            # Splice all template fasta files
            getAllTemplateFasta(self._lTemplateFasta,
                                self._sAllTemplateFastFile)
            # get the quaternary sequence
            self._lTemplateAll = getQuaStr(self._sAllTemplateFastFile)
            # Get the 5' 6 base and 3' reverse complementary six bases of the template sequence
            self._setTemp6 = self._getTemplate6()
            # use blast to generate template library
            log.info(runBlastDB(self._sAllTemplateFastFile))
        except Exception as e:
            log.error("preProcess failed : {}".format(e) + str(format_exc()))
            raise

    def _getTemplate6(self):
        '''
        Get template 5' six bases and 3' reverse complementary six bases, save as decimal numbers
        @return: list of numbers
        '''
        setTemp6 = set([temp[:6] for temp in self._lTemplateAll] +
                       [getOther(temp[-6:]) for temp in self._lTemplateAll])
        self._setTemp6 = {int(temp, 4) for temp in setTemp6}
        return self._setTemp6

    def _primerDesign(self, iRandTempNumInCycle, iRandTempLen, iPrimer3PairNum, iBreakNum, iBreakSec):
        '''
        Get 3-terminal verified primer pairs
        :param iRandTempNumInCycle: The number of random sequence templates used in each cycle. Because each template will use primer3 to design a set of primers,
                the size of the value directly affects the time of each cycle, and it does not need to be set too large
        :param iRandTempLen: The random sequence length used to design primers, the default minimum of primer3 is 100 bases
        :param iPrimer3PairNum: The number of primer pairs designed by Primer3 for each template (the maximum value, it may not be able to design so many)
        :param iBreakNum: primer pair to break out of the loop (limited with time)
        :param iBreakSec: Number of seconds to run to break out of the loop (shared with primer pairs)
        :return: list of primer pairs
        '''
        try:
            # Generate and verify primer pairs using random sequences
            self._lPrimerPair = list()
            tm = datetime.now()
            cycle = 0
            while(1):
                cycle += 1

                # Each cycle uses 5 random sequences, each generating up to 10 pairs of primers, for a total of 50 pairs
                lTempPrimerPairTemp = getPrimerPair(
                    iRandTempNumInCycle, iRandTempLen, self._iPrimerLen, iPrimer3PairNum)

                # Check primer length
                if self._bLenStrict:
                    log.debug("check length")
                    lTempPrimerPair = list()
                    for (lPrimer,rPrimer) in  lTempPrimerPairTemp:
                        if len(lPrimer)==self._iPrimerLen and len(rPrimer)==self._iPrimerLen:
                            lTempPrimerPair.append((lPrimer,rPrimer))
                else:
                    lTempPrimerPair = lTempPrimerPairTemp
                # Check primer 3' six bases
                self._lPrimerPair += getPrimerPairList(
                    lTempPrimerPair, self._setTemp6)

                iSec = (datetime.now() - tm).seconds
                log.info('In the {} cycle, {} pairs of primers have been obtained, and it has taken {}s'.format(
                    cycle, len(self._lPrimerPair), iSec))
                if len(self._lPrimerPair) >= iBreakNum or iSec >= iBreakSec:
                    # Use the last set of results to fill when there are no design results
                    if len(self._lPrimerPair) == 0:
                        self._lPrimerPair += lTempPrimerPair
                        self._succeed = False
                    break

            if len(self._lPrimerPair) == 0:
                log.error("No primer")
            return self._lPrimerPair

        except Exception as e:
            log.error("primerDesign failed : {}".format(e) + str(format_exc()))
            raise

    def _getPrimerFasta(self):
        '''
        Generate primer fasta files for blast alignment
        :return:
        '''
        try:
            lFasta = list()
            for index, (sLeft, sRight) in enumerate(self._lPrimerPair):
                lFasta.append('>primer_{}\n'.format(index * 2))
                lFasta.append(sLeft + '\n')
                lFasta.append('>primer_{}\n'.format(index * 2 + 1))
                lFasta.append(sRight + '\n')

            with open(self._sPrimerFastaFile, 'w') as file_output:
                file_output.writelines(lFasta)
            return True
        except Exception as e:
            log.error("getPrimerFasta failed : {}".format(e) + str(format_exc()))
            raise

    def _blastCheck(self):
        '''
        Generate primer fasta files and compare them using blast
        :return:
        '''
        try:
            # blast compares files and generates alignment results
            log.info(runBlastCompair(self._sPrimerFastaFile,
                     self._sAllTemplateFastFile, self._sBlastOutputFile))
            # Analyze the alignment results and get the recommended primer information
            res = pd.read_csv(self._sBlastOutputFile, sep=',',
                              names=['primerID', 'templateID', 'identity', 'length', 'primerStart', 'primerEnd',
                                     'templateStart', 'templateEnd'], index_col=['primerID'])

            lOut = list()
            for index, (sLeft, sRight) in enumerate(self._lPrimerPair):
                leftTm, leftGC, rightTm, rightGC = calcTM(
                    sLeft), calcGC(sLeft), calcTM(sRight), calcGC(sRight)
                leftMatch, rightMatch = list(), list()
                leftMatchLine, leftMaxMatch, rightMatchLine, rightMaxMatch = 0, iWordSize-1, 0, iWordSize-1
                leftID = 'primer_{}'.format(index * 2)
                rightID = 'primer_{}'.format(index * 2 + 1)
                if leftID in res.index:
                    leftMatch = res.loc[leftID].values.tolist()
                    leftMatchLine = len(res.loc[leftID])
                    leftMaxMatch = res.loc[leftID].length.max()
                if rightID in res.index:
                    rightMatch = res.loc[rightID].values.tolist()
                    rightMatchLine = len(res.loc[rightID])
                    rightMaxMatch = res.loc[rightID].length.max()
                lOut.append(
                    [sLeft, leftGC, leftTm, sRight, rightGC, rightTm, leftMatch, leftMatchLine, leftMaxMatch, rightMatch, rightMatchLine, rightMaxMatch])

            self._dfOutPut = pd.DataFrame(lOut, columns=(
                ['left', 'leftGC', 'leftTM', 'right', 'rightGC', 'rightTM', 'leftMatch', 'leftMatchLine', 'leftMaxMatch', 'rightMatch', 'rightMatchLine',
                 'rightMaxMatch']))
            # Calculate the total maximum number of matched bases and the total number of matched templates for the left and right primers
            self._dfOutPut['MaxMatchSum'] = self._dfOutPut.leftMaxMatch + \
                self._dfOutPut.rightMaxMatch
            self._dfOutPut['MatchLineSum'] = self._dfOutPut.leftMatchLine + \
                self._dfOutPut.rightMatchLine
            self._dfOutPut = self._dfOutPut.sort_values(
                ['MaxMatchSum', 'MatchLineSum']).reset_index(drop=True)
            return self._dfOutPut

        except Exception as e:
            log.error("blastCheck failed : {}".format(e) + str(format_exc()))
            raise


##############################################################################
# general function
##############################################################################

def getPrimerPair(iRandomSeqNum=50, iRandomSeqLen=200, iPrimerLen=20, iPrimer3PairNum=5):
    '''
    Design primer pairs using primer3
    @param iRandomSeqNum: Number of random sequences used to generate primers
    @param iPrimer3Num: Number of primer pairs generated using primer3 per random sequence
    @return: List of primer pairs
    '''
    lPrimerPair = list()
    tm = datetime.now()
    lTemplate = getRandomSeq(iRandomSeqNum, iRandomSeqLen)
    primer3Setting.setMaxNumber(iPrimer3PairNum)
    primer3Setting.setPrimerLen(iPrimerLen)

    for index,sTemplateSeq in enumerate(lTemplate):
        log.debug("{}, {}, {}\n".format(iPrimerLen,index,sTemplateSeq))
        primer3Setting.setTemplateSeq(sTemplateSeq)
        # primer3 output result
        result = primer3.bindings.designPrimers(
            primer3Setting.seq_args, primer3Setting.global_args)
        # The number of primer pairs generated by primer3
        iPairNum = result['PRIMER_PAIR_NUM_RETURNED']

        if iPairNum < 1:
            log.warning("no results")
        else:
            for index in range(iPairNum):
                sLeft = result['PRIMER_LEFT_{}_SEQUENCE'.format(index)]
                sRight = result['PRIMER_RIGHT_{}_SEQUENCE'.format(index)]
                if check3End(sLeft) and check3End(sRight):
                    lPrimerPair.append((sLeft, sRight))
                    log.debug("{},{}".format(sLeft, sRight))
    log.debug((datetime.now()-tm).total_seconds())
    return lPrimerPair


def check3End(sSeq, iLen=8):
    '''
    Check 3' to avoid repeats of more than three consecutive bases
    The rule is: check iLen bases at the end of 3', if there are 4 consecutive repeating bases, return False, otherwise return True
    '''
    if len(sSeq) < iLen or iLen < 4:
        return True

    def f(a): return a[0] == a[1] == a[2] == a[3]
    for i in range(iLen-4):
        if f(sSeq[-iLen+i: -iLen+i+4]):
            return False
    if f(sSeq[-4:]):
        return False
    return True


# get random dna/rna sequence
def getRandomSeq(iNumberOfSeq, iRandLen):
    '''
    get random sequence
    @param iNumberOfSeq: number of sequences
    @param iRandLen: sequence length
    @return list of random sequences
    '''
    lList = list()
    for iNum in range(iNumberOfSeq):
        each_Seq = ""
        for k in range(iRandLen - 1):
            iBaseNum = random.randint(0, 3)
            each_Seq += dNumBase[str(iBaseNum)]
        lList.append(each_Seq)
    return lList


def getAllTemplateFasta(lFile, sOutputFile):
    """
    Merge fasta files
    @param lFile: list of file paths
    @param sOutputFile: Merged fasta file
    @return:
    """
    lAll = list()
    allIndex = 0
    lChar = ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G']
    for sFile in lFile:
        with open(sFile, 'r') as f:
            lLine = f.readlines()
        index = 1
        for line in lLine:
            if line[0] in lChar:
                lAll.append('>join{}\n'.format(index))
                lAll.append(line)
                index+=1
    with open(sOutputFile, 'w') as f:
        f.writelines(lAll)


def getOther(sSeq):
    '''
    Get the reverse complement
    '''
    sComp = ''.join([dMap[x] for x in sSeq])  # complementary
    sRevComp = sComp[::-1]  # reverse complement
    return sRevComp


def getQuaStr(sFile):
    '''
    Get quaternary list from fasta file
    @param sFile: file path
    @return: quaternary string
    '''
    lNum = list()
    with open(sFile, 'r') as f:
        lLine = f.readlines()
    lChar = ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G']
    for line in lLine:
        if line[0] in lChar:
            sSeq = line.replace('\n', '')
            lNum.append(seqToQua(sSeq))
    return lNum


def quaToSeq(num, iLimitLen=0):
    '''
    quaternary transbase sequence
    @param num: quaternary number or string
    @param iLimitLen: base length（Used to supplement the base corresponding to the high digit zero）
    @return: base string
    '''
    # A:0 C:1 G:2 T:3
    num = str(num)
    if len(num) < iLimitLen:
        num = '0'*(iLimitLen-len(num)) + num
    res = ''.join([dNumBase[i] for i in num])
    return res


def seqToQua(sSeq):
    '''
    base sequence to quaternary string
    '''
    res = ''.join([dBaseNum[i] for i in sSeq])
    return res


def getPrimerPairList(lPrimerPair, setTemp6):
    '''
    Get a list of primer pairs
    '''
    lPrimer = list()
    for sLeft, sRight in lPrimerPair:
        if int(seqToQua(sLeft[-6:]), 4) in setTemp6 or int(seqToQua(sRight[-6:]), 4) in setTemp6:
            continue
        lPrimer.append((sLeft, sRight))
    return lPrimer


def runBlastDB(sDB):
    '''
    Run the balst command to generate the template database
    :param sDB: Template fasta file
    :return:
    '''
    cmd = "makeblastdb -in {}  -parse_seqids -hash_index -dbtype nucl".format(
        sDB)
    out = popen(cmd, 'r')
    return str(out.read())

# cmd:blastn -task blastn -query primer.fasta -db blast/allTemplate.fasta -out blast_out.txt -word_size 11  -outfmt "10 delim=, qseqid sseqid pident nident qstart qend sstart send"


def runBlastCompair(sPrimer, sDB, sOutPutFile):
    '''
    Run the blast command to compare primer.fasta with the database
    :param sPrimer: Primer fasta file
    :param sDB: Template fasta file (need to have generated library)
    :param sOutPutFile: output
    :return:
    '''
    cmd = 'blastn -task blastn -query {} -db {} -out {} -word_size {}' \
          ' -outfmt "10 delim=, qseqid sseqid pident nident qstart qend sstart send"'.format(
              sPrimer, sDB, sOutPutFile, iWordSize)
    out = popen(cmd, 'r')
    return str(log.info(out.read()))


def calcGC(listSeq):
    if type(listSeq) == type(""):
        listSeq = list(listSeq)
    gc = (listSeq.count('G') + listSeq.count('C')) / len(listSeq)
    return gc


def calcTM(listSeq):
    if type(listSeq) == type(""):
        listSeq = list(listSeq)
    if len(listSeq) < 25:
        tm = 4 * (listSeq.count('G') + listSeq.count('C')) + \
            2 * (listSeq.count('A') + listSeq.count('T'))
    else:
        tm = primer3.calcTm(''.join(listSeq))
    return tm

def mergeTemPrimer(lPrimer, lTemFile):
    for tPrimer,sFilePath in zip(lPrimer, lTemFile):
        sLeftAdd = tPrimer[0]
        sRightAdd = getOther(tPrimer[1])
        print(sLeftAdd, sRightAdd, sFilePath)
        lLine = list()
        with open(sFilePath, 'r') as file:
            for line in file.readlines():
                line = line.strip('\n')
                # add encode parameter
                if line[0]==';':
                    lLine.append(line+'\n')
                elif line[0] == ">":
                    lLine.append(line + '\n')
                elif line[0] in ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G']:
                    lLine.append(sLeftAdd+line+sRightAdd + '\n')
        with open(sFilePath+'.merge', 'w') as f:
            f.writelines(lLine)
