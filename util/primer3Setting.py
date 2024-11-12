#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 16:39:02 2021

@author: cathy
@brief: primer3 setting
"""

iPrimerNumer = 100 #Maximum number of returned primer pairs
iPrimerLen = 20

def setTemplateSeq(sSeq):
    '''
    Set the nucleotide sequence for generating primers
    '''
    global seq_args
    seq_args['SEQUENCE_TEMPLATE'] = sSeq
    

def setMaxNumber(iNumber = 50):
    '''
    Change the optimal number of generated primer pairs
    '''
    global iPrimerNumer,global_args
    iPrimerNumber = iNumber
    global_args['PRIMER_NUM_RETURN'] = iPrimerNumber
    
def setPrimerLen(iLen = 20):
    '''
    Set optimal primer length
    '''
    global iPrimerLen,global_args
    iPrimerLen = iLen
    global_args['PRIMER_OPT_SIZE'] = iPrimerLen

# ************************************************************
# "Sequence" input tags start with SEQUENCE_... and describe a particular input sequence to Primer3. 
# They are reset after every Boulder record. 
# Errors in "Sequence" input tags invalidate the current record, but Primer3 will continue to process additional records.
seq_args = {
'SEQUENCE_ID' : 'example',
'SEQUENCE_TEMPLATE': ""
 }


# "Global" input tags start with PRIMER_... and describe the general parameters that Primer3 should use in its searches. 
# The values of these tags persist between input Boulder records until or unless they are explicitly reset. 
# Errors in "Global" input tags are fatal because they invalidate the basic conditions under which primers are being picked.

# Because the laboratory detection step using internal oligos is independent of the PCR amplification procedure, 
# internal oligo tags have defaults that are independent of the parameters that govern the selection of PCR primers. 
# For example, the melting temperature of an oligo used for hybridization might be considerably lower than that used as a PCR primer. 
global_args= {
     ### Base sequence setting: ###
    'PRIMER_TASK' : 'generic',
    'PRIMER_LOWERCASE_MASKING' : 0,
    'PRIMER_PICK_LEFT_PRIMER' : 1,

    'PRIMER_PICK_RIGHT_PRIMER' : 1,
    'PRIMER_NUM_RETURN' : iPrimerNumer,
    'PRIMER_OPT_SIZE' : iPrimerLen,
    'PRIMER_MIN_SIZE' : 18,
    'PRIMER_MAX_SIZE' : 24,
    'PRIMER_MAX_NS_ACCEPTED' : 0,
    'PRIMER_MAX_POLY_X' : 5,
    
    ###  GC  ###
    'PRIMER_OPT_GC_PERCENT' : 50,
    'PRIMER_MIN_GC' : 40,
    'PRIMER_MAX_GC' : 60,
    'PRIMER_WT_GC_PERCENT_LT' : 0.0,
    'PRIMER_WT_GC_PERCENT_GT' : 0.0,
    
    ###  TM  ###
    'PRIMER_TM_FORMULA' : 1,
    'PRIMER_MIN_TM' : 58,
    'PRIMER_MAX_TM' : 70,
    'PRIMER_OPT_TM' : 60,
    'PRIMER_PAIR_MAX_DIFF_TM' : 4.0,
    'RIMER_WT_TM_LT' : 0.0,
    'PRIMER_WT_TM_GT' : 0.0,
    'PRIMER_PAIR_WT_DIFF_TM' : 0.0,

    'PRIMER_EXPLAIN_FLAG' : 1,
    

    'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT' : 1,

    'PRIMER_MAX_SELF_ANY' : 8.00,
    'PRIMER_WT_SELF_ANY' : 0.0,
    'PRIMER_MAX_SELF_ANY_TH' : 45.00,
    'PRIMER_WT_SELF_ANY_TH' : 123.2,

    'PRIMER_MAX_SELF_END' : 3.00,
    'PRIMER_WT_SELF_END' : 0.0,
    'PRIMER_MAX_SELF_END_TH' : 35.00,
    'PRIMER_WT_SELF_END_TH' : 302.4,

    'PRIMER_PAIR_MAX_COMPL_ANY' : 8.00,
    'PRIMER_PAIR_WT_COMPL_ANY' : 0.0,
    'PRIMER_PAIR_MAX_COMPL_ANY_TH' : 45.00,
    'PRIMER_PAIR_WT_COMPL_ANY_TH' : 123.2,

    'PRIMER_PAIR_MAX_COMPL_END' : 3.00,
    'PRIMER_PAIR_WT_COMPL_END' : 0.0,
    'PRIMER_PAIR_MAX_COMPL_END_TH' : 35.00,
    'PRIMER_PAIR_WT_COMPL_END_TH' : 302.4,

    'PRIMER_MAX_HAIRPIN_TH' : 24.00,
    'PRIMER_WT_HAIRPIN_TH' : 672,

    'PRIMER_MAX_END_STABILITY' : 9.0,
    'PRIMER_WT_END_STABILITY' : 1,

    'PRIMER_MUST_MATCH_THREE_PRIME' : 'nnnnb',

    'PRIMER_PAIR_WT_PR_PENALTY' : 1.0,
    }