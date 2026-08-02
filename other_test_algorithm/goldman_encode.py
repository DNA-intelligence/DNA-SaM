# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 10:37:21 2021

@author: Wei QWiang
"""

import sys, os, argparse
from tqdm import tqdm
from time import time



# read_file_path = "../input_file/3390.txt"
# dna_path = "../goldman_test/3390.dna"


# read_file_path = "../input_file/big_fac.jpg"
# dna_path = "../goldman_test/big_fac_modify.dna"

## BGI list
goldman_list = ["22201", "00100", "11220", "00211", "20222", "00222", "02211", "222110",
                "22002", "02100", "22001", "222122", "12001", "02021", "10100", "02010",
                "20101", "12211", "12120", "11111", "21211", "21221", "20220", "00122",
                "20022", "12121", "21111", "00221", "00202", "222202", "222102", "00010",
                "02212", "10011", "22011", "02221", "21212", "21021", "11211", "10111",
                "12220", "22110", "22101", "11122", "22022", "01210", "00210", "02122",
                "10122", "01011", "11101", "01102", "22112", "12122", "11012", "222112",
                "02201", "02011", "20021", "222021", "00022", "222200", "222120", "21010",
                "00121", "02022", "20100", "10211", "21001", "21210", "10212", "222212",
                "20110", "20010", "21220", "21022", "21000", "01211", "10220", "12002",
                "12011", "11212", "21100", "12210", "20112", "22200", "22102", "21222",
                "21012", "12101", "10120", "01202", "10200", "02210", "222211", "11201",
                "00102", "01112", "22010", "00012", "22100", "20001", "20202", "02102",
                "20200", "20210", "20012", "11100", "02101", "11021", "00021", "02110",
                "12102", "01012", "10101", "10222", "10221", "10002", "01120", "00201",
                "10020", "222111", "222220", "02111", "222222", "00000", "10112", "22121",
                "02000", "10000", "20111", "00212", "22021", "21112", "11022", "01220",
                "11102", "20011", "22111", "10021", "12212", "11202", "10201", "02200",
                "02002", "11120", "20102", "11110", "11002", "22000", "21002", "21102",
                "222221", "11020", "20221", "01002", "11001", "00120", "02202", "10202",
                "10012", "22012", "20211", "21201", "00220", "11222", "21011", "10110",
                "20002", "20122", "22122", "20201", "10022", "21101", "12110", "12222",
                "00200", "21202", "10210", "10010", "02012", "12221", "12022", "02222",
                "01100", "02121", "01122", "00112", "01020", "222100", "01222", "21020",
                "01201", "00001", "12021", "12010", "20121", "21120", "00002", "222201",
                "00011", "01010", "12112", "11112", "02120", "11010", "01110", "01212",
                "20120", "12000", "12100", "11210", "11011", "21200", "12200", "01111",
                "01200", "12012", "10121", "10102", "222210", "00020", "01000", "20020",
                "11121", "10001", "02001", "01101", "222121", "21121", "02220", "01001",
                "222101", "01022", "20212", "00101", "222022", "01021", "00111", "11200",
                "12201", "11000", "02112", "01221", "00110", "11221", "01121", "12111",
                "12020", "02020", "22020", "20000", "21110", "22120", "12202", "21122", "222020"]

# rotate_codes_dic = {'A': ['C', 'G', 'T'], 'C': ['G', 'T', 'A'], 'G': ['T', 'A', 'C'], 'T': ['A', 'C', 'G']}
# last_nucleotide = "A"


def linuxCommand():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help = "input path")
    parser.add_argument("-o", help = "output path")
    parser.add_argument("-l", help = "ideal length", type=int)
    args = parser.parse_args()
    
    return args.i, args.o, args.l

def transSystem(number: int, system: int):
    trans_num_list = []
    while True:
        quot = number//system
        remiander = number%system
        trans_num_list.append(str(remiander))
        if quot == 0:
            break
        else:
            number = quot
    trans_num = "".join(trans_num_list[::-1])
    return trans_num


def transHuffman(bin_num):
    # 5位无法全部表示，这样处理可唯一 详见BGI list
    decimal_num = int(bin_num, 2)
    if decimal_num > 235:
        decimal_num += 472

    # transform to ternary number
    ternary_num = transSystem(decimal_num, 3)

    if len(ternary_num) < 5:
        ternary_num = "0"*(5-len(ternary_num)) + ternary_num
    
    return ternary_num



def readAsBin(read_file_path):
    f_path = read_file_path
    bin_str = ""
    tmp_list = []
    with open(f_path, "rb") as f:
        for _line in f:
            for _bytes in _line:
                tmp_list.append(_bytes)
                _bin_str = bin(_bytes)[2:]
                _bin_str = "0"*(8-len(_bin_str)) + _bin_str
                bin_str += _bin_str
    # print(tmp_list)
    return bin_str

def transTernaryNum(bin_str):
    ternary_str = ""
    for i in range(0, len(bin_str), 8):
        _bin_str = bin_str[i: i+8]
        _terary_num = transHuffman(_bin_str)

        ternary_str += _terary_num

    return ternary_str

def encodeNt(ternary_str):
    last_nt="A"
    rotate_codes_dic = {'A': ['C', 'G', 'T'], 'C': ['G', 'T', 'A'], 'G': ['T', 'A', 'C'], 'T': ['A', 'C', 'G']}

    nt_seq = ""
    for _terary_num in ternary_str:
        _nt = rotate_codes_dic[last_nt][int(_terary_num)]
        nt_seq += _nt
        last_nt = _nt

    return nt_seq


def segment(ternary_str, ideal_len = 100):
    total_len = len(ternary_str)

    n = 0
    while True:
        seg_len = int(ideal_len // 4 * 4-4 * n)
        if total_len%seg_len == 0:
            seg_num = int(total_len//(seg_len/4)-3)
        else:
            seg_num = int(total_len//(seg_len/4)-2)

        idnex_len = len(transSystem(seg_num, 3))
        if idnex_len > (ideal_len-seg_len):
            n += 1
        else:
            break
    

    ternary_seg_list = []
    step = int(seg_len/4)
    _index = 0
    add_len = 0
    for i in range(0, len(ternary_str), step):
        tag = False
        seg_seq = ternary_str[i: i+seg_len]

        # no remiander nt
        if len(seg_seq) == 3*step:
            return ternary_seg_list, idnex_len, seg_len-len(seg_seq)
        # remainder nt
        elif len(seg_seq) < seg_len:
            tag = True
            add_len = seg_len-len(seg_seq)
            seg_seq += "0"*add_len
            # for _i in range(seg_len-len(seg_seq)):
            #     seg_seq += list({"A", "T", "C", "G"}-set(seg_seq[-1]))[0]

        _index_num = transSystem(_index, 3)
        _index_num = "0"*(idnex_len-len(_index_num))+_index_num
        ternary_seg = _index_num + seg_seq
        ternary_seg_list.append(ternary_seg)

        _index += 1
        if tag == True:
            return ternary_seg_list, idnex_len, add_len


def outputResult(output_path, nt_seq_list):
    with open(output_path, "w") as f:
        for nt_seq in nt_seq_list:
            _output_line = nt_seq + "\n"
            f.write(_output_line)


def goldmanEncode(bin_str, ideal_len = 100):
    ternary_str = transTernaryNum(bin_str)
    # print("ternary_str", len(ternary_str)) #bkp
    ternary_seg_list, idnex_len, add_len = segment(ternary_str, ideal_len=ideal_len)
    nt_seq_list = []
    for ternary_seg in ternary_seg_list:
        nt_seq = encodeNt(ternary_seg)
        nt_seq_list.append(nt_seq)

    return nt_seq_list, idnex_len, add_len, ternary_seg_list, ternary_str


def goldmanMain(read_file_path, output_path, ideal_len=100):
    bin_str = readAsBin(read_file_path)
    nt_seq_list, idnex_len, add_len, ternary_seg_list, ternary_str = goldmanEncode(bin_str, ideal_len=ideal_len)
        
    outputResult(output_path, nt_seq_list)

    return nt_seq_list, idnex_len, add_len



if __name__ == "__main__":
    read_file_path, output_path, ideal_len = linuxCommand()
    t1 = time()
    code_r = goldmanMain(read_file_path, output_path, ideal_len)
    print('[time use]: ', time()-t1)