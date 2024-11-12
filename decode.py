from util.DNAException import *
from collections import Counter
import single_edit_code as sec
from operator import itemgetter
import json
import argparse
import math
import sys
import os
import re

complementary_table = {"A":"T", "C":"G", "T":"A", "G":"C", "N":"N"}
base_xor = [[0, 1, 2, 3], [1, 0, 3, 2],
            [2, 3, 0, 1], [3, 2, 1, 0]]
file_types = {"txt":"ATCG", 'docx':'ATGC', 'log':'TACG', 'jpg':'TCGA', 
             'png':'ACTG', 'gif':'ACGT', 'bmp':'AGCT', 'tif':'AGTC',
             'mp3':'TAGC', 'av':'TCAG', 'aac':'TGAC', 'flac':'TGCA',
             'ogg':'CATG', 'xlsx':'CAGT', 'csv':'CTAG', 'zip':'CTGA',
             'rar':'CGAT', '7z':'CGTA', 'gz':'GATC', 'exe':'GACT',
             'app':'GTCA', 'html':'GTAC', 'css':'GCTA', 'js':'GCAT',
             'db':'GATG', 'sql':'GTAG', 'c':'CATC', 'py':'CTAC',
             'java':'ACGA', 'js':'AGCA', 'ppt':'TCGT', 'pdf':'TGCT',
             'md':'ACAC', 'epub':'CACA','mp4':'GCAC'}
units_dict = {
    'B':'TATA', 'KB':'ATAT', 'MB':'AGAG', 'GB':'GAGA', 'TB':'TCTC'
}

def get_opts():
    opts = argparse.ArgumentParser()
    opts.add_argument('-i','--input',required=True,help='input file')
    opts.add_argument('-o','--output',required=True,help='output file')
    opts.add_argument('-d', '--dir', required=True, help='directory of decode file')
    opts.add_argument('-c','--config',required=True,help='config file of primer')
    
    return opts.parse_args()
    

class BaseDecoder(object):
    def __init__(self, curr_dir, input_data, config, decode_file_dir) -> None:
        try:
            self.code = sec.SingleEditCode(int(48 * math.log(140, 4)))  # 可改
            self.curr_dir = curr_dir
            self.input = input_data
            self.primer1 = config['primer1']
            self.primer2 = config['primer2']
            self.decode_file_dir = decode_file_dir
            self.slice = self._read_in_slice()
            self.slice_map, self.slice_inv_map = self._read_in_slice_map()
            self.slice_inv_index = self._build_inv_index()
            self.record = self._read_in_record()
            self.forbidden_list = self._read_in_forbidden_list()
            self.ecc = self._read_in_ecc_data()
            self.match_seq_count_dict = self._read_in_base_data()
        except ConfigFileException as ce:
            print(ce)
            exit(0)

    def _read_in_slice(self, slice_file='slice.txt'):
        # 示例：12300131\n
        try:
            sf = open(self.curr_dir + os.sep + 'util' + os.sep + slice_file, 'r')
            slice = sf.readlines()
            sf.close()
            slice = [s[:-1] for s in slice]
        except:
            raise ConfigFileException('slice file error\n')
        else:
            return slice

    def _read_in_slice_map(self, slice_map_file='slice_map.txt'):
        # 示例：31995 37\n
        try:
            mf = open(self.curr_dir + os.sep + 'util' + os.sep + slice_map_file, 'r')
            slice_map = mf.readlines()
            mf.close()
            slice_map = [int(s.split(' ')[0]) for s in slice_map]

            # 只有哈希表才为O(1)
            slice_inv_map = dict()
            for id, s in enumerate(slice_map):
                slice_inv_map[s] = id
        except:
            raise ConfigFileException('slice map file error\n')
        else:
            return slice_map, slice_inv_map

    def _build_inv_index(self):
        try:
            slice_inv_index = dict()
            for id, s in enumerate(self.slice):
                slice_inv_index[s] = id
        except:
            raise ConfigFileException('slice file error\n')
        else:
            return slice_inv_index

    def _read_in_record(self, record_file='record.txt'):
        # 示例，00100011-2,8,4-\n，地址-禁忌序列id,替换串长度,索引\n
        try:
            print(self.decode_file_dir)
            rf = open(self.decode_file_dir + os.sep + record_file, 'r')
            raw_record = rf.readlines()
            raw_record = [r.split('-')[:-1] for r in raw_record]

            record = dict()

            for row in raw_record:
                single_record = list()

                current_addr = 0
                for id, item in enumerate(row):
                    if id == 0:
                        current_addr = item
                        continue

                    single_record.append([int(i) for i in item.split(',')])

                record[current_addr] = single_record
        except:
            raise ConfigFileException('record file error\n')
        else:
            return record

    def _read_in_forbidden_list(self, forbidden_list_file='forbidden_list.txt'):
        # 示例：12300131\n
        try:
            ff = open(self.curr_dir + os.sep + 'util' + os.sep + forbidden_list_file, 'r')
            forbidden_list = ff.readlines()
            ff.close()
            forbidden_list = [f[:-1] for f in forbidden_list]
        except:
            raise ConfigFileException('forbidden list file error\n')
        else:
            return forbidden_list

    def base_quad(self, base):
        quads = ''
        base_dict = {'A':'0', 'T':'1', 'C':'2', 'G':'3'} 
        for i in range(len(base)):
            quads += base_dict[base[i]]
        return quads

    def ecc_decode(self,code_data, addr):
        n, N, l = self.ecc[addr][0], self.ecc[addr][1], self.ecc[addr][2]
        code_data = sec.QaryString(4, list(code_data))
        raw_data = self.code.decode(code_data, n, N, l, verbose=False)
        raw_data = ''.join([str(r.val[0]) for r in raw_data])
        return raw_data


    def check_seq(self)->dict:
        '''
        check sequences whether pass error correct, base to quad, quad to bit
        '''
        id_addr_info = dict()               # dict[id_number] = [addr_quad, info_seq_quad]
        pass_ecc_seq = dict()
        for k,v in self.match_seq_count_dict.items():
            r = self.base_quad(k)
            addr = r[:8]
            info = r[8:]
            try:
                id = self.slice_inv_map[self.slice_inv_index[addr]]
            except:
                continue
            if addr in self.record.keys():
                info = self.revert_forbidden(
                    info, self.record[addr])
            try:
                info = self.ecc_decode(info, addr)
            except:
                continue
            if id not in id_addr_info.keys():
                id_addr_info[id] = [[addr, info, v]]
                pass_ecc_seq[id] = [self.primer1 + k + self.primer2]
            else:
                id_addr_info[id].append([addr, info, v])
                pass_ecc_seq[id].append(self.primer1 + k + self.primer2)

        pass_decode_bit = dict()    # seq dict which pass quad to bit
        for id in range(len(id_addr_info)):
            for lst in id_addr_info[id]:
                if self.decode_frag(lst[1]):
                    revert_bits = self.decode_frag(lst[1])
                    if id not in pass_decode_bit.keys():
                        pass_decode_bit[id] = [[revert_bits,lst[2]]]
                    else:
                        pass_decode_bit[id].append([revert_bits,lst[2]])
                        
        return pass_decode_bit

    def decode_frag(self, base_frag):
        # to_bytes function
        base_len = 8
        binary_segment_len = 15
        try:
            # print('{:0>15b}'.format(43))
            revert_binary_frag = ['{:0>{}b}'.format(self.slice_inv_map[self.slice_inv_index[base_frag[i:i+base_len]]],
                                                    binary_segment_len) for i in range(0, len(base_frag), base_len)]
            revert_binary_frag = ''.join(revert_binary_frag)
            revert_binary_frag = [revert_binary_frag[i*base_len:(
                i+1)*base_len] for i in range(0, int(len(revert_binary_frag) / base_len))]

            revert_binary_frag = [int(f, 2).to_bytes(1, 'big')
                                for f in revert_binary_frag]
            revert_binary_frag = b''.join(revert_binary_frag)
        except:
            return False
        return revert_binary_frag


    def _read_in_base_data(self):
        # 示例：16碱基地址|128碱基\n
        seqs = [i.strip() for i in self.input]
        primer1_r = ''.join([complementary_table[i] for i in self.primer1[::-1]])
        primer2_r = ''.join([complementary_table[i] for i in self.primer2[::-1]])
        pattern_1 = re.compile("(\.*{})(\\S+){}\.*".format(self.primer1, self.primer2))
        pattern_2 = re.compile("(\.*{})(\\S+){}\.*".format(primer2_r,primer1_r))

        merged_seq_count = Counter(seqs)
        merged_seq_count = dict(merged_seq_count)
        match_seq_count_dict = dict()
        for seq in set(seqs):
            match_seq_1 = re.finditer(pattern_1, seq)
            match_seq_2 = re.finditer(pattern_2, seq)

            seq1_list = [[i.group(2), i.span(2)] for i in match_seq_1]
            seq2_list = [[i.group(2), i.span(2)] for i in match_seq_2]
            if seq1_list != []:
                try:
                    _seq1 = seq1_list[0][0]
                    match_seq_count_dict[_seq1] += merged_seq_count[seq]
                except:
                    _seq1 = seq1_list[0][0]
                    match_seq_count_dict[_seq1] = merged_seq_count[seq]
            elif seq2_list != []:
                try:
                    _seq2 = "".join([complementary_table[i] for i in seq2_list[0][0][::-1]])
                    match_seq_count_dict[_seq2] += merged_seq_count[seq]
                except:
                    _seq2 = "".join([complementary_table[i] for i in seq2_list[0][0][::-1]])
                    match_seq_count_dict[_seq2] = merged_seq_count[seq]
            else:
                continue
        
        return match_seq_count_dict

    def _read_in_ecc_data(self, ecc_file='ecc.txt'):
        try:
            ef = open(self.decode_file_dir + os.sep + ecc_file, 'r')
            ecc_list = ef.readlines()
            ef.close()

            ecc_dict = dict()
            for ecc in ecc_list:
                ecc = ecc[:-1].split('-')
                ecc_dict[ecc[0]] = [int(ecc[1]), int(ecc[2]), int(ecc[3])]
        except:
            raise ConfigFileException('ecc file error\n')
        else:
            return ecc_dict

    def _ecc_decode(self, code_data, addr):
        n, N, l = self.ecc[addr][0], self.ecc[addr][1], self.ecc[addr][2]
        code_data = sec.QaryString(4, list(code_data))
        raw_data = self.code.decode(code_data, n, N, l, verbose=False)
        raw_data = ''.join([str(r.val[0]) for r in raw_data])
        return raw_data

    def quad_base(self, aa):
        bases = ''
        quad_dict = {'0':'A', '1':'T', '2':'C', '3':'G'}
        for i in range(len(aa)):
            bases += quad_dict[aa[i]]
        return bases

    def base_quad(self, base):
        quads = ''
        base_dict = {'A':'0', 'T':'1', 'C':'2', 'G':'3'} 
        for i in range(len(base)):
            if base[i] not in base_dict.keys():
                continue
            else:
                quads += base_dict[base[i]]
        return quads

    def seq_xor(self, bases):
        result = ''
        for i in range(len(bases[0])):
            x = 0
            for b in bases:
                x = base_xor[int(b[i])][x]
            result = result + str(x)
        return result

    def revert_forbidden(self, base_data, record):
        record_id_forbidden_id, record_id_replace_len, record_id_replace_index = 0, 1, 2
        for r in record:
            forbidden_id, replace_len, replace_index = r[record_id_forbidden_id], r[
                record_id_replace_len], r[record_id_replace_index]
            base_data = base_data[:replace_index] + \
                self.forbidden_list[forbidden_id] + \
                base_data[replace_len+replace_index:]
        return base_data

    def redundancy(self, bases):
        if len(bases) == 0:
            return None

        ini = bytearray(b''.join([b'\x00'] * len(bases[0])))
        for b in bases:
            ini = bytearray(x ^ y for x, y in zip(ini, bytearray(b)))

        return bytes(ini)

    def decode_file_attr(self, sum_n):
        primer_addr = ['GCAGCGTGTGAGATTCATGG','GCTGTACGTCGCACCTCTAG']
        primer_file_info = ['GGCCGCTTTCGTCACATAAC','CGAGACCCGCAACTTGACTG']
        primer_file_xor = ['CGATCGTTGCACCTCTTGAC','TCTGTGACTGCGACCCAAGC']
        other_primer_lst = [primer_addr, primer_file_info, primer_file_xor]
        other_match_seq = dict()

        other_match_list = [[],[],[]]
        for seq in set(self.match_seq_count_dict.keys()):
            for n,(r1, r2) in enumerate(other_primer_lst):
                r2_r = ''.join([complementary_table[i] for i in r2[::-1]])
                r1_r = ''.join([complementary_table[i] for i in r1[::-1]])
                pattern_1 = re.compile("(\.*{})(\\S+){}\.*".format(r1,r2))
                pattern_2 = re.compile("(\.*{})(\\S+){}\.*".format(r2_r,r1_r))
                match_seq_1 = re.finditer(pattern_1, seq)
                match_seq_2 = re.finditer(pattern_2, seq)

                seq1_list = [[i.group(2), i.span(2)] for i in match_seq_1]
                seq2_list = [[i.group(2), i.span(2)] for i in match_seq_2]
                if seq1_list != []:
                    try:
                        _seq1 = seq1_list[0][0]
                        other_match_seq[_seq1] += self.match_seq_count_dict[seq]
                    except:
                        _seq1 = seq1_list[0][0]
                        other_match_seq[_seq1] = self.match_seq_count_dict[seq]
                    other_match_list[n].append(_seq1)
                    _seq1 = ''
                elif seq2_list != []:
                    try:
                        _seq2 = "".join([complementary_table[i] for i in seq2_list[0][0][::-1]])
                        other_match_seq[_seq2] += self.match_seq_count_dict[seq]
                    except:
                        _seq2 = "".join([complementary_table[i] for i in seq2_list[0][0][::-1]])
                        other_match_seq[_seq2] = self.match_seq_count_dict[seq]
                    other_match_list[n].append(_seq2)
                    _seq2 = ''
                else:
                    continue
        file_addr_pass = []
        file_info_pass = []
        file_xor_pass = []
        for seq in other_match_list[0]:
            r = self.base_quad(seq)
            addr = r[:8]
            add_lst = r[8:]
            if addr in self.record.keys():
                add_lst = self.revert_forbidden(
                    add_lst, self.record[addr])
            try:
                add_lst = self.ecc_decode(add_lst, addr)
            except:
                continue
            file_addr_pass.append([addr, add_lst, other_match_seq[seq]])

        for seq in other_match_list[1]:
            r = self.base_quad(seq)
            addr = r[:8]
            info_seq = r[8:]
            if addr in self.record.keys():
                info_seq = self.revert_forbidden(
                    info_seq, self.record[addr])
            try:
                info_seq = self.ecc_decode(info_seq, addr)
            except:
                continue
            file_info_pass.append([addr, info_seq, other_match_seq[seq]])

        for seq in other_match_list[2]:
            r = self.base_quad(seq)
            addr = r[:8]
            xor_info = r[8:]
            if addr in self.record.keys():
                xor_info = self.revert_forbidden(
                    xor_info, self.record[addr])
            try:
                xor_info = self.ecc_decode(xor_info, addr)
            except:
                continue
            file_xor_pass.append([addr, xor_info, other_match_seq[seq]])

        sorted_addr_count = sorted(file_addr_pass,key=itemgetter(2),reverse=True)
        sorted_info_count = sorted(file_info_pass,key=itemgetter(2),reverse=True)
        sorted_xor_count = sorted(file_xor_pass,key=itemgetter(2),reverse=True)
        com_addr_bit = sorted_addr_count[0][1]
        com_info_bit = sorted_info_count[0][1]
        com_xor_bit = sorted_xor_count[0][1]
        com_lst = [com_addr_bit, com_info_bit, com_xor_bit]
        
        if self.seq_xor([com_addr_bit[:-17], com_info_bit[:-17]]) != com_xor_bit:
            for i in range(3):
                for tmp_seq in other_match_list[i]:
                    measuer_lst = [com_addr_bit[:-17], com_info_bit[:-17], com_xor_bit]
                    tmp_list = measuer_lst.pop(i)
                    xor_bit = self.seq_xor(tmp_list)
                    sam_n = 0
                    for index in range(len(xor_bit)):
                        if xor_bit[index] == tmp_seq[index]:
                            sam_n += 1
                    if sam_n > len(xor_bit) * 0.5:
                        r1 = xor_bit + tmp_seq[:-17]
                        com_lst[i] = r1
        addr_lst_bit = com_lst[0]
        file_info_bit = com_lst[1]
        addr_lst = []
        for i in range(0,len(addr_lst_bit),8):
            addr_lst.append(
            str(self.slice_inv_map[self.slice_inv_index[addr_lst_bit[i:i+8]]]))

        # file attribution decode
        file_type = [k for k,v in file_types.items() if v == self.quad_base(file_info_bit[16:20])][0]
        date_y, date_m, date_d = file_info_bit[20:28], file_info_bit[28:36], file_info_bit[36:44]
        year = str(self.slice_inv_map[self.slice_inv_index[date_y]])
        m_d = str(self.slice_inv_map[self.slice_inv_index[date_m]])
        month = m_d[:-2]
        day = m_d[-2:]
        min_sec = str(self.slice_inv_map[self.slice_inv_index[date_d]])
        mine = min_sec[:-2]
        sec = min_sec[-2:]
        size_pre, size_suf, size_unit = file_info_bit[44:52],file_info_bit[52:60],file_info_bit[60:64]
        unit =  [k for k,v in units_dict.items() if v == self.quad_base(size_unit)]
        file_perm = file_info_bit[-20:]
        if file_perm == '32120231120011233232':
            perm = 'user'
        elif file_perm == '20030232330302033122':
            perm = 'admin'
        else:
            perm = 'other'
        date = year + '. ' + month + '. '+ day+'\t'+mine+':'+sec
        file_size = str(self.slice_inv_map[self.slice_inv_index[size_pre]]) \
                    + '.'+str(self.slice_inv_map[self.slice_inv_index[size_suf]]) \
                    + ' ' +unit[0]
        file_addr = ','.join(addr_lst)
        total_size = str(self.slice_inv_map[self.slice_inv_index[size_pre]]) \
                    + '.'+str(self.slice_inv_map[self.slice_inv_index[size_suf]])
        exp_num = 10 * list(units_dict.keys()).index(unit[0])
        total_size = round(float(total_size) * (2 ** exp_num))
        return file_type, date, file_size, perm, file_addr, total_size

    def decode(self)->str:
        result = b''
        pass_decode_seq = self.check_seq()
        all_bits_frag = []
        for i in range(len(pass_decode_seq)):
            tmp_count_dict = dict()
            seq_count = pass_decode_seq[i]
            for i_1, i_2 in seq_count:
                if i_1 not in tmp_count_dict.keys():
                    tmp_count_dict[i_1] = i_2
                tmp_count_dict[i_1] += i_2
            seq_count_lst = [[k,v] for k,v in tmp_count_dict.items()]
            sorted_count = sorted(seq_count_lst,key=itemgetter(1),reverse=True)
            s, n = sorted_count[0]
            all_bits_frag.append(s)

        tmp_list = []
        for i in range(len(all_bits_frag)):
            if (i + 1) % (5 + 1) != 0:
                tmp_list.append(all_bits_frag[i])
            else:
                if self.redundancy(tmp_list) == all_bits_frag[i]:
                    for j in tmp_list:
                        result += j
                    tmp_list = []
                    continue
                else:
                    print("redundancy: ",all_bits_frag[i])
                    for k in range(1,5):
                        for l in pass_decode_seq[i-k]:
                            if k == 1:
                                merge_list = list(tmp_list[:-k])+[l[0]]
                            else:
                                merge_list = list(tmp_list[:-k])+[l[0]]+list(tmp_list[-k-1:])
                            if self.redundancy(merge_list) == all_bits_frag[i]:
                                tmp_list = list(merge_list)
                            else:
                                redun = self.redundancy(merge_list)
                                sam_n = 0
                                for index in range(len(redun)):
                                    if redun[index] == all_bits_frag[i][index]:
                                        sam_n += 1
                                if sam_n > len(redun) *0.5:
                                    r1 = self.redundancy(list(tmp_list[:-k])+[redun]+list(tmp_list[-k:]))
                                    tmp_list = list(tmp_list[:-k])+[r1]+list(tmp_list[-k:])
                                else:
                                    continue
                    for j in tmp_list:
                        result += j
                    tmp_list = []

        for j in tmp_list:
            result += j
        return result


def main():
    opts = get_opts()
    input = opts.input
    output = opts.output
    dir = opts.dir
    with open(opts.config, 'r') as f:
        config = json.load(f)
    with open(input, 'r') as f:
        input_data = f.readlines()
         
    decode_file_dir = dir
    curr_dir = os.getcwd()
    workdir = os.getcwd() + os.sep + output
    try:
        os.mkdir(workdir)
    except:
        pass
    os.chdir(workdir)
    decoder = BaseDecoder(curr_dir, input_data, config, decode_file_dir)
    result = decoder.decode()
    sum_n = len(result)
    file_type, date, file_size, perm, file_addr, total_size = decoder.decode_file_attr(sum_n)
    attr_file = 'file_attribution.txt'
    
    with open(attr_file,'w') as f:
        f.write("file name:\t" + 'recover.' + file_type + '\n')
        f.write("file creation date:\t"+ date+'\n')
        f.write('file size:\t' + file_size+'\n' )
        f.write('file permission:\t' + perm + '\n')
        f.write('file address:\t'+file_addr)
    
    recover_file = 'recover.' + file_type
    with open(recover_file, 'wb') as f:
        f.write(result[:total_size])
    
if __name__ == '__main__':
    main()
