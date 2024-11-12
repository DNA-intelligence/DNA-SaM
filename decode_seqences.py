from collections import Counter
from util.DNAException import *
import single_edit_code as sec
from operator import itemgetter
import os
import argparse
import json
import re
import time
import math

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
    group = argparse.ArgumentParser(description="Cluster and Filter sequences from PE sequencing reads")
    group.add_argument("-f", "--filter", help="value of proportion of low quality bases in the sequence", default=0.2, type=float)
    group.add_argument("-c", "--config", help="config path", required=True)
    group.add_argument("-o", "-output", help="output file", required=True)
    return group

class Decoder():
    def __init__(self):
        self.code = sec.SingleEditCode(int(48 * math.log(140, 4)))
        self.slice = self.read_in_slice()
        self.slice_map, self.slice_inv_map = self.read_in_slice_map()
        self.slice_inv_index = self.build_inv_index()
        self.record = self.read_in_record()
        self.forbidden_list = self.read_in_forbidden_list()
        self.ecc = self.read_in_ecc_data()   

    def read_in_forbidden_list(self,forbidden_list_file='forbidden_list.txt'):
        # 示例：12300131\n
        try:
            ff = open(forbidden_list_file, 'r')
            forbidden_list = ff.readlines()
            ff.close()
            forbidden_list = [f[:-1] for f in forbidden_list]
        except:
            raise ConfigFileException('forbidden list file error\n')
        else:
            return forbidden_list

    def read_in_ecc_data(self,ecc_file='ecc.txt'):
        try:
            ef = open(ecc_file, 'r')
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

    def ecc_decode(self,code_data, addr):
        n, N, l = self.ecc[addr][0], self.ecc[addr][1], self.ecc[addr][2]
        code_data = sec.QaryString(4, list(code_data))
        raw_data = self.code.decode(code_data, n, N, l, verbose=False)
        raw_data = ''.join([str(r.val[0]) for r in raw_data])
        return raw_data

    def revert_forbidden(self,base_data, record):
        record_id_forbidden_id, record_id_replace_len, record_id_replace_index = 0, 1, 2
        for r in record:
            forbidden_id, replace_len, replace_index = r[record_id_forbidden_id], r[
                record_id_replace_len], r[record_id_replace_index]
            base_data = base_data[:replace_index] + \
                self.forbidden_list[forbidden_id] + \
                base_data[replace_len+replace_index:]
        return base_data
    
    def read_in_record(self,record_file='record.txt'):
        # 示例，00100011-2,8,4-\n，地址-禁忌序列id,替换串长度,索引\n
        try:
            rf = open(record_file, 'r')
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
        
    def read_in_slice(self,slice_file='slice.txt'):
        # 示例：12300131\n
        try:
            sf = open(slice_file, 'r')
            slice = sf.readlines()
            sf.close()
            slice = [s[:-1] for s in slice]
        except:
            raise ConfigFileException('slice file error\n')
        else:
            return slice

    def read_in_slice_map(self,slice_map_file='slice_map.txt'):
        # 示例：31995 37\n
        try:
            mf = open(slice_map_file, 'r')
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

    def build_inv_index(self):
        try:
            slice_inv_index = dict()
            for id, s in enumerate(self.slice):
                slice_inv_index[s] = id
        except:
            raise ConfigFileException('slice file error\n')
        else:
            return slice_inv_index

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

    def check_seq(self, count_dict: dict, out_dir: str)->dict:
        '''
        check sequences whether pass error correct, base to quad, quad to bit
        '''
        match_seq_count_dict = count_dict   # match seq and count 
        id_addr_info = dict()               # dict[id_number] = [addr_quad, info_seq_quad]
        pass_ecc_seq = dict()
        pass_ecc_file = out_dir + os.sep + 'pass_ecc_file.txt'
        ecc_out = open(pass_ecc_file, 'w')
        for k,v in match_seq_count_dict.items():
            r = base_quad(k)
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
                pass_ecc_seq[id] = [primer1 + k + primer2]
            else:
                id_addr_info[id].append([addr, info, v])
                ecc_out.write(str(id) + '\t' + info + '\n')
                pass_ecc_seq[id].append(primer1 + k + primer2)
        ecc_out.close()
        pass_decode_bit = dict()    # seq dict which pass quad to bit
        for id in range(len(id_addr_info)):
            for lst in id_addr_info[id]:
                if self.decode_frag(lst[1]):
                    revert_bits = self.decode_frag(lst[1])
                    if id not in pass_decode_bit.keys():
                        pass_decode_bit[id] = [[revert_bits,lst[2]]]
                    else:
                        pass_decode_bit[id].append([revert_bits,lst[2]])
                        
        return pass_decode_bit, pass_ecc_seq
                
    def decode(self, pass_seq: dict)->str:
        result = b''
        pass_decode_seq = pass_seq
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
                if redundancy(tmp_list) == all_bits_frag[i]:
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
                            if redundancy(merge_list) == all_bits_frag[i]:
                                tmp_list = list(merge_list)
                            else:
                                redun = redundancy(merge_list)
                                sam_n = 0
                                for index in range(len(redun)):
                                    if redun[index] == all_bits_frag[i][index]:
                                        sam_n += 1
                                if sam_n > len(redun) *0.5:
                                    r1 = redundancy(list(tmp_list[:-k])+[redun]+list(tmp_list[-k:]))
                                    tmp_list = list(tmp_list[:-k])+[r1]+list(tmp_list[-k:])
                                else:
                                    continue
                    for j in tmp_list:
                        result += j
                    tmp_list = []

        for j in tmp_list:
            result += j
        return result
    
    def decode_file_attr(self, match_seq_count_dict,sum_n):
        primer_addr = ['GCAGCGTGTGAGATTCATGG','GCTGTACGTCGCACCTCTAG']
        primer_file_info = ['GGCCGCTTTCGTCACATAAC','CGAGACCCGCAACTTGACTG']
        primer_file_xor = ['CGATCGTTGCACCTCTTGAC','TCTGTGACTGCGACCCAAGC']
        other_primer_lst = [primer_addr, primer_file_info, primer_file_xor]
        other_match_seq = dict()

        other_match_list = [[],[],[]]
        for seq in set(match_seq_count_dict.keys()):
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
                        other_match_seq[_seq1] += match_seq_count_dict[seq]
                    except:
                        _seq1 = seq1_list[0][0]
                        other_match_seq[_seq1] = match_seq_count_dict[seq]
                    other_match_list[n].append(_seq1)
                    _seq1 = ''
                elif seq2_list != []:
                    try:
                        _seq2 = "".join([complementary_table[i] for i in seq2_list[0][0][::-1]])
                        other_match_seq[_seq2] += match_seq_count_dict[seq]
                    except:
                        _seq2 = "".join([complementary_table[i] for i in seq2_list[0][0][::-1]])
                        other_match_seq[_seq2] = match_seq_count_dict[seq]
                    other_match_list[n].append(_seq2)
                    _seq2 = ''
                else:
                    continue
        file_addr_pass = []
        file_info_pass = []
        file_xor_pass = []
        for seq in other_match_list[0]:
            r = base_quad(seq)
            if seq[:8] != 'CTAGAATC':
                continue
            else:
                addr = '21330020'
            if addr in self.record.keys():
                add_lst = self.revert_forbidden(
                    r, self.record[addr])
            try:
                info = self.ecc_decode(r, addr)
            except:
                continue
            addr = info[:8]
            add_lst = info[8:]
            file_addr_pass.append([addr, add_lst, other_match_seq[seq]])

        for seq in other_match_list[1]:
            r = base_quad(seq)
            addr = r[:8]
            
            try:
                id = self.slice_inv_map[self.slice_inv_index[addr]]
                if id != sum_n+1:
                    continue
            except:
                continue
            if addr in self.record.keys():
                info = self.revert_forbidden(
                    r, self.record[addr])
            try:
                info = self.ecc_decode(r, addr)
            except:
                continue
            info = info[8:]
            file_info_pass.append([addr, info, other_match_seq[seq]])

        for seq in other_match_list[2]:
            r = base_quad(seq)
            addr = self.slice[self.slice_map[sum_n+2]]
            info = r[8:]
            if addr in self.record.keys():
                info = self.revert_forbidden(
                    info, self.record[addr])
            try:
                info = self.ecc_decode(r, addr)
            except:
                continue
            info = info[8:]
            file_xor_pass.append([addr, info, other_match_seq[seq]])

        sorted_addr_count = sorted(file_addr_pass,key=itemgetter(2),reverse=True)
        sorted_info_count = sorted(file_info_pass,key=itemgetter(2),reverse=True)
        sorted_xor_count = sorted(file_xor_pass,key=itemgetter(2),reverse=True)
        com_addr_bit = sorted_addr_count[0][1]
        com_info_bit = sorted_info_count[0][1]
        com_xor_bit = sorted_xor_count[0][1]
        com_lst = [com_addr_bit, com_info_bit, com_xor_bit]
        
        if seq_xor([com_addr_bit[:-17], com_info_bit[:-17]]) != com_xor_bit:
            for i in range(3):
                for tmp_seq in other_match_list[i]:
                    measuer_lst = [com_addr_bit[:-17], com_info_bit[:-17], com_xor_bit]
                    tmp_list = measuer_lst.pop(i)
                    xor_bit = seq_xor(tmp_list)
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
        for i in range(0,len(addr_lst_bit)-1,8):
            addr_lst.append(
            str(self.slice_inv_map[self.slice_inv_index[addr_lst_bit[i:i+8]]]))

        # file attribution decode
        file_type = [k for k,v in file_types.items() if v == quad_base(file_info_bit[16:20])][0]
        date_y, date_m, date_d = file_info_bit[20:28], file_info_bit[28:36], file_info_bit[36:44]
        year = str(self.slice_inv_map[self.slice_inv_index[date_y]])
        m_d = str(self.slice_inv_map[self.slice_inv_index[date_m]])
        month = m_d[:-2]
        day = m_d[-2:]
        min_sec = str(self.slice_inv_map[self.slice_inv_index[date_d]])
        mine = min_sec[:-2]
        sec = min_sec[-2:]
        size_pre, size_suf, size_unit = file_info_bit[44:52],file_info_bit[52:60],file_info_bit[60:64]
        unit =  [k for k,v in units_dict.items() if v == quad_base(size_unit)]
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
        return file_type, date, file_size, perm, file_addr

def quad_base(aa):
    bases = ''
    quad_dict = {'0':'A', '1':'T', '2':'C', '3':'G'}
    for i in range(len(aa)):
        bases += quad_dict[aa[i]]
    return bases

def base_quad(base):
    quads = ''
    base_dict = {'A':'0', 'T':'1', 'C':'2', 'G':'3'} 
    for i in range(len(base)):
        if base[i] not in base_dict.keys():
            continue
        else:
            quads += base_dict[base[i]]
    return quads
    
def getTime():
    """
    get current time
    """
    cur_time = time.asctime(time.localtime(time.time()))
    return cur_time

def DictToFile(_dict, file):
    """
    save dict to file
    """
    f = open(file,'w')
    for key,value in _dict.items():
        f.write(key + "\t" + str(value) + "\n")

def seq_xor(bases):
    result = ''
    for i in range(len(bases[0])):
        x = 0
        for b in bases:
            x = base_xor[int(b[i])][x]
        result = result + str(x)
    return result

def readToDict(file):
    """
    read file to dict
    """
    _dict = dict()
    f = open(file)
    for i in f:
        i = i.strip().split("\t")
        key, value = i[0], i[1]
        _dict[key] = value
    
    return _dict

def redundancy(bases):
    if len(bases) == 0:
        return None

    ini = bytearray(b''.join([b'\x00'] * len(bases[0])))
    for b in bases:
        ini = bytearray(x ^ y for x, y in zip(ini, bytearray(b)))

    return bytes(ini)

def AssembleReads(pe1, pe2, result_dir):
    """ 
    assemble pair-end reads with pear software
    """
    # run pear command
    command = "pear -f {} -r {} -o {} -y 4G -j 24\n".format(pe1, pe2, result_dir+os.sep + 'reads')

    log.write("step1: run Reads assemble\t{}\n".format(getTime()))
    
    try:
        os.system(command)
    except:
        raise("pear fail run!\n")
    log.write("step1: done!\t{}\n".format(getTime()))

def FilterSeq(primer1, primer2,merged_list):
    """
    select sequences between left primer and right primer
    """
    log.write("step2: filter sequence start\t{}\n".format(getTime()))

    # Match primer to sequences
    primer1_r = ''.join([complementary_table[i] for i in primer1[::-1]])
    primer2_r = ''.join([complementary_table[i] for i in primer2[::-1]])
    pattern_1 = re.compile("(\.*{})(\\S+){}\.*".format(primer1,primer2))
    pattern_2 = re.compile("(\.*{})(\\S+){}\.*".format(primer2_r,primer1_r))

    merged_seq_count = Counter(merged_list)
    merged_seq_count = dict(merged_seq_count)
    match_seq_count_dict = dict()
    for seq in set(merged_list):
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
    
    log.write("step2: filter sequence end\t{}\n".format(getTime()))
    return match_seq_count_dict


if __name__ == "__main__":
    group = get_opts()
    opts = group.parse_args()
    proportion = opts.filter
    config = opts.config
    output = opts.output
    
    with open(opts.config, 'r') as f:
        config = json.load(f)
    pe1 = config['pe1']
    pe2 = config['pe2']
    primer1 = config['primer1']
    primer2 = config['primer2']
    p = os.path.basename(pe1).split('.')[0]
    workdir = os.getcwd() + os.sep + output
    if not os.path.exists(workdir):
        print(f'{workdir} already exist')
    else:
        os.system("mkdir {}".format(workdir))
    
    log = open("result.log","w")

    # ******* 1. merge pair-end reads with pear software *******
    result_dir = workdir + os.sep + "1_Assemble_Reads"
    os.system("mkdir -p {}".format(result_dir))

    AssembleReads(pe1, pe2, result_dir)



    # ****** 2. filter sequences between left primer and right primer *******
    
    log.write("step2: filter sequence\t{}\n".format(getTime()))
    # output file
    filter_dir = workdir + os.sep + "2_FilterSeq"
    if os.path.exists(filter_dir):
        print(f'{filter_dir} already exist')
    else:
        os.system("mkdir {}".format(filter_dir))

    filter_low_qual = filter_dir + os.sep + "low_qual_seq.tab"
    filter_low_qual_f = open(filter_low_qual,"w")
    filter_low_qual_f.write("sequence\tcount\tstart_pos\tend_pos\n")

    # mark duplication sequence and count duplication number
    fastq_file = workdir + os.sep + "1_Assemble_Reads" + os.sep + "reads.assembled.fastq"
    command = "awk '{{if(NR%4 == 2) print}}' {} > {}/assemble_seqs.txt".format(
                fastq_file,filter_dir)
    os.system(command)
    
    merged_list = []
    flag = 0
    with open(filter_dir+os.sep+"assemble_seqs.txt") as f:
        for i in f.readlines():
            i = i.strip()
            merged_list.append(i)
            flag += 1
            
    match_seq_count_dict = FilterSeq(primer1, primer2, merged_list)
    seq_count_file = filter_dir + os.sep + "matched_seq_count.txt"
    DictToFile(match_seq_count_dict, seq_count_file)
    log.write("step2: end filter sequence\t{}\n".format(getTime()))
    
    
    # ******* 3. select pass error correct sequences *******
    log.write("step3: check filtered sequence\t{}\n".format(getTime()))
    check_dir = workdir + os.sep + "3_select"
    if os.path.exists(check_dir):
        print(f'{check_dir} already exist')
    else:
        os.system("mkdir {}".format(check_dir))
    match_seq_count_dict = dict()
    with open(seq_count_file,'r') as f:
        data = [i.strip() for i in f.readlines()]
    for i in data:
        k,v = i.split('\t')
        match_seq_count_dict[k] = int(v)
    decoder = Decoder()
    pass_decode_bit, pass_ecc_seq = decoder.check_seq(match_seq_count_dict, check_dir)
    pass_seq_file = check_dir + os.sep + 'pass_ecc.seq'                
    with open(pass_seq_file,'w') as f:
        # table: id \t revert_bit \t count
        for k,v in pass_ecc_seq.items():
            f.write(str(k) + '\t'+'\t'.join(v)+'\n')
    log.write("step3: end\t{}\n".format(getTime()))
    
    
    
    # ******* 4. decode sequences *******
    log.write("step4: decode sequence\t{}\n".format(getTime()))
    decode_dir = workdir + os.sep + "4_decode"
    if os.path.exists(decode_dir):
        print(f'{decode_dir} already exist')
    else:
        os.system("mkdir {}".format(decode_dir))
    
    with open(pass_seq_file,'r') as f:
        data = [i.strip() for i in f.readlines()]
    sum_n = len(data)
    result = decoder.decode(pass_decode_bit)
    file_type, date, file_size, perm, file_addr = decoder.decode_file_attr(
        match_seq_count_dict, sum_n)
    attr_file = decode_dir + os.sep + 'file_attribution.txt'
    
    with open(attr_file,'w') as f:
        f.write("file name:\t" + 'recover.' + file_type + '\n')
        f.write("file creation date:\t"+ date+'\n')
        f.write('file size:\t' + file_size+'\n' )
        f.write('file permission:\t' + perm + '\n')
        f.write('file address:\t'+file_addr)
    
    recover_file = decode_dir + os.sep + 'recover.' + file_type
    with open(recover_file, 'wb') as f:
        f.write(result)
    f.close()
    log.write("step4: end\t{}\n".format(getTime()))
    
    
    