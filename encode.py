from util import SAM
import numpy as np
from util.DNAException import *
import single_edit_code as sec
import argparse
import math
import time
import sys
import os
import random
from datetime import datetime
import json


file_type = {"txt":"ATCG", 'docx':'ATGC', 'log':'TACG', 'jpg':'TCGA', 
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

primer_addr = ['GCAGCGTGTGAGATTCATGG','GCTGTACGTCGCACCTCTAG']
primer_file_info = ['GGCCGCTTTCGTCACATAAC','CGAGACCCGCAACTTGACTG']
primer_file_xor = ['CGATCGTTGCACCTCTTGAC','TCTGTGACTGCGACCCAAGC']

def get_opts():
    opts = argparse.ArgumentParser()
    opts.add_argument('-i','--input',required=True,help='input file')
    opts.add_argument('-o','--output',required=True,help='output file')
    opts.add_argument('-c','--config',required=True,help='config file of primer')
    
    return opts.parse_args()

class BaseCoder(object):
    def __init__(self, forbidden_list, file_name, config, file_len):
        # self.msg = msg
        try:
            self.slice = self._read_in_slice()
            self.slice_map, self.slice_inv_map = self._read_in_slice_map()
            self.slice_inv_index = self._build_inv_index()
            self.forbidden_list = self._read_in_forbidden_list(forbidden_list)
            self.replace_list = self._read_in_replace_list()
            self.base_xor = [[0, 1, 2, 3], [
                1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]
            self.max_address_index = len(self.slice_map) - 1
            self.file_id, self.file_suffix = file_name.split(".")
            self.config = config
            self.file_len = file_len

        except ConfigFileException as ce:
            print(ce)
            exit(0)

        self.replace_index = 0
        self.bases = ''
        self.replaced_list = []  # replace -> forbidden
        self.sam = SAM.SAM()
        self.code = sec.SingleEditCode(int(48 * math.log(140, 4)))  # 可改

    def _read_in_slice(self, slice_file='slice.txt'):
        # 示例：12300131\n
        try:
            print('util' + os.sep + slice_file)
            sf = open('util' + os.sep + slice_file, 'r')
            slice = sf.readlines()
            sf.close()
            slice = [s.strip() for s in slice]
            # slice = [s[10000] for s in slice]
        except:
            raise ConfigFileException('slice file error\n')
        else:
            return slice

    def _read_in_forbidden_list(self, forbidden_list_file='forbidden_list.txt'):
        # 示例：12300131\n
        try:
            ff = open(forbidden_list_file, 'r')
            forbidden_list = ff.readlines()
            ff.close()
            forbidden_list = [f[:-1] for f in forbidden_list]
        except:
            raise ConfigFileException('forbidden list file error\n')

        return forbidden_list

    def _read_in_replace_list(self, replace_list_file='replace_list.txt'):
        # 示例：31995 37\n
        try:
            rf = open('util' + os.sep + replace_list_file, 'r')
            replace_list = rf.readlines()
            rf.close()
            replace_list = [int(r.split(' ')[0]) for r in replace_list]
        except:
            raise ConfigFileException('replace list file error\n')

        return replace_list

    def _build_inv_index(self):
        try:
            slice_inv_index = dict()
            for id, s in enumerate(self.slice):
                slice_inv_index[s] = id
        except:
            raise ConfigFileException('slice file error\n')
        else:
            return slice_inv_index

    # 加地址时要加到最大的那个
    # 每30字节编码为128个碱基
    def encode_frag(self, bs):  #输入bs为30字节
        bit_len = 8
        base_len = 30
        frag_num = 2
        frag_len = int(base_len / frag_num)
        assert len(bs) == base_len
        msgs = [bs[i:i + frag_len] for i in range(0, len(bs), frag_len)]

        result = ''

        for msg in msgs:    #msgs为包含两个15字节的列表
            bits = list()
            for c in msg:   #将每个字节转换成8位二进制
                bits_per_byte = list()
                for i in range(bit_len):
                    bits_per_byte.append('1' if c & (1 << i) else '0')
                bits_per_byte.reverse()
                bits.extend(bits_per_byte)

            binary_to_base = list()
            for i in range(int(len(bits) / frag_len)):  #将15个字节变成8个15位二进制组成的列表
                binary_to_base.append(bits[i * frag_len: (i + 1) * frag_len])

            binary_to_base = [''.join(c) for c in binary_to_base]   #15位二进制组成字符串
            binary_to_base = [int(c, 2) for c in binary_to_base]    #字符串变成十进制
            binary_to_base = [self.slice[self.slice_map[i]]         #十进制映射成碱基
                              for i in binary_to_base]
            binary_to_base = ''.join(binary_to_base)
            result = result + binary_to_base                        #最终返回16个8 bit序列

        return result

    # 将新增加的碱基序列加入SAM
    def sam_extend(self, str):
        for c in str:
            self.sam.extend(c)

    # 检查当前的碱基序列有无禁忌序列，有则替换
    def check_forbidden(self):
        cur_bases = 'a'
        # forbidden replace index
        replace_record = list()
        while cur_bases != self.bases:
            cur_bases = self.bases
            # forbidden_extra = []
            for i, f in enumerate(self.forbidden_list):     #后缀自动机遍历forbidden
                found, idx = self.sam.traverse(f)
                if found:
                    replace_base = ''
                    if f not in self.replaced_list:           #生成替换列表并保存
                        if self.replace_index == len(self.replace_list):
                            self.replace_index = 0
                        replace_base = self.slice[self.replace_list[self.replace_index]]
                        self.replaced_list.append(f)
                    else:
                        if self.replace_index + 1 == len(self.replace_list)-2:
                            self.replace_index = 0
                        else:
                            self.replace_index += 1
                        replace_base = self.slice[self.replace_list[self.replace_index]]
                        # forbidden_extra.append(replace_base)

                    # TODO 替换的位置有待商榷
                    new_base = self.bases[:idx] + \
                        replace_base + self.bases[idx + len(f):]
                    self.bases = new_base
                    replace_record.append(
                        [str(i), str(len(replace_base)), str(idx)])

                    self.sam = SAM.SAM()
                    self.sam_extend(self.bases)
                    break

            # self.forbidden_list.extend(forbidden_extra)

        return replace_record

    def base_quad(self, base):
        quads = ''
        base_dict = {'A':'0', 'T':'1', 'C':'2', 'G':'3'} 
        for i in range(len(base)):
            quads += base_dict[base[i]]
        return quads
    
    def convert_size(self, size):
        units = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
        index = 0
        while size >= 1024 and index < len(units) -1:
            size /= 1024
            index += 1
        num = f"{size:.4f}"
        num_pre, num_suf = num.split('.')
        pre_nt = self.slice[self.slice_map[int(num_pre)]]
        suf_nt = self.slice[self.slice_map[int(num_suf)]]
        
        return pre_nt + suf_nt + self.base_quad(units_dict[units[index]])

    def add_ecc(self, raw_data=''):
        try:
            raw_data_string = sec.QaryString(4, list(raw_data))
            # code = sec.SingleEditCode(int(36.0 * math.log(len(raw_data_string), 4)))
            code_data, n, N, l = self.code.encode(raw_data_string)
            # r_data = self.code.decode(code_data, n, N, l, verbose=False)
            # print(r_data.__eq__(raw_data_string))
            code_data = ''.join([str(c.val[0]) for c in code_data])
        except Exception:
            raise EncodeException('encode ecc error\n')
        else:
            return code_data, [n, N, l]

    def redundancy(self, bases):
        if len(bases) == 0:
            return None

        ini = bytearray(b''.join([b'\x00'] * len(bases[0])))
        for b in bases:
            ini = bytearray(x ^ y for x, y in zip(ini, bytearray(b)))

        return bytes(ini)
        
    def replace_forb(self, base, addr, records):
        self.bases = base
        replace_record = self.check_forbidden()
        self.sam = SAM.SAM()
        self.sam_extend(base)
        replace_record = self.check_forbidden()
        if replace_record:
            replace_record.append(addr)
            replace_record.reverse()
            records.append(replace_record)
        return records

    def _read_in_slice_map(self, slice_map_file='slice_map.txt'):
        # 示例：31995 37\n
        try:
            mf = open('util' + os.sep + slice_map_file, 'r')
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
    
    def quad_base(self, aa):
        bases = ''
        quad_dict = {'0':'A', '1':'T', '2':'C', '3':'G'}
        for i in range(len(aa)):
            bases += quad_dict[aa[i]]
        return bases
            
    def seq_xor(self,bases):
        result = ''
        if len(bases[0]) > len(bases[1]):
            bases[1] = bases[1].ljust(len(bases[0]), '0')
        elif len(bases[0]) < len(bases[1]):
            bases[0] = bases[0].ljust(len(bases[1]), '0')
        for i in range(len(bases[0])):
            x = 0
            for b in bases:
                x = self.base_xor[int(b[i])][x]
            result = result + str(x)
        return result
    
    def encode(self, message=b''):
        # 切分成每30个字符一份，不足30个的会补齐0，所以最后解码的结果可能会多几个byte
        try:
            frag_len = 30
            frags = [message[i:i + frag_len].ljust(frag_len, b'\x00')
                     for i in range(0, len(message), frag_len)]
            result = list()
            records = list()
            addr_idx = 0
            ecc = dict()  # [n, N, l]

            code_frags_redundancy = list()
            for i, f in enumerate(frags):
                # 通过校验码设置冗余
                if i != 0 and i % 5 == 0:
                    xor_base = self.redundancy(code_frags_redundancy)
                    # addr_high = math.floor(addr_idx / self.max_address_index)
                    # addr_low = addr_idx % self.max_address_index
                    addr = self.slice[self.slice_map[addr_idx]]
                    addr_idx += 1
                    new_base = self.encode_frag(xor_base)
                    # records = self.replace_forb(new_base, addr, records)        #替换禁忌序列
                    self.bases, ecc[addr] = self.add_ecc(new_base)
                    records = self.replace_forb(self.bases, addr, records)
                    xor_base = addr + self.bases
                    xor_base = self.quad_base(xor_base)
                    xor_base = self.config['primer1'] + xor_base + self.config['primer2']
                    result.append(xor_base)
                    code_frags_redundancy.clear()

                # 替换禁忌序列
                addr = self.slice[self.slice_map[addr_idx]]
                new_base = self.encode_frag(f)      #字节转换成碱基
                self.bases = new_base
                # records = self.replace_forb(new_base, addr, records)
                # self.sam_extend(new_base)
                # replace_record = self.check_forbidden()

                # if replace_record:
                #     replace_record.append(addr)
                #     replace_record.reverse()
                #     records.append(replace_record)

                code_frags_redundancy.append(f)

                # 添加纠错码
                self.bases, ecc[addr] = self.add_ecc(new_base)
                records = self.replace_forb(self.bases, addr, records)
                data_base = addr + self.bases
                data_base = self.quad_base(data_base)
                data_base = self.config['primer1'] + data_base + self.config['primer2']
                result.append(data_base)

                addr_idx += 1

            # merge_seq = []
            
            # 地址链
            rand_list = [random.randint(1,1024) for i in range(11)]
            rand_list = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            print(rand_list)
            pre_addr = ''
            addr = self.slice[self.slice_map[self.max_address_index]]
            addr_idx += 1
            for _addr in rand_list:
                pre_addr += self.slice[self.slice_map[_addr]]
            # pre_addr += '01230'
            g_addr = pre_addr
            g_addr_base, ecc[addr] = self.add_ecc(g_addr)
            records = self.replace_forb(g_addr_base, addr, records)
            seqs = [self.bases[:-17]]
            g_addr = self.slice[self.slice_map[self.max_address_index]]+self.bases
            
            
            g_addr_base = self.quad_base(g_addr)
            g_addr_base = self.config['primer1'] + primer_addr[0] + g_addr_base + primer_addr[-1] + self.config['primer2']
            # merge_seq.append(g_addr_base)
            result.append(g_addr_base)

            # 文件信息链
            fileType_base= file_type[self.file_suffix]
            fileType_base = self.base_quad(fileType_base)       #文件类型

            size_base = self.convert_size(self.file_len)        #文件大小
            
            ms_a = 16                                       #文件ID 1
            ms_b = 6179                                        #文件ID 2   ID=1+2随机数组成
            ms_a_base = self.slice[self.slice_map[int(ms_a)]]
            ms_b_base = self.slice[self.slice_map[int(ms_b)]]
            ms_base = ms_a_base + ms_b_base
                
            date = datetime.now()
            date = str(date)
            date, tim = date.split(' ')
            tim = tim.split('.')[0]
            y, m, d = date.split('-')
            h, mi = tim.split(':')[:2]
            year = int(y)
            md = int(m+d)
            hm = int(h+mi)
            y_base = self.slice[self.slice_map[int(year)]]
            md_base = self.slice[self.slice_map[int(md)]]
            hm_base = self.slice[self.slice_map[int(hm)]]
            date_base = y_base + md_base + hm_base          #文件日期
            
            if self.config['permission'] == 'user':                         #文件权限
                # perm_base = "GCTCACGTTCAATTCGGCGC"
                perm_base = '32120231120011233232'
            elif self.config['permission'] == 'admin':
                # perm_base = "CAAGACGCGGAGACAGGTCC"
                perm_base = '20030232330302033122'
            else:
                # perm_base = "ATTCCGGCTGTGAGGTGTGG"
                perm_base = '01122332131303313133'
                
            # rand_base = "TACG"
            rand_base = '1023'
            file_base = ms_base + fileType_base + date_base + size_base + rand_base + perm_base
            
            addr = self.slice[self.slice_map[addr_idx]]
            addr_idx += 1
            file_base, ecc[addr] = self.add_ecc(file_base)
            records = self.replace_forb(file_base, addr, records)
            file_base = self.bases
            seqs.append(file_base[:-17])
            file_base = addr + file_base
            
            file_base = self.quad_base(file_base)
            file_base = self.config['primer1'] + primer_file_info[0] + file_base + primer_file_info[-1] + self.config['primer2']
            result.append(file_base)
            
            # 文件信息异或
            addr = self.slice[self.slice_map[addr_idx]]
            file_xor_base = self.seq_xor(seqs)
            file_xor_base, ecc[addr] = self.add_ecc(file_xor_base)
            records = self.replace_forb(file_xor_base, addr, records)
            
            file_xor_base = addr + self.bases
            
            file_xor_base = self.quad_base(file_xor_base)   
            file_xor_base = self.config['primer1'] + primer_file_xor[0] + file_xor_base + primer_file_xor[-1] +self.config['primer2']   #反向引物错了与正向相同
            result.append(file_xor_base)
            
        except Exception:
            raise EncodeException('encode error\n')
        else:
            return result, records, ecc


def main():
    opts = get_opts()
    input = opts.input
    output = opts.output
    with open(opts.config, 'r') as f:
        config = json.load(f)
    try:
        os.mkdir(output)
    except:
        pass
    input_forbidden_list = 'util' + os.sep + 'forbidden_list.txt'

    bs = open(input, 'rb').read()
    encoder = BaseCoder(input_forbidden_list, input, config, len(bs))
    T1 = time.time()
    coded_msg, record, ecc = encoder.encode(bs)
    T2 = time.time()
    print(T2 - T1)

    record_save = list()
    for r in record:
        s = ''
        for id, item in enumerate(r):
            if id == 0:
                s = item
            else:
                s = s + ','.join(item)

            s = s + '-'

        record_save.append(s)

    ecc_save = list()
    for e in ecc:
        s = str(e)
        for l in ecc[e]:
            s = s + '-' + str(l)

        ecc_save.append(s)

    a = '\n'.join(ecc_save)

    wf = open(output + os.sep + 'result.txt', 'w')
    rf = open(output + os.sep + 'record.txt', 'w')
    ef = open(output + os.sep + 'ecc.txt', 'w')

    wf.write('\n'.join(coded_msg) + '\n')
    rf.write('\n'.join(record_save) + '\n')
    ef.write('\n'.join(ecc_save) + '\n')

    wf.close()
    rf.close()
    ef.close()

    # print(coded_msg, record)


if __name__ == '__main__':
    main()
