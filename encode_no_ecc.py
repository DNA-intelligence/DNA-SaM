import numpy as np
from DNAException import *
import math
import time
import sys
import random
from datetime import datetime


file_type = {"txt":"ATCG", 'docx':'ATGC', 'log':'TACG', 'jpg':'TCGA', 
             'png':'ACTG', 'gif':'ACGT', 'bmp':'AGCT', 'tif':'AGTC',
             'mp3':'TAGC', 'av':'TCAG', 'aac':'TGAC', 'flac':'TGCA',
             'ogg':'CATG', 'xlsx':'CAGT', 'csv':'CTAG', 'zip':'CTGA',
             'rar':'CGAT', '7z':'CGTA', 'gz':'GATC', 'exe':'GACT',
             'app':'GTCA', 'html':'GTAC', 'css':'GCTA', 'js':'GCAT',
             'db':'GATG', 'sql':'GTAG', 'c':'CATC', 'py':'CTAC',
             'java':'ACGA', 'js':'AGCA', 'ppt':'TCGT', 'pdf':'TGCT',
             'md':'ACAC', 'epub':'CACA'}
units_dict = {
    'B':'TATA', 'KB':'ATAT', 'MB':'AGAG', 'GB':'GAGA', 'TB':'TCTC'
}

# primer_file = ['GCCACAGATTTCGCGATACC','GGAGTCATGTGGGGCCTTTC']  # google file system
# primer_file = ['TACCTCGCACGCACATTAGC','GAACCCAGCCTCAGACCTTG']      # top500
# primer_file = ['CAGAATTCGCGGGTGGTCAC','ATCCTCGCCGAGTGTAAAGG']       # hard disk.png
# primer_file = ['GACCGATACACAGGGCCAAC','TGGGGATCGACTGCACTATG']       # HTML
# primer_file = ['GTTTGCAACGCCGCTCTTTC','GGAGAGGAGCTTGTCGAACC']   # 3D_DNA
primer_file = ['CAGAGTACCCGAGGCGAAAG','GCGCAGGTGTGACAACTTAC']       # css
primer_addr = ['GCAGCGTGTGAGATTCATGG','GCTGTACGTCGCACCTCTAG']
primer_file_info = ['GGCCGCTTTCGTCACATAAC','CGAGACCCGCAACTTGACTG']
primer_file_xor = ['CGATCGTTGCACCTCTTGAC','TCTGTGACTGCGACCCAAGC']
primer_merge = ['AGGATCATCTGAGGCGCAAC','AAGACAAACGGGCAGCGTCC']

class BaseCoder(object):
    def __init__(self, forbidden_list, GAddress, file_name, perm, file_len):
        # self.msg = msg
        try:
            self.GAddress = bytes(GAddress, encoding='utf-8')
            self.slice = self._read_in_slice()
            self.slice_map, self.slice_inv_map = self._read_in_slice_map()
            self.slice_inv_index = self._build_inv_index()
            self.base_xor = [[0, 1, 2, 3], [1, 0, 3, 2],
                             [2, 3, 0, 1], [3, 2, 1, 0]]
            self.max_address_index = len(self.slice_map) - 1
            self.file_id, self.file_suffix = file_name.split(".")
            self.perm = perm
            self.file_len = file_len

        except ConfigFileException as ce:
            print(ce)
            exit(0)

        self.bases = ''
        
    def _read_in_slice(self, slice_file='slice.txt'):
        # 示例：12300131\n
        try:
            sf = open(slice_file, 'r')
            slice = sf.readlines()
            sf.close()
            slice = [s[:-1] for s in slice]
            # slice = [s[10000] for s in slice]
        except:
            raise ConfigFileException('slice file error\n')
        else:
            return slice

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
        base_len = 60
        frag_num = 4
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

    def encode_GAddress(self):
        frag_len = 30
        frag = self.GAddress.ljust(frag_len, b'\x00')
        # addr = self.slice[self.slice_map[self.max_address_index]] + self.slice[self.slice_map[self.max_address_index]]
        new_base = self.encode_frag(frag)
        return self.add_ecc(new_base)
        

    def base_quad(self, base):
        quads = ''
        base_dict = {'A':'0', 'T':'1', 'C':'2', 'G':'3'} 
        for i in range(len(base)):
            quads += base_dict[base[i]]
        return quads
    
    def redundancy(self, bases):
        if len(bases) == 0:
            return None

        ini = bytearray(b''.join([b'\x00'] * len(bases[0])))
        for b in bases:
            ini = bytearray(x ^ y for x, y in zip(ini, bytearray(b)))

        return bytes(ini)
        
    def _read_in_slice_map(self, slice_map_file='slice_map.txt'):
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
    
    def quad_base(self, aa):
        bases = ''
        quad_dict = {'0':'A', '1':'T', '2':'C', '3':'G'}
        for i in range(len(aa)):
            bases += quad_dict[aa[i]]
        return bases
            
    def seq_xor(self,bases):
        result = ''
        for i in range(len(bases[0])):
            x = 0
            for b in bases:
                x = self.base_xor[int(b[i])][x]
            result = result + str(x)
        return result
    
    def encode(self, message=b'', output=''):
        # 切分成每30个字符一份，不足30个的会补齐0，所以最后解码的结果可能会多几个byte
        frag_len = 60
        frags = [message[i:i + frag_len].ljust(frag_len, b'\x00')
                    for i in range(0, len(message), frag_len)]
        result = list()
        data_seq = list()
        addr_idx = 0
        ecc = dict()  # [n, N, l]

        data_block = list()
        tmp_redundancy = list()
        
        # 数据链
        for i, f in enumerate(frags):
            if len(data_block) == 3:       # 每三个frags作为一组进行封装编码
                addr = self.slice[self.slice_map[addr_idx]]
                tmp_redundancy.append(data_block)
                ecc[addr] = []
                data_base = ''
                for data in data_block:
                    new_data = self.encode_frag(data)
                    data_seq.append(self.quad_base(new_data))
                    self.bases = new_data
                    
                data_base = self.quad_base(addr + data_base)
                result.append(data_base)
                # addr_idx += 1
                if len(tmp_redundancy) == 5:
                    addr = self.slice[self.slice_map[addr_idx]]
                    # addr_idx += 1
                    ecc[addr] = []
                    for i in range(len(tmp_redundancy[0])):
                        tmp_list = [lst[i] for lst in tmp_redundancy]
                        xor_base = self.redundancy(tmp_list)
                        # print(len(xor_base))
                        xor_base = self.encode_frag(xor_base)
                        data_seq.append(self.quad_base(xor_base))
                    result.append(xor_base)
                    tmp_redundancy = []
                data_block = []
                
            data_block.append(f)                

        with open(output, 'w') as f:
            for i, seq in enumerate(data_seq):
                f.write(seq+'\n')
                    

def main():
    # print(sys.argv)
    # if len(sys.argv) != 3:
    #     print('args error')
    #     exit(0)

    # input_data = sys.argv[1]
    # input_forbidden_list = sys.argv[2]
    # input_data = '98873709_p0.jpg'
    input_data = sys.argv[1]
    output = sys.argv[2]
    # input_data = 'dna.png'
    input_forbidden_list = 'forbidden_list.txt'

    perm = "user"
    bs = open(input_data, 'rb').read()
    encoder = BaseCoder(input_forbidden_list, '111111111', input_data, perm, len(bs))
    # encoder.encode_frag(b'\x0a\x0a\x0a\x0a\x0a\x0a\x0a\x0a\x0a\x0a\x0a\x0a\x0a\x0a\x0a')
    T1 = time.time()
    encoder.encode(bs, output)
    T2 = time.time()
    print('[time use]: ',T2 - T1)
    # print('程序运行时间:%s毫秒' % ((T2 - T1)*1000))



    # print(coded_msg, record)


if __name__ == '__main__':
    main()
