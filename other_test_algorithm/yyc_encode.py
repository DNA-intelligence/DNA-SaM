'''
FilePath: 11_yy_encode.py
Author: wang yu
Date: 2024-04-25 19:09:38
LastEditTime: 2024-12-11 17:03:51
'''
from yyc import pipeline
from yyc import scheme
import sys
import time

t1 = time.time()
pipeline.encode(method=scheme.YYC(support_bases="A", 
                                  base_reference=[0, 1, 0, 1], 
                                  current_code_matrix=[[1, 1, 0, 0], 
                                                       [1, 0, 0, 1], 
                                                       [1, 1, 0, 0], 
                                                       [1, 1, 0, 0]], 
                                  search_count=100, max_homopolymer=4, 
                                  max_content=0.6, time=t1), 
                input_path=sys.argv[1], 
                output_path=sys.argv[2], 
                model_path="./yyc.pkl", 
                need_index=False, need_log=False, time=t1)
print('[time use]: ', time.time()-t1)