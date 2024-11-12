# DNA-SaM
DNA-SaM: a robust system for large-scale data storage

## Background
This is a linear encoding algorithm for DNA information storage, used for DNA encoding of large-scale data.
**----**

## Requirement
* numpy: 1.26.0
* scipy: 1.10.1

## Usage
* **encode**  
// python encode.py -i The_Google_File_System.pdf -o test -c config.json  
> -i: file for DNA encode  
> -o: output file  
> -c: config file which include primer sequence and other information  

* **decode**  
// python decode.py -i .\test\result.txt -o out -c config.json -d /your_path/test  
> -i: decode file include assembled oligo sequences  
> -o: output directory  
> -c: config file which include primer sequence and other information  
> -d: directory of record file and ecc file  

* **decode from fastq file**  
// python decode_sequences.py -c config.json -o out -f 0.2  
> -c: config file which include primer sequence, fastq file and other information  
> -o: output directory  
> -f: value of proportion of low quality bases in the sequence  
