# fastq_subsetter
A lightweight and high-performant C++ tool  for downsampling fastq reads.<br>
Compilation (as done with Ubuntu Linux): <code>g++ -std=c++11 fastq_trimmer.cpp -o fastq_trimmer -lz -pthread -O2</code><br>
Usage: <code>./fastq_trimmer --in INPUT_DIRECTORY --out OUTPUT_DIRECTORY --N TRIM_VALUE</code>
