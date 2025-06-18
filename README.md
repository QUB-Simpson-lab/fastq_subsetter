# 🧬 fastq_subsetter 🧬 
A fast, multithreaded tool for downsampling gzipped (and soon plain) FASTQ files to normalize sequencing depth and streamline comparative analyses.

---

## 🚧 Development Status  
**This tool is under active development.** It currently supports gzipped FASTQ files only. Output concurrency and file-locking are handled via mutexes, but users may observe slight overlaps in logging or edge cases under heavy threading. These will be addressed in future releases.

---

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Usage](#usage)
- [Installation](#installation)
- [Example](#example)
- [License](#license)

## Latest Version
v0.4 (18 June 2025)

## Introduction  
The `fastq_subsetter` tool is designed for efficient subsampling of FASTQ files. It allows users to extract specific numbers of reads or generate multiple subsampled versions of the same file. This is especially useful for **normalizing sequencing depth** across samples, enabling **equitable comparisons** — for example, avoiding bias when analyzing datasets with variable read counts, or evaluating how read depth affects outcomes like variant calling or taxonomic resolution.

Written in C++ for performance, the tool supports multithreaded execution and customizable options for read counts, ratios, and reproducibility.

---

## Features  
- 🧬 Efficient subsampling of gzipped FASTQ files  
- ⚙️ Multithreaded execution using POSIX threads  
- 🔍 Regular expression-based filename filtering  
- 🧪 Paired-end and single-end read support  
- ♻️ Skip or overwrite existing files (`--force`)  
- 🎲 Seeded randomness for reproducible results  
- 🔢 Support for read-count ranges (`--start`, `--stop`, `--step`)  
- 🔻 Ratio-based downsampling (e.g. 1/4 or 25%)  
- 📦 Caching of read counts for performance boost when repeatedly subsampling
- 📄 (Planned) Support for plain (non-gzipped) FASTQ files  
- 📊 (Planned) Custom read level input via text file or CLI  

---


## Usage
To use the `fastq_subsetter` tool, follow these steps:

1. Clone the repository or download the code.
2. Build the executable using a C++ compiler.
3. Run the executable with appropriate command-line arguments.

The available command-line options are as follows:


- `--in`, `-i <input_dir>`: Directory containing input gzipped FASTQ files.
- `--out`, `-o <output_dir>`: Directory to save output files.
- `--regex`, `-p <pattern>`: Regular expression to match input filenames (default: `.*_R[12]_001\.fastq(\.gz)?`).
- `--start`, `-b <start>`: Start value for read count subsampling (default: 0).
- `--stop`, `-e <stop>`: Stop value for read count subsampling (default: 0).
- `--step`, `-c <step>`: Step size for read count subsampling (default: 0).
- `--ratio`, `-r <ratio>`: Ratio of reads to retain (e.g., 0.25 for 25%, or `4` for 1/4).
- `--seed`, `-s <seed>`: Seed for random number generator (default: 0).
- `--force`, `-f`: Overwrite output files if they already exist.
- `--unpaired`, `-u`: Process files individually, not as read pairs.

Note that incorrect or missing selection of `start`, `stop`, or `step` will cause the program to use default subsampling levels:
```c++
reads = {100, 200, 300, 400, 500, 1000, 1500, 2000, 3000, 4000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000,
                 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000, 100000, 110000, 120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000};
```

## Installation
Before compiling, ensure you have the necessary dependencies installed:
- g++ compiler (if not installed: `sudo apt install build-essential`)
- zlib library (`-lz` (if not installed: `sudo apt install zlib1g-dev`))
- POSIX threads (`-pthread` (if not installed: `sudo apt install libpthread-stubs0-dev`))

Compile the tool using the provided command after navigating to a directory containing the source file:
```sh
g++ -std=c++11 fastq_subsetter.cpp -o fastq_subsetter -lz -pthread -O2
```
Alternatively, you may follow this recipe as a guide:
To build the `fastq_subsetter` executable, you need a C++ compiler and the necessary dependencies (e.g., zlib and POSIX threads). Follow this recipe:
```sh
# Clone the repo:
git clone https://github.com/QUB-Simpson-lab/fastq_subsetter.git
# Navigate to the repository directory:
cd fastq_subsetter
# Compile the code:
g++ -o fastq_subsetter fastq_subsetter.cpp -lz -pthread -O2
# The `fastq_subsetter` executable will be generated in the same directory, and can be installed:
sudo cp fastq_subsetter /usr/local/bin/.
```

## Example

Suppose you have a directory named `input_fastq` containing your FASTQ files, and you want to subsample them with varying read numbers. You also want to store the subsampled files in the `output_fastq` directory. Here's how you can use the `fastq_subsetter` tool:

```bash
./fastq_subsetter --in input_fastq --out output_fastq --start 100 --stop 1000 --step 100
```

This command will subsample the input FASTQ files, generating multiple subsampled versions with read counts ranging from `100`_start_ to `1000` _stop_ (inclusive), incrementing by `100` _step_ in between. The resulting files will be stored in the `output_fastq` directory

## Subsampling with Ratios
When using the `--ratio` option, the tool will instead process each file a single time to achieve a fractional value of the original number of reads. Users can specify this as either a decimal fraction or as an integer that will be converted into a fraction. For example:

Inputting 16 will yield 1/16, equivalent to 0.0625.
Inputting 4 will yield 1/4, equivalent to 0.25:

### Downsample each file to 25% of its reads:
```sh
./fastq_subsetter -i input_fastq -o output_fastq -r 0.25
```
### Equivalent alternative:
```sh
./fastq_subsetter -i input_fastq -o output_fastq -r 4  # Interpreted as 1/4
```

## Planned future updates
- Non-Uniform Subsampling Levels: Enable users to specify non-uniform range-distributed subsampling levels through the command line or a text file, offering greater flexibility in subsampling strategies.
- Input Format Support: Allow the use of non-gzipped FASTQ files as input, broadening compatibility with various data sources.

## License
GNU General Public License v3 (GPL-3.0) - Simpson Lab @ Queen’s University Belfast