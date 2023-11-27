# fastq_subsetter
The `fastq_subsetter` is a lightweight and high-performant  C++ command-line tool that allows you to efficiently subsample FASTQ files. It can be particularly useful for downsampling large FASTQ files for downstream analysis.

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Usage](#usage)
- [Installation](#installation)
- [Example](#example)
- [License](#license)

## Introduction
The `fastq_subsetter` tool is designed to perform subsampling on FASTQ files, allowing you to extract a specific number of reads or create multiple subsampled versions of a file. The tool supports parallel processing for faster execution and provides options to customize the subsampling process.

## Features
- Efficient subsampling of FASTQ (and gzipped-FASTQ) files.
- Parallel processing using multiple threads for improved performance.
- Customizable subsampling levels.
- Regular expression-based file name filtering.
- Option to skip existing output files to avoid redundant processing. (`force = false`)
- (Future) Ability to set random seeds for reproducible downsampling.

## Usage
To use the `fastq_subsetter` tool, follow these steps:

1. Clone the repository or download the code.
2. Build the executable using a C++ compiler.
3. Run the executable with appropriate command-line arguments.

The available command-line options are as follows:

- `--in/-i <input_dir>`: Specify the input directory containing the FASTQ files.
- `--out/-o <output_dir>`: Specify the output directory for storing the subsampled files.
- `[--regex/r <pattern>]`: Specify the regular expression pattern for matching input file names (optional, defaults to `.*_R[12]_001\.fastq(\.gz)?`).
- `[--start <start>]`: Specify the starting number of reads for subsampling _per fastq file_ (optional, defaults to `0`).
- `[--stop <stop>]`: Specify the stopping number of reads for subsampling _per fastq file_ (optional, defaults to `0`).
- `[--step <step>]`: Specify the step size for subsampling _per fastq file_ (optional, defaults to `0`).
- `[--force]`: Force subsampling even if the output file already exists (optional, defaults to `false`).

Note that incorrect or missing selection of `start`, `stop`, or `step` will cause the program to use default subsampling levels:
```c++
reads = {100, 200, 300, 400, 500, 1000, 1500, 2000, 3000, 4000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000,
                 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000, 10000, 110000, 120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000};
```

## Future updates
In the future, I will incorporate support for the user to specify the random seed (presently, it is static at `42`, so you will get the same results each time), and will permit the specification of non-uniformly distributed subsampling levels from either the command-line or a text file.


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
