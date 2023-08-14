#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <random>
#include <algorithm>
#include <thread>
#include <future>
#include <regex>
#include <dirent.h>
#include <zlib.h>

std::size_t count_reads(const std::string& file_path) {
    gzFile file = gzopen(file_path.c_str(), "rb");
    if (file == nullptr) {
        throw std::runtime_error("Error opening file: " + file_path);
    }

    std::size_t line_count = 0;
    char buffer[1024];
    while (gzgets(file, buffer, sizeof(buffer)) != nullptr) {
        ++line_count;
    }

    gzclose(file);

    return line_count / 4;
}

void iterate_fastq(const std::string& input_path, const std::set<std::size_t>& indices_to_extract, const std::string& output_path) {
    gzFile file = gzopen(input_path.c_str(), "rb");
    if (file == nullptr) {
        throw std::runtime_error("Error opening file: " + input_path);
    }

    gzFile output_file = gzopen(output_path.c_str(), "wb");  // Open the output file for writing
    if (output_file == nullptr) {
        gzclose(file);
        throw std::runtime_error("Error opening output file: " + output_path);
    }

    char buffer[1024];  // Define the buffer for reading lines
    std::string line;
    std::size_t current_index = 0;
    std::size_t record_length = 0;
    std::vector<std::string> current_record;

    while (gzgets(file, buffer, sizeof(buffer)) != nullptr) {
        line = buffer;
        current_record.push_back(line);
        if (record_length == 3) {
            if (indices_to_extract.find(current_index) != indices_to_extract.end()) {
                for (const auto& record_line : current_record) {
                    gzwrite(output_file, record_line.c_str(), record_line.length());  // Write to the output file
                }
            }
            current_record.clear();
            ++current_index;
            record_length = 0;
        } else {
            ++record_length;
        }
    }

    gzclose(file);
    gzclose(output_file);  // Close the output file when done
}

// Function to subsample a FASTQ file
void subsample_fastq(const std::string& input, const std::string& output, std::size_t num_reads, unsigned int random_seed = 0) {
    std::size_t total_records = count_reads(input);
    
    if (total_records <= num_reads) {
        throw std::runtime_error("Requested number of reads:" + std::to_string(total_records) + " meets or exceeds available number of reads");
    }
    
    std::vector<std::size_t> indices(total_records);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(random_seed));
    indices.resize(num_reads);

    iterate_fastq(input, std::set<std::size_t>(indices.begin(), indices.end()), output);
    //iterate_fastq(input, std::set<std::size_t>(indices.begin(), indices.end()), current_record);
}

std::vector<std::string> getFilesInDirectory(const std::string& directory) {
    std::vector<std::string> files;
    DIR* dir;
    struct dirent* ent;
    if ((dir = opendir(directory.c_str())) != nullptr) {
        while ((ent = readdir(dir)) != nullptr) {
            std::string file_name = ent->d_name;
            if (file_name != "." && file_name != "..") {
                files.push_back(file_name);
            }
        }
        closedir(dir);
    }
    return files;
}
int main() {
    std::string input_dir = "trial2/raw_fastqs/";
    std::string output_dir = "trial2/downsampled_fastqs2/";
    //std::vector<std::size_t> reads = {10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000, 100000};
    std::vector<std::size_t> reads = {500, 1000, 1500, 2000, 3000, 4000, 5000, 110000, 120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000};
    // ... Populate file_names with appropriate file names ...
    // Regular expression pattern for matching file names
    std::regex pattern(".*_R[12]_001\\.fastq(\\.gz)?");

    // Get list of files in the input directory
    std::vector<std::string> file_names = getFilesInDirectory(input_dir);

    // Parallel processing using std::async and std::thread
    std::vector<std::future<void>> futures;

    for (std::size_t num_reads_to_downsample : reads) {
        for (const std::string& file_name : file_names) {
            futures.push_back(std::async(std::launch::async, [input_dir, output_dir, num_reads_to_downsample, file_name]() {
                try {
                    subsample_fastq(input_dir + file_name, output_dir + std::to_string(num_reads_to_downsample) + "_" + file_name, num_reads_to_downsample, 42);
                } catch (const std::exception& e) {
                    std::cerr << "Error processing " << file_name << ": " << e.what() << std::endl;
                }
            }));
        }
    }

    // Wait for all tasks to complete
    for (auto& future : futures) {
        future.wait();
    }

    return 0;
}

