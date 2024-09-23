#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <random>
#include <algorithm>
#include <thread>
#include <future>
#include <regex>
#include <unordered_map>
#include <dirent.h>
#include <zlib.h>
#include <unistd.h>
#include <getopt.h>

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

void subsample_fastq(const std::string& input, const std::string& output, const std::size_t& num_reads, std::unordered_map<std::string, std::size_t>& cache, const std::size_t& random_seed = 0) {
    std::size_t total_records = count_reads(input);

    if (total_records < num_reads) {
        throw std::runtime_error("Requested number of reads:" + std::to_string(num_reads) + " exceeds available number of reads: " + std::to_string(total_records));
    }

    std::vector<std::size_t> indices(total_records);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(random_seed));
    indices.resize(num_reads);

    iterate_fastq(input, std::set<std::size_t>(indices.begin(), indices.end()), output);
}

void process_file_range(const std::string& input_path, const std::string& output_dir, const std::size_t& num_reads_to_downsample, std::unordered_map<std::string, std::size_t>& cache, const std::size_t& seed) {
    std::string file_name = input_path.substr(input_path.find_last_of('/') + 1);
    std::string output_path = output_dir + std::to_string(num_reads_to_downsample) + "_" + file_name;

    try {
        subsample_fastq(input_path, output_path, num_reads_to_downsample, cache, seed);
    } catch (const std::exception& e) {
        std::cerr << "Error processing " << file_name << ": " << e.what() << std::endl;
    }
}

void subsample_fastq_ratio(const std::string& input, const std::string& output, const double& ratio, const std::size_t& random_seed = 0) {
    std::size_t total_records = count_reads(input);

    std::vector<std::size_t> indices(total_records);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(random_seed));
    indices.resize(static_cast<std::size_t>(total_records * ratio));

    iterate_fastq(input, std::set<std::size_t>(indices.begin(), indices.end()), output);
}

void process_file_ratio(const std::string& input_path, const std::string& output_dir, const double& ratio, const std::size_t& seed ) {
    std::string file_name = input_path.substr(input_path.find_last_of('/') + 1);
    std::string output_path = output_dir + "ratio=" + std::to_string(ratio) + "_" + file_name;

    try {
        subsample_fastq_ratio(input_path, output_path, ratio, seed);
    } catch (const std::exception& e) {
        std::cerr << "Error processing " << file_name << ": " << e.what() << std::endl;
    }
}


std::vector<std::string> getFilesInDirectory(const std::string& directory, const std::regex& pattern) {
    std::vector<std::string> files;
    DIR* dir;
    struct dirent* ent;
    if ((dir = opendir(directory.c_str())) != nullptr) {
        while ((ent = readdir(dir)) != nullptr) {
            std::string file_name = ent->d_name;
            if (file_name != "." && file_name != ".." && std::regex_match(file_name, pattern)) {
                files.push_back(directory + file_name);
            }
        }
        closedir(dir);
    }
    return files;
}


int main(int argc, char* argv[]) {
    std::string input_dir, output_dir;
    std::vector<std::size_t> reads;
    std::string pattern_str = ".*_R[12]_001\\.fastq(\\.gz)?";
    std::size_t start = 0, stop = 0, step = 0;
    std::size_t seed = 0;
    double ratio = 0.0;
    bool force = false;

    // Define long-named options
    static struct option long_options[] = {
        {"in", required_argument, 0, 'i'},
        {"out", required_argument, 0, 'o'},
        {"regex", required_argument, 0, 'p'},
        {"start", required_argument, 0, 'b'},
        {"stop", required_argument, 0, 'e'},
        {"step", required_argument, 0, 'g'},
        {"seed", required_argument, 0,'s'},
        {"ratio", required_argument, 0, 'r'},
        {"force", no_argument, 0, 'f'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "i:o:p::b::e::g::s::r::f", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'i':
                input_dir = optarg;
                break;
            case 'o':
                output_dir = optarg;
                break;
            case 'p':
                pattern_str = optarg;
                break;
            case 'b':
                start = std::stoul(optarg);
                break;
            case 'e':
                stop = std::stoul(optarg);
                break;
            case 'g':
                step = std::stoul(optarg);
                break;
            case 's':
                seed = std::stoul(optarg);
                break;
            case 'r':
                ratio = std::stod(optarg);
                break;
            case 'f':
                force = true;
                break;
            default:
                std::cerr << "Usage: " << argv[0] << " --in input_dir --out output_dir [--regex pattern] [--start start] [--stop stop] [--step step]\n";
                return 1;
        }
    }

    std::cout << "Input directory: " << input_dir << std::endl;
    std::cout << "Output directory: " << output_dir << std::endl;
    std::cout << "Regex pattern: " << pattern_str << std::endl;


    if(ratio == 0)  {
        // Generate downsampling levels
        if (start > 0 && stop > 0 && step > 0) {
            for (std::size_t i = start; i <= stop; i += step) {
                reads.push_back(i);
            }
        } else {
            std::cerr << "Warning: using default reads." << std::endl;
            reads = {100, 200, 300, 400, 500, 1000, 1500, 2000, 3000, 4000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000,
                     55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000, 100000, 110000, 120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000};
        }
    } else {
        if (ratio > 1) {
            std::cout << "--ratio/-r passed:" << ratio << " - setting to fractional form: " << 1/ratio << std::endl;
            ratio = 1/ratio;
        }
    }

    // Regular expression for file matching
    std::regex pattern(pattern_str);
    auto file_names = getFilesInDirectory(input_dir, pattern);

    // Parallel processing
    std::vector<std::future<void>> futures;
    std::size_t num_cores = std::thread::hardware_concurrency();
    std::size_t completed_iterations = 0;

    if (ratio == 0) {
        std::size_t total_iterations = reads.size() * file_names.size();
        // Cache for original read counts for 'range' mode.
        std::unordered_map<std::string, std::size_t> read_counts;

        for (std::size_t num_reads_to_downsample : reads) {
            for (const std::string& file_name : file_names) {
                std::string output_file_name = output_dir + std::to_string(num_reads_to_downsample) + "_" + file_name.substr(file_name.find_last_of('/') + 1);
                bool file_exists = (access(output_file_name.c_str(), F_OK) == 0);

                if (!force && file_exists) {
                    std::cout << "Output file already exists for " << file_name << ", skipping..." << std::endl;
                    ++completed_iterations;
                    continue; // Skip processing this iteration
                }

                if (futures.size() >= num_cores) {
                    futures.front().wait();
                    futures.erase(futures.begin());
                }

                std::cout << "Processing: " << file_name << " (" << num_reads_to_downsample << " reads)" << std::endl;
                futures.push_back(std::async(std::launch::async, [&]() {
                    process_file_range(file_name, output_dir, num_reads_to_downsample, read_counts, seed);
                    ++completed_iterations;
                    std::cout << "Progress: " << completed_iterations << "/" << total_iterations << " iterations" << std::endl;
                }));
            }
        }

        // Wait for all remaining tasks to complete
        for (auto& future : futures) {
            future.wait();
        }
    } else {
        std::size_t total_iterations = file_names.size();

        for (const std::string& file_name : file_names) {
            std::string output_file_name = output_dir + "ratio=" + std::to_string(ratio) + "_" + file_name.substr(file_name.find_last_of('/') + 1);
            bool file_exists = (access(output_file_name.c_str(), F_OK) == 0);

            if (!force && file_exists) {
                std::cout << "Output file already exists for " << file_name << ", skipping..." << std::endl;
                ++completed_iterations;
                continue; // Skip processing this iteration
            }
            if (futures.size() >= num_cores) {
                futures.front().wait();
                futures.erase(futures.begin());
            }

            std::cout << "Processing: " << file_name << " (" << ratio << " ratio)" << std::endl;
            futures.push_back(std::async(std::launch::async, [&]() {
                process_file_ratio(file_name, output_dir, ratio, seed);
                ++completed_iterations;
                std::cout << "Progress: " << completed_iterations << "/" << total_iterations << " iterations" << std::endl;
            }));
        }

        // Wait for all remaining tasks to complete
        for (auto& future : futures) {
            future.wait();
        }
    }
    return 0;
}