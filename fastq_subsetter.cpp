#include <algorithm>
#include <atomic>
#include <future>
#include <iostream>
#include <mutex>
#include <random>
#include <regex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
#include <dirent.h>
#include <zlib.h>
#include <unistd.h>
#include <getopt.h>

// Coders: Dr Evan Troendle & Dr Stephen J Bridgett
// Software version: v0.4 (18 June 2025)

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

std::size_t get_or_cache_read_count(const std::string& file,
                                    std::unordered_map<std::string, std::size_t>& cache,
                                    std::mutex& cache_mutex) {
    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        auto it = cache.find(file);
        if (it != cache.end()) {
            return it->second;
        }
    }

    // unlock before expensive operation
    std::size_t count = count_reads(file);

    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        cache[file] = count;
    }

    return count;
}


bool should_extract(const std::vector<std::size_t>& sorted_indices, std::size_t index) {
    return std::binary_search(sorted_indices.begin(), sorted_indices.end(), index);
}

void iterate_fastq(const std::string& input_path, const std::vector<std::size_t>& indices_to_extract, const std::string& output_path) {
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
            if (should_extract(indices_to_extract, current_index)) {
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

void subsample_fastq(const std::string& input, const std::string& output,
                     const std::size_t& num_to_sample, 
                     std::unordered_map<std::string, std::size_t>& cache,
                     std::mutex& cache_mutex,
                     const std::size_t& random_seed = 0) {
    std::size_t total_records = get_or_cache_read_count(input, cache, cache_mutex);
    
    if (total_records < num_to_sample) {
        throw std::runtime_error("Requested number of reads:" + std::to_string(num_to_sample) + " exceeds available number of reads: " + std::to_string(total_records));
        return;
    }
    
    std::vector<std::size_t> indices(total_records);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(random_seed));
    indices.resize(num_to_sample);
    std::sort(indices.begin(), indices.end());

    iterate_fastq(input, indices, output);
}

void subsample_paired_fastq(const std::string& input_R1, const std::string& input_R2,
                            const std::string& output_R1, const std::string& output_R2,
                             std::size_t num_to_sample, std::unordered_map<std::string, std::size_t>& cache,
                            std::mutex& cache_mutex, std::size_t random_seed = 0) {
    std::size_t total_records = get_or_cache_read_count(input_R1, cache, cache_mutex);

    if (total_records < num_to_sample) {
        throw std::runtime_error("Requested number of reads:" + std::to_string(num_to_sample) + " exceeds available number of reads: " + std::to_string(total_records));
        return;
    }
    std::vector<std::size_t> indices(total_records);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(random_seed));
    indices.resize(num_to_sample);
    std::sort(indices.begin(), indices.end());

    iterate_fastq(input_R1, indices, output_R1);
    iterate_fastq(input_R2, indices, output_R2);
}


void process_file_range(const std::string& input_path, const std::string& output_dir, const std::size_t& num_reads_to_downsample, std::unordered_map<std::string, std::size_t>& cache, std::mutex& cache_mutex, const std::size_t& seed) {
    std::string file_name = input_path.substr(input_path.find_last_of('/') + 1);
    std::string output_path = output_dir + std::to_string(num_reads_to_downsample) + "_" + file_name;

    try {
        subsample_fastq(input_path, output_path, num_reads_to_downsample, cache, cache_mutex, seed);
    } catch (const std::exception& e) {
        std::cerr << "Error processing " << file_name << ": " << e.what() << std::endl;
    }
}

void process_paired_file_range(const std::string& input_R1, const std::string& input_R2, const std::string& output_dir, const std::size_t& num_reads_to_downsample, std::unordered_map<std::string, std::size_t>& cache, std::mutex& cache_mutex, const std::size_t& seed) {
    std::string file_name_R1 = input_R1.substr(input_R1.find_last_of('/') + 1);
    std::string file_name_R2 = input_R2.substr(input_R2.find_last_of('/') + 1);

    std::string output_R1 = output_dir + std::to_string(num_reads_to_downsample) + "_" + file_name_R1;
    std::string output_R2 = output_dir + std::to_string(num_reads_to_downsample) + "_" + file_name_R2;

    try {
        subsample_paired_fastq(input_R1, input_R2, output_R1, output_R2, num_reads_to_downsample, cache, cache_mutex, seed);
    } catch (const std::exception& e) {
        std::cerr << "Error processing pair " << input_R1 << ' ' << input_R2 << ": " << e.what() << std::endl;
    }
	return;
}

void subsample_fastq_ratio(const std::string& input, const std::string& output,
                     const double& ratio,
                     std::unordered_map<std::string, std::size_t>& cache,
                     std::mutex& cache_mutex,
                     const std::size_t& random_seed = 0) {
    std::size_t total_records = get_or_cache_read_count(input, cache, cache_mutex);
    std::size_t num_to_sample = total_records * ratio;
    if (total_records < num_to_sample) {
        throw std::runtime_error("Requested number of reads:" + std::to_string(num_to_sample) + " exceeds available number of reads: " + std::to_string(total_records));
        return;
    }
    std::vector<std::size_t> indices(total_records);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(random_seed));
    indices.resize(static_cast<std::size_t>(num_to_sample));
    std::sort(indices.begin(), indices.end());

    iterate_fastq(input, indices, output);
    return;
}

void subsample_paired_fastq_ratio (const std::string& input_R1, const std::string& input_R2,
                            const std::string& output_R1, const std::string& output_R2,
                            double ratio, std::unordered_map<std::string, std::size_t>& cache,
                            std::mutex& cache_mutex, std::size_t random_seed = 0) {

    std::size_t total_records = get_or_cache_read_count(input_R1, cache, cache_mutex);
    std::size_t num_to_sample = total_records * ratio;
    if (total_records < num_to_sample) {
        throw std::runtime_error("Requested number of reads:" + std::to_string(num_to_sample) + " exceeds available number of reads: " + std::to_string(total_records));
        return;
    }
    std::vector<std::size_t> indices(total_records);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(random_seed));
    indices.resize(static_cast<std::size_t>(num_to_sample));
    std::sort(indices.begin(), indices.end());

    iterate_fastq(input_R1, indices, output_R1);
    iterate_fastq(input_R2, indices, output_R2);
    return;
}


void process_file_ratio(const std::string& input_path, const std::string& output_dir, const double& ratio, std::unordered_map<std::string, std::size_t>& cache, std::mutex& cache_mutex, const std::size_t& seed ) {
    std::string file_name = input_path.substr(input_path.find_last_of('/') + 1);
    std::string output_path = output_dir + "ratio_" + std::to_string(ratio) + "_" + file_name;

    try {
        subsample_fastq_ratio(input_path, output_path, ratio, cache, cache_mutex, seed);
    } catch (const std::exception& e) {
        std::cerr << "Error processing " << file_name << ": " << e.what() << std::endl;
    }
}

void process_paired_file_ratio(const std::string& input_R1, const std::string& input_R2, const std::string& output_dir, const double& ratio, std::unordered_map<std::string, std::size_t>& cache, std::mutex& cache_mutex, const std::size_t& seed) {
    std::string file_name_R1 = input_R1.substr(input_R1.find_last_of('/') + 1);
    std::string file_name_R2 = input_R2.substr(input_R2.find_last_of('/') + 1);

    std::string output_R1 = output_dir + "ratio_" + std::to_string(ratio) + "_" + file_name_R1;
    std::string output_R2 = output_dir + "ratio_" + std::to_string(ratio) + "_" + file_name_R2;

    try {
        subsample_paired_fastq_ratio(input_R1, input_R2, output_R1, output_R2, ratio, cache, cache_mutex, seed);
    } catch (const std::exception& e) {
        std::cerr << "Error processing pair " << input_R1 << ' ' << input_R2 << ": " << e.what() << std::endl;
    }
	return;
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

std::vector<std::pair<std::string, std::string>> getPairedFiles(const std::string& input_dir, const std::regex& pattern) {
    std::map<std::string, std::pair<std::string, std::string>> paired_map;

    auto files = getFilesInDirectory(input_dir, pattern);

    for (const auto& file : files) {
        std::string basename = file.substr(file.find_last_of('/') + 1);
        std::smatch match;
        std::regex r("(.*)(_R[12]_001\\.fastq(\\.gz)?)");
        if (std::regex_match(basename, match, r)) {
            std::string sample_prefix = match[1];
            std::string read_part = match[2];
            if (read_part.find("_R1_") != std::string::npos) {
                paired_map[sample_prefix].first = file;
            } else if (read_part.find("_R2_") != std::string::npos) {
                paired_map[sample_prefix].second = file;
            }
        }
    }

    std::vector<std::pair<std::string, std::string>> paired_files;
    for (const auto& kv : paired_map) {
        if (!kv.second.first.empty() && !kv.second.second.empty()) {
            paired_files.push_back(kv.second);
        }
    }

    return paired_files;
}

std::vector<std::string> getUnpairedFiles(const std::string& input_dir, const std::regex& pattern) {
    return getFilesInDirectory(input_dir, pattern);
}

void process_unpaired_files(const std::vector<std::string>& file_names,
                            const std::string& output_dir,
                            const std::vector<std::size_t>& reads,
                            double ratio,
                            bool force,
                            std::size_t seed) {
    std::unordered_map<std::string, std::size_t> cache;
    std::mutex cache_mutex;
    std::mutex progress_mutex;
    std::atomic<std::size_t> completed_iterations{0};

    std::vector<std::future<void>> futures;
    std::size_t num_cores = std::thread::hardware_concurrency();

    if (ratio == 0) {
        std::size_t total_iterations = reads.size() * file_names.size();

        for (std::size_t num_reads_to_downsample : reads) {
            for (const std::string& file_name : file_names) {
                std::string base_name = file_name.substr(file_name.find_last_of('/') + 1);
                std::string output_file_name = output_dir + std::to_string(num_reads_to_downsample) + "_" + base_name;
                bool file_exists = (access(output_file_name.c_str(), F_OK) == 0);

                if (!force && file_exists) {
                    std::lock_guard<std::mutex> lock(progress_mutex);
                    std::cout << "Output file already exists for " << base_name << ", skipping..." << std::endl;
                    completed_iterations++;
                    continue;
                }

                if (futures.size() >= num_cores) {
                    futures.front().wait();
                    futures.erase(futures.begin());
                }

                futures.push_back(std::async(std::launch::async,
                    [file_name, output_dir, num_reads_to_downsample, &cache, &cache_mutex, seed, &progress_mutex,
                     &completed_iterations, total_iterations]() {
                        try {
                            process_file_range(file_name, output_dir, num_reads_to_downsample, cache, cache_mutex, seed);
                        } catch (const std::exception& e) {
                            std::lock_guard<std::mutex> lock(progress_mutex);
                            std::cerr << "Error: " << e.what() << std::endl;
                        }

                        completed_iterations++;
                        {
                            std::lock_guard<std::mutex> lock(progress_mutex);
                            std::cout << "Progress: " << completed_iterations << "/" << total_iterations << " iterations" << std::endl;
                        }
                    }
                ));
            }
        }
    } else {
        std::size_t total_iterations = file_names.size();

        for (const std::string& file_name : file_names) {
            std::string base_name = file_name.substr(file_name.find_last_of('/') + 1);
            std::string output_file_name = output_dir + "ratio=" + std::to_string(ratio) + "_" + base_name;
            bool file_exists = (access(output_file_name.c_str(), F_OK) == 0);

            if (!force && file_exists) {
                std::lock_guard<std::mutex> lock(progress_mutex);
                std::cout << "Output file already exists for " << base_name << ", skipping..." << std::endl;
                completed_iterations++;
                continue;
            }

            if (futures.size() >= num_cores) {
                futures.front().wait();
                futures.erase(futures.begin());
            }

            futures.push_back(std::async(std::launch::async,
                [file_name, output_dir, ratio, &cache, &cache_mutex,
                 seed, &progress_mutex, &completed_iterations, total_iterations]() {
                    try {
                        process_file_ratio(file_name, output_dir, ratio, cache, cache_mutex, seed);
                    } catch (const std::exception& e) {
                        std::lock_guard<std::mutex> lock(progress_mutex);
                        std::cerr << "Error: " << e.what() << std::endl;
                    }

                    completed_iterations++;
                    {
                        std::lock_guard<std::mutex> lock(progress_mutex);
                        std::cout << "Progress: " << completed_iterations << "/" << total_iterations << " iterations" << std::endl;
                    }
                }
            ));
        }
    }

    for (auto& fut : futures) {
        fut.wait();
    }
}

void process_paired_files(const std::vector<std::pair<std::string, std::string>>& paired_files,
                           const std::string& output_dir,
                           const std::vector<std::size_t>& reads,
                           double ratio,
                           bool force,
                           std::size_t seed) {
    std::unordered_map<std::string, std::size_t> cache;
    std::mutex cache_mutex;
    std::mutex progress_mutex;
    std::atomic<std::size_t> completed_iterations{0};  // thread-safe

    std::vector<std::future<void>> futures;
    std::size_t num_cores = std::thread::hardware_concurrency();

    if (ratio == 0) {
        std::size_t total_iterations = reads.size() * paired_files.size();
        std::unordered_map<std::string, std::size_t> read_counts;

        for (std::size_t num_reads_to_downsample : reads) {
            for (const auto& pair : paired_files) {
                const auto& file_R1 = pair.first;
                const auto& file_R2 = pair.second;

                std::string base_name_R1 = file_R1.substr(file_R1.find_last_of('/') + 1);
                std::string output_file_name_R1 = output_dir + std::to_string(num_reads_to_downsample) + "_" + base_name_R1;
                std::string base_name_R2 = file_R2.substr(file_R2.find_last_of('/') + 1);
                std::string output_file_name_R2 = output_dir + std::to_string(num_reads_to_downsample) + "_" + base_name_R2;

                bool file_exists = (access(output_file_name_R1.c_str(), F_OK) == 0) &&
                                   (access(output_file_name_R2.c_str(), F_OK) == 0);

                if (!force && file_exists) {
                    std::lock_guard<std::mutex> lock(progress_mutex);
                    std::cout << "Output files already exist for pair: " << base_name_R1 << " and " << base_name_R2 << ", skipping..." << std::endl;
                    completed_iterations++;
                    continue;
                }

                if (futures.size() >= num_cores) {
                    futures.front().wait();
                    futures.erase(futures.begin());
                }

                futures.push_back(std::async(std::launch::async,
                    [file_R1, file_R2, output_dir, num_reads_to_downsample, seed, &cache, &cache_mutex, &progress_mutex, &completed_iterations, total_iterations]() {
                        try {
                            process_paired_file_range(file_R1, file_R2, output_dir, num_reads_to_downsample, cache, cache_mutex, seed);
                        } catch (const std::exception& e) {
                            std::lock_guard<std::mutex> lock(progress_mutex);
                            std::cerr << "Error: " << e.what() << std::endl;
                        }

                        completed_iterations++;
                        {
                            std::lock_guard<std::mutex> lock(progress_mutex);
                            std::cout << "Progress: " << completed_iterations << "/" << total_iterations << " iterations" << std::endl;
                        }
                    }
                ));
            }
        }
    } else {
        std::size_t total_iterations = paired_files.size();

        for (const auto& pair : paired_files) {
            const auto& file_R1 = pair.first;
            const auto& file_R2 = pair.second;

            std::string base_name_R1 = file_R1.substr(file_R1.find_last_of('/') + 1);
            std::string output_file_name_R1 = output_dir + "ratio=" + std::to_string(ratio) + "_" + base_name_R1;
            std::string base_name_R2 = file_R2.substr(file_R2.find_last_of('/') + 1);
            std::string output_file_name_R2 = output_dir + "ratio=" + std::to_string(ratio) + "_" + base_name_R2;

            bool file_exists = (access(output_file_name_R1.c_str(), F_OK) == 0) &&
                               (access(output_file_name_R2.c_str(), F_OK) == 0);

            if (!force && file_exists) {
                std::lock_guard<std::mutex> lock(progress_mutex);
                std::cout << "Output files already exist for pair: " << base_name_R1 << " and " << base_name_R2 << ", skipping..." << std::endl;
                completed_iterations++;
                continue;
            }

            if (futures.size() >= num_cores) {
                futures.front().wait();
                futures.erase(futures.begin());
            }

            futures.push_back(std::async(std::launch::async,
                [file_R1, file_R2, output_dir, ratio, seed, &cache, &cache_mutex, &progress_mutex, &completed_iterations, total_iterations]() {
                    try {
                         process_paired_file_ratio(file_R1, file_R2, output_dir, ratio, cache, cache_mutex, seed);
                    } catch (const std::exception& e) {
                        std::lock_guard<std::mutex> lock(progress_mutex);
                        std::cerr << "Error: " << e.what() << std::endl;
                    }

                    completed_iterations++;
                    {
                        std::lock_guard<std::mutex> lock(progress_mutex);
                        std::cout << "Progress: " << completed_iterations << "/" << total_iterations << " iterations" << std::endl;
                    }
                }
            ));
        }
    }

   // Wait for remaining tasks
    for (auto& fut : futures) {
        fut.wait();
    }
}

int main(int argc, char* argv[]) {
    std::string input_dir, output_dir;
    std::vector<std::size_t> reads;
    std::string pattern_str = ".*_R[12]_001\\.fastq(\\.gz)?";
    std::size_t start = 0, stop = 0, step = 0;
    std::size_t seed = 0;
    double ratio = 0.0;
    bool force = false;
    bool unpaired = false;

    // Define long-named options
    static struct option long_options[] = {
        {"in", required_argument, 0, 'i'},
        {"out", required_argument, 0, 'o'},
        {"regex", required_argument, 0, 'p'},
        {"start", required_argument, 0, 'b'},
        {"stop", required_argument, 0, 'e'},
        {"step", required_argument, 0, 'c'},
        {"seed", required_argument, 0,'s'},
        {"ratio", required_argument, 0, 'r'},
        {"force", no_argument, 0, 'f'},
        {"unpaired", no_argument, 0, 'u'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "i:o:p::b::e::c::s::r::f", long_options, nullptr)) != -1) {
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
            case 'c':
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
            case 'u':
                unpaired = true;
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
    
    if (unpaired) {
        std::vector<std::string> unpaired_files = getUnpairedFiles(input_dir, pattern);
        process_unpaired_files(unpaired_files, output_dir, reads, ratio, force, seed);
    } else {
        //
        std::vector<std::pair<std::string, std::string>> paired_files = getPairedFiles(input_dir, pattern);
        process_paired_files(paired_files, output_dir, reads, ratio, force, seed);
    }
    return 0;
}