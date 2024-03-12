#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <zlib.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>
#include <filesystem>

namespace fs = std::filesystem;

// Helper function to split a string by a delimiter
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::istringstream tokenStream(s);
    std::string token;
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to read k-mers from the first column of a CSV, skipping the header line
std::vector<std::string> readKmersFromCSV(const fs::path& filePath) {
    std::ifstream file(filePath);
    std::vector<std::string> kmers;

    if (!file.is_open()) {
        throw std::runtime_error("Unable to open k-mer file: " + filePath.string());
    }

    std::string line;
    getline(file, line); // Skip the header line

    while (std::getline(file, line)) {
        auto tokens = split(line, ',');
        if (!tokens.empty()) {
            kmers.push_back(tokens.front());
        }
    }

    file.close();
    return kmers;
}

// Function to remove the extension from a filename
std::string stripExtension(const std::string& filename) {
    size_t lastdot = filename.find_last_of('.');
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
}

// Global variables for worker threads
std::mutex queueMutex;
std::condition_variable cvQueue;
std::queue<std::string> workQueue;
std::atomic<bool> done(false);

// Worker thread function to count k-mers
void workerThread(const std::vector<std::string>& kmers, std::map<std::string, std::atomic<int>>& kmerCounts) {
    while (true) {
        std::string sequence;
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            cvQueue.wait(lock, [] { return !workQueue.empty() || done.load(); });
            if (done && workQueue.empty()) break; // No more work to be done
            sequence = std::move(workQueue.front());
            workQueue.pop();
        }
        
        // Count k-mers in the given sequence
        for (const auto& kmer : kmers) {
            size_t pos = 0;
            while ((pos = sequence.find(kmer, pos)) != std::string::npos) {
                kmerCounts[kmer]++;
                pos += kmer.length(); // Move past this k-mer
            }
        }
    }
}

int main(int argc, char* argv[]) {
    // Check command-line arguments
    if (argc != 9) {
        std::cerr << "Usage: " << argv[0]
                  << " --kmer_dir <kmer_csv_dir>"
                  << " --fastq_file <fastq_gz_file>"
                  << " --output_dir <output_csv_dir>"
                  << " --threads <number_of_threads>\n";
        return EXIT_FAILURE;
    }

    // Parse command-line arguments
    fs::path kmerDirPath;
    fs::path outputDirPath;
    std::string fastqFilePath;
    int numThreads = 1;

    for (int i = 1; i < argc; i += 2) {
        std::string arg = argv[i];
        if (arg == "--kmer_dir") {
            kmerDirPath = fs::path(argv[i + 1]);
        } else if (arg == "--fastq_file") {
            fastqFilePath = argv[i + 1];
        } else if (arg == "--output_dir") {
            outputDirPath = fs::path(argv[i + 1]);
        } else if (arg == "--threads") {
            numThreads = std::stoi(argv[i + 1]);
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            return EXIT_FAILURE;
        }
    }

    // Verify that the k-mer directory exists and is a directory
    if (!fs::exists(kmerDirPath) || !fs::is_directory(kmerDirPath)) {
        std::cerr << "Provided path for k-mer directory is not valid: " << kmerDirPath << std::endl;
        return EXIT_FAILURE;
    }

    // Create the output directory if it doesn't exist
    if (!fs::exists(outputDirPath)) {
        std::cout << "Creating output directory: "<< outputDirPath << std::endl;
        fs::create_directories(outputDirPath);
    }

    // Read k-mers and process FASTQ for each CSV file in the k-mer directory
    for (const auto& dirEntry : fs::directory_iterator(kmerDirPath)) {
        // Skip if not a CSV file
        if (dirEntry.path().extension() != ".csv") {
            continue;
        }

        std::vector<std::string> kmers = readKmersFromCSV(dirEntry.path());
        std::map<std::string, std::atomic<int>> kmerCounts;

        // Initialize counts to 0
        for (const auto& kmer : kmers) {
            kmerCounts[kmer] = 0;
        }

        // Set the flag done to false before starting threads
        done = false;

        // Start worker threads
        std::vector<std::thread> threads;
        for (int i = 0; i < numThreads; ++i) {
            threads.emplace_back(workerThread, std::cref(kmers), std::ref(kmerCounts));
        }

        // Open the FASTQ file with gzip compression
        gzFile gzfp = gzopen(fastqFilePath.c_str(), "rb");
        if (gzfp == Z_NULL) {
            std::cerr << "Failed to open FASTQ file: " << fastqFilePath << std::endl;
            done = true;
            cvQueue.notify_all();
            for (auto& t : threads) {
                t.join();
            }
            // Skip to the next file
            continue;
        }

        // Process the FASTQ file
        char buffer[1024];
        std::string line;
        int lineCounter = 0;
        while (gzgets(gzfp, buffer, sizeof(buffer))) {
            line = buffer;
            if (++lineCounter % 4 == 2) { // Sequence line
                line.pop_back(); // Remove new line character at the end
                std::lock_guard<std::mutex> lock(queueMutex);
                workQueue.push(line);
                cvQueue.notify_one();
            }
        }
        gzclose(gzfp);

        // Signal worker threads that there is no more data
        done = true;
        cvQueue.notify_all();

        // Wait for worker threads to finish
        for (auto& t : threads) {
            t.join();
        }

        // Write the counts to the output CSV, maintaining order
        fs::path outputFilePath = outputDirPath / (stripExtension(dirEntry.path().filename().string()) + "_counts.csv");
        std::ofstream outputFile(outputFilePath);
        if (!outputFile.is_open()) {
            std::cerr << "Failed to open output file for writing: " << outputFilePath << std::endl;
            continue;
        }

        outputFile << "K-mer,Count\n";
        for (const auto& kmer : kmers) {
            outputFile << kmer << "," << kmerCounts[kmer].load() << "\n";
        }
        outputFile.close();

        // Clear work queue and kmer counts for the next file
        std::queue<std::string> emptyQueue;
        std::swap(workQueue, emptyQueue);
        kmerCounts.clear();
    }

    return EXIT_SUCCESS;
}

