/**
 * Copyright (c) 2020 joverbeck8@gmail.com
 *
 * Author: Josh Overbeck
 * Date Created: 5/26/20
 * Description: This is a simple program for processing the genome in a fasta file.
 *              It will take in a path to a fasta file.  It will then provide the 
 *              genome count for A, T, G, C, & N.
 *
 * Compile: g++ -Wall -O3 -std=c++14 -pthread -o bio-util Bio_Util.cpp
 */

#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <exception>
#include <vector>
#include <algorithm>
#include <thread>
#include <mutex>  // Lock guard

using BigInt = uint64_t;  // Change here if experiencing overflow
using ThrdVec = std::vector<std::thread>;
using BigIVec = std::vector<BigInt>;

// Globals to have on the heap
int numThreads;
std::mutex mute;
std::fstream ofile;

/**
 * This is a struct to help manage getting descriptions of genomes from the 
 * file.
 */
struct Description {
    std::string desc;  // To hold the description of the genome
    BigInt ending;     // To hold the index that the description ended at
};  // End of the 'Description' struct

/**
 * This is a helper function that will prompt out the usage to the user.
 */
void usage() {
    std::cerr << "Usage: ./<EXECUTABLE> <PATH_TO_FASTA_FILE> <NUM_THREADS>\n";
}  // End of the 'usage' function

/**
 * This is the function that will print out the stats collected from the file.
 *
 * @param desc The description of the genome from the file.
 * @param G The count of G nucleotides in the genome.
 * @param C The count of C nucleotides in the genome.
 * @param T The count of T nucleotides in the genome.
 * @param A The count of A nucleotides in the genome.
 * @param N The count of N necleotides in the genome.
 * @param total The total count of nucleotides in the genome.
 */
void printStats(const std::string& desc, size_t G, size_t C, size_t A, 
                                        size_t T, size_t N, size_t total) {
//    {  // Critical section
//    std::lock_guard<std::mutex> lock(mute);
//    ofile << '\n' << desc << '\n' << '\n';
//    ofile << "G: " << G << '\n';
//    ofile << "C: " << C << '\n';
//    ofile << "A: " << A << '\n';
//    ofile << "T: " << T << '\n';
//    ofile << "N: " << N << '\n';
//    ofile << std::setw(35) << std::setfill('-');
//    ofile << "\nTotal: " << total << '\n';
//    }  // End of the critical section

    std::stringstream table;
    table << "\n" << desc << "\n\n";
    table << "G: " << G << "\n";
    table << "C: " << C << "\n";
    table << "A: " << A << "\n";
    table << "T: " << T << "\n";
    table << "N: " << N << "\n";
    table << "-----------------------------------";
    table << "\nTotal: " << total << "\n";
    {  // Critical section
    std::lock_guard<std::mutex> lock(mute);
    ofile << table.str();
    }
}  // End of the 'printStats' function

/**
 * This is the function that will get the description of a genome from the file
 * that is opened.
 *
 * @param mem The file that is mmaped.
 * @param start The index of where the description starts.
 * @param size The size of the file.
 * @returns The description struct that contains the name from the file.  And
 *          the index that the description ends at.
 */
Description getDescription(const char* mem, BigInt start, BigInt size) {
    Description des;
    std::string name = "";
    for (BigInt i = start; i < size; i++) {
        if (mem[i] == '\n') {
            des.ending = i;
            break;
        }
        name += mem[i];
    }
    des.desc = name;
    return des;
}  // End of the 'getDescription' function

/**
 * This is the function that will read each nucleotide from the file and 
 * collect thier counts.  Then it will invoke the function that prints out the 
 * counts.  It will be run as a task so that it can be parallelized.
 *
 * @param desc The description of the genome from the file. 
 * @param start The starting index after the description of the genome.
 * @param end The ending index of the individual genome.
 * @param mem The char array of the file.
 */
void collectCounts(std::string desc, BigInt start, BigInt end, const char* mem) {
    // Declare size_ts to hold counts
    size_t G = 0; size_t C = 0; size_t A = 0; 
    size_t T = 0; size_t N = 0; size_t total = 0;
    
    // Loop through nucleotides in the file
    char nucleotide;
    for (BigInt i = start; i < end; i++) {
        nucleotide = mem[i];
        switch(nucleotide) {
            case 'G':
                G ++;
                total++;
                break;
            case 'C':
                C ++;
                total++;
                break;
            case 'A':
                A ++;
                total++;
                break;
            case 'T':
                T ++;
                total++;
                break;
            case 'N':
                N ++;
                total++;
                break;
        }  // End of the switch/case block
    }  // End of the loop block
    printStats(desc, G, C, A, T, N, total);
}  // End of the 'collectCounts' function

/**
 * This is a helper function that will manage threads for counting nucleotides 
 * in each genome.
 *
 * @param indicies The vector containing each index each genome falls in.
 * @param mem The char array that contains the FASTA file.
 * @param size The size of the the char array.
 */
void stageCollections(BigIVec& indicies, const char* mem, BigInt size) {
    std::cout << "Counting nucleotides...\n";
    ThrdVec threads;
    int threadsUsed;
    for (BigInt i = 0; i < (indicies.size() - 1); (i += threadsUsed)) {
        threadsUsed = 0;
        for (int t = 0; t < numThreads; t++) {
            // Calculate  the index
            BigInt index = i + t;
            // Prevent overflow past indicies vector's size
            if (index < (indicies.size() - 1)) {
                // Get the desciption
                Description des = getDescription(mem, indicies[index], size);
                // Collect counts of the nucleotides
                threads.push_back(std::thread(collectCounts, des.desc, 
                            des.ending, indicies[index + 1], mem));
                threadsUsed++;
            }
        }
        for (auto &th : threads) {
            th.join();
        }
        threads.clear();
    }
    std::cout << "Done counting nucleotides...\n";
}  // End of the 'stageCollections' function

/**
 * This is a helper function for the for each thread to execute.  It will 
 * check for indecies where a new genome starts.  Then it will append  them 
 * to the vector.
 *
 * @param indicies The vector that contains all of the indicies of 
 *                 beginning/ending of each genome.
 * @param start The beginning index of the file for this thread to begin at.
 * @param end The ending indes of the file for this thread to end at.
 * @param mem The char array that holds the contents of the FASTA file.
 */
void processChunks(BigIVec& indicies, BigInt start, BigInt end, const char* mem) {
    for (BigInt i = start; i < end; i++) {
        if (mem[i] == '>') {
            {  // Critical section
            std::lock_guard<std::mutex> lock(mute);
            indicies.push_back(i);
            }
        }
    }
}  // End of the 'getIndicies' function

/**
 * This is a helper method to get all the start end endpoints of each 
 * genome in the file.
 *
 * @param mem The FASTA file.
 * @param size The size of the file.
 *
 * @returns A vector that contains all of the indicies.
 */
void getIndicies(BigIVec& indicies, const char* mem, BigInt size) {
    // Create a vector to hold the threads
    ThrdVec threads;
    // Find out the size of the chunk of the file each thread gets
    BigInt chunkSize = size / numThreads;

    // Launch threads to process chunks of the file
    for (int chunk = 0; chunk < numThreads; chunk++) {
        BigInt start = chunk * chunkSize;
        BigInt end   = (chunk == (numThreads -1)) ? size : (start + chunkSize);
        threads.push_back(std::thread(processChunks, std::ref(indicies), 
                    start, end, mem));
    }

    // Block until all the threads are finished
    for (auto& t : threads) {
        t.join();
    }

    threads.clear();
    
    // Include the end of the file
    indicies.push_back(size);

    // Sort the indicies in ascending order 
    std::sort(indicies.begin(), indicies.end(), [](BigInt a, BigInt b)
            {return a < b;});

//    // TMP
//    for (auto i : indicies) {
//        std::cout << i << '\n';
//    }
}  // End of the 'getIndicies' function

/**
 * This is the function that will open the file.  Then it will get the size of 
 * the file and put it on the heap as a char array.  Then it will grab the 
 * description from the file.  And then invoke the function that will get the 
 * counts for the nucleotides.
 *
 * @param file The path to the fasta file.
 */
void readFile(std::string file) {
    // Declaring for the file and file descriptor
    const char* mem;
    int fd;

    // Necessary for getting the size
    struct stat sb;

    // Open the file and get the stats
    fd = open(file.c_str(), O_RDONLY);
    fstat(fd, &sb);

    // MMap the file and make sure it was successful
    mem = static_cast<char*>(mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
    if (mem == MAP_FAILED) {
        std::cerr << "mmap failed" << std::endl;
        exit(-1);
    }
    
    std::cout << "Pre-processing..." << std::endl;

    // Get all of the indicies of each  genome in the file
    BigIVec indicies; 
    getIndicies(indicies, mem, sb.st_size);

    std::cout << "Done pre-processing..." << std::endl;

    // Stage the threads for nucleotide counting
    stageCollections(indicies, mem, sb.st_size);

    // Close the file
    close(fd);
}  // End of the 'readFile' function

/**
 * The main function.
 */
int main (int argc, char** argv) {
    // Make sure that the user enter the right number of args
    if (argc != 3) {
        // Prompt the usage
        usage();
    } else {
        // Get the file path supplied by the user
        std::string filePath;
        filePath = argv[1];
        try {
            // Get the number of threads to use
            numThreads = std::stoi(argv[2]);
            // Open a file to put the output in
            ofile.open("out.txt", std::ios::out);
            // Invoke the function that will read the file
            readFile(filePath);
            std::cout << "Output is stored in file named out.txt" << std::endl;
            ofile.close();
        } catch (std::exception& e) {
            // If things go wrong, prompt the usage 
            usage();
            // And let the user know the program crashed
            std::cerr << "Something went wrong..." << std::endl;
            std::cerr << e.what() << std::endl;
        }
    }
    return 0;
}  // End of the 'main' function


// END OF FILE
