#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <ctime>

int main()
{
    // Specify the output file path (replace 'ADMIN' with your username)
    std::string outputFilePath = "C:\\Users\\ADMIN\\Desktop\\sample.fasta";//"C:\Users\ADMIN\Desktop\sample.fasta"
    std::ofstream outFile(outputFilePath);
    if (!outFile.is_open())
    {
        std::cerr << "Unable to open output file: " << outputFilePath << std::endl;
        return 1;
    }

    // Set up random number generation
    std::string bases = "ACGT"; // Possible DNA bases
    std::default_random_engine generator(static_cast<unsigned int>(time(0)));
    std::uniform_int_distribution<int> lengthDist(270, 450); // Random length between 50 and 150
    std::uniform_int_distribution<int> baseDist(0, 3);      // For selecting A, C, G, T

    // Generate 300 sequences
    for (int i = 1; i <= 600; ++i)
    {
        outFile << ">sequence" << i << std::endl; // Write identifier
        int length = lengthDist(generator);       // Random length for this sequence
        for (int j = 0; j < length; ++j)
        {
            outFile << bases[baseDist(generator)]; // Write random base
        }
        outFile << std::endl; // New line after sequence
    }

    outFile.close();
    std::cout << "Generated 600 sequences in " << outputFilePath << std::endl;
    return 0;
}