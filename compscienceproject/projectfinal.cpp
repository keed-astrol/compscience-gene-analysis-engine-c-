#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <utility>
#include <iomanip>

class DNASequence {
private:
    std::string sequence;
    // Static map for restriction enzyme recognition sequences
    static const std::map<std::string, std::string> restrictionEnzymes;

    // Static map for genetic code (codon to amino acid)
    static const std::map<std::string, char> geneticCode;
    static const std::map<std::string, std::string> proteins;

    // Helper function to validate sequence
    void validateSequence(const std::string& seq) {
        for (char c : seq) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
                throw std::invalid_argument("Invalid nucleotide: " + std::string(1, c));
            }
        }
    }

    // Helper function to compute KMP pi table for pattern matching
    static std::vector<size_t> computePi(const std::string& pattern) {
        size_t m = pattern.length();
        std::vector<size_t> pi(m, 0);
        size_t j = 0;
        for (size_t i = 1; i < m; ++i) {
            while (j > 0 && pattern[i] != pattern[j]) {
                j = pi[j - 1];
            }
            if (pattern[i] == pattern[j]) {
                ++j;
            }
            pi[i] = j;
        }
        return pi;
    }

    // Helper function to convert string to uppercase
    static std::string toUpperCase(const std::string& str) {
        std::string upperStr = str;
        std::transform(upperStr.begin(), upperStr.end(), upperStr.begin(),
                       [](unsigned char c) { return std::toupper(c); });
        return upperStr;
    }

public:
    // Constructor
    DNASequence(const std::string& seq) {
        sequence = toUpperCase(seq);
        validateSequence(sequence);
    }

    // Static method to load multiple sequences from FASTA file
    static std::vector<std::pair<std::string, DNASequence>> loadFromFASTA(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Unable to open file: " + filename);
        }
        std::vector<std::pair<std::string, DNASequence>> sequences;
        std::string currentID;
        std::string currentSeq;
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!currentID.empty()) {
                    sequences.emplace_back(currentID, DNASequence(currentSeq));
                    currentSeq.clear();
                }
                currentID = line.substr(1); // Remove '>'
            } else {
                currentSeq += line;
            }
        }
        if (!currentID.empty()) {
            sequences.emplace_back(currentID, DNASequence(currentSeq));
        }
        file.close();
        // Remove whitespace from sequences
        for (auto& pair : sequences) {
            pair.second.sequence.erase(std::remove_if(pair.second.sequence.begin(), pair.second.sequence.end(), ::isspace), pair.second.sequence.end());
        }
        return sequences;
    }

    // Base counting
    size_t countBase(char base) const {
        char upperBase = std::toupper(base);
        if (upperBase != 'A' && upperBase != 'C' && upperBase != 'G' && 
            upperBase != 'T' && upperBase != 'N') {
            throw std::invalid_argument("Invalid base: " + std::string(1, base));
        }
        return std::count(sequence.begin(), sequence.end(), upperBase);
    }

    // Base frequency
    double frequency(char base) const {
        size_t total = sequence.length();
        if (total == 0) return 0.0;
        return static_cast<double>(countBase(base)) / total;
    }

    // GC content percentage
    double gcContent() const {
        size_t total = sequence.length();
        if (total == 0) return 0.0;
        size_t gcCount = countBase('G') + countBase('C');
        return static_cast<double>(gcCount) / total * 100.0;
    }

    // Pattern matching using KMP algorithm
    std::vector<size_t> findMotif(const std::string& motif) const {
        std::string upperMotif = toUpperCase(motif);
        std::vector<size_t> pi = computePi(upperMotif);
        std::vector<size_t> occurrences;
        size_t n = sequence.length();
        size_t m = upperMotif.length();
        size_t j = 0;
        for (size_t i = 0; i < n; ++i) {
            while (j > 0 && sequence[i] != upperMotif[j]) {
                j = pi[j - 1];
            }
            if (sequence[i] == upperMotif[j]) {
                ++j;
            }
            if (j == m) {
                occurrences.push_back(i - m + 1);
                j = pi[j - 1];
            }
        }
        return occurrences;
    }

    // Get complementary strand
    std::string getComplementary() const {
        std::string comp;
        for (char c : sequence) {
            switch (c) {
                case 'A': comp += 'T'; break;
                case 'T': comp += 'A'; break;
                case 'C': comp += 'G'; break;
                case 'G': comp += 'C'; break;
                case 'N': comp += 'N'; break;
                default: comp += 'N'; break;
            }
        }
        return comp;
    }

    // Get reverse complementary strand
    std::string getReverseComplementary() const {
        std::string comp = getComplementary();
        std::reverse(comp.begin(), comp.end());
        return comp;
    }

    // Find restriction enzyme cut sites
    std::vector<size_t> findRestrictionSites(const std::string& enzyme) const {
        auto it = restrictionEnzymes.find(enzyme);
        if (it == restrictionEnzymes.end()) {
            throw std::invalid_argument("Unknown enzyme: " + enzyme);
        }
        std::string pattern = toUpperCase(it->second);
        return findMotif(pattern);
    }

    // Translate DNA sequence to protein
    std::string translateToProtein() const {
        std::string protein;
        size_t seqLen = sequence.length();
        for (size_t i = 0; i + 2 < seqLen; i += 3) {
            std::string codon = sequence.substr(i, 3);
            auto it = geneticCode.find(codon);
            if (it != geneticCode.end()) {
                char aa = it->second;
                if (aa == '*') break; // Stop codon
                protein += aa;
            } else {
                protein += 'X'; // Unknown codon
            }
        }
        return protein;
    }

    // Accessor for sequence length
    size_t length() const {
        return sequence.length();
    }

    // Accessor for sequence
    std::string getSequence() const {
        return sequence;
    }
};

// Definition of static members
const std::map<std::string, std::string> DNASequence::restrictionEnzymes = {
    {"EcoRI", "GAATTC"},
    {"HindIII", "AAGCTT"},
    {"BamHI", "GGATCC"}
};

const std::map<std::string, char> DNASequence::geneticCode = {
    {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
    {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
    {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'},
    {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
    {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
    {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
    {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
    {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
    {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'},
    {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
    {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
    {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
    {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'},
    {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
    {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
    {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
};
// Note: 'proteins' map is declared but not defined in the original code; assuming it's unused for now
const std::map<std::string, std::string> DNASequence::proteins = {};

// Define the codon map (unused in main functionality but kept as per original code)
static const std::unordered_map<std::string, std::string> codon_map = {
    {"UUU", "Phenylalanine"}, {"UUC", "Phenylalanine"},
    {"UUA", "Leucine"}, {"UUG", "Leucine"},
    {"CUU", "Leucine"}, {"CUC", "Leucine"}, {"CUA", "Leucine"}, {"CUG", "Leucine"},
    {"AUU", "Isoleucine"}, {"AUC", "Isoleucine"}, {"AUA", "Isoleucine"},
    {"AUG", "Methionine (Start Codon)"},
    {"GUU", "Valine"}, {"GUC", "Valine"}, {"GUA", "Valine"}, {"GUG", "Valine"},
    {"UCU", "Serine"}, {"UCC", "Serine"}, {"UCA", "Serine"}, {"UCG", "Serine"},
    {"CCU", "Proline"}, {"CCC", "Proline"}, {"CCA", "Proline"}, {"CCG", "Proline"},
    {"ACU", "Threonine"}, {"ACC", "Threonine"}, {"ACA", "Threonine"}, {"ACG", "Threonine"},
    {"GCU", "Alanine"}, {"GCC", "Alanine"}, {"GCA", "Alanine"}, {"GCG", "Alanine"},
    {"UAU", "Tyrosine"}, {"UAC", "Tyrosine"},
    {"UAA", "Stop Codon"}, {"UAG", "Stop Codon"}, {"UGA", "Stop Codon"},
    {"CAU", "Histidine"}, {"CAC", "Histidine"},
    {"CAA", "Glutamine"}, {"CAG", "Glutamine"},
    {"AAU", "Asparagine"}, {"AAC", "Asparagine"},
    {"AAA", "Lysine"}, {"AAG", "Lysine"},
    {"GAU", "Aspartic Acid"}, {"GAC", "Aspartic Acid"},
    {"GAA", "Glutamic Acid"}, {"GAG", "Glutamic Acid"},
    {"UGU", "Cysteine"}, {"UGC", "Cysteine"},
    {"UGG", "Tryptophan"},
    {"CGU", "Arginine"}, {"CGC", "Arginine"}, {"CGA", "Arginine"}, {"CGG", "Arginine"},
    {"AGU", "Serine"}, {"AGC", "Serine"},
    {"AGA", "Arginine"}, {"AGG", "Arginine"},
    {"GGU", "Glycine"}, {"GGC", "Glycine"}, {"GGA", "Glycine"}, {"GGG", "Glycine"}
};

// Main function with protein composition analysis
int main() {
    try {
        std::string fastaFilePath = "C:\\Users\\ADMIN\\Desktop\\sample.fasta";

        // Load sequences from the FASTA file
        auto sequences = DNASequence::loadFromFASTA(fastaFilePath);

        // Output CSV file path on the desktop
        std::string csvFilePath = "C:\\Users\\ADMIN\\Desktop\\analysis_results.csv";
        std::ofstream csvFile(csvFilePath);
        if (!csvFile.is_open()) {
            throw std::runtime_error("Unable to open output CSV file");
        }

        // Define amino acids and their full names
        const std::vector<char> aminoAcids = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'};
        const std::unordered_map<char, std::string> aaNames = {
            {'A', "Alanine"}, {'C', "Cysteine"}, {'D', "Aspartic Acid"}, {'E', "Glutamic Acid"},
            {'F', "Phenylalanine"}, {'G', "Glycine"}, {'H', "Histidine"}, {'I', "Isoleucine"},
            {'K', "Lysine"}, {'L', "Leucine"}, {'M', "Methionine"}, {'N', "Asparagine"},
            {'P', "Proline"}, {'Q', "Glutamine"}, {'R', "Arginine"}, {'S', "Serine"},
            {'T', "Threonine"}, {'V', "Valine"}, {'W', "Tryptophan"}, {'Y', "Tyrosine"},
            {'X', "Unknown"}
        };

        // Write CSV header
        csvFile << "Sequence ID,Length,GC Content (%),Motif ATGC Count,EcoRI Site Count,Protein Length";
        for (char aa : aminoAcids) {
            csvFile << "," << aaNames.at(aa) << " Count," << aaNames.at(aa) << " Percent";
        }
        csvFile << "\n";

        // Set precision for percentages
        csvFile << std::fixed << std::setprecision(2);

        // Analyze each sequence and write results to CSV
        for (const auto& pair : sequences) {
            const auto& id = pair.first;
            const auto& seq = pair.second;
            size_t length = seq.length();
            double gc = seq.gcContent();
            size_t motifCount = seq.findMotif("ATGC").size();
            size_t ecoRICount = seq.findRestrictionSites("EcoRI").size();
            std::string protein = seq.translateToProtein();
            size_t proteinLength = protein.length();

            // Count amino acids
            std::map<char, size_t> counts;
            for (char aa : protein) {
                counts[aa]++;
            }

            // Write data to CSV
            csvFile << "\"" << id << "\"," << length << "," << gc << "," << motifCount << "," << ecoRICount << "," << proteinLength;
            for (char aa : aminoAcids) {
                size_t count = counts[aa]; // 0 if aa not present
                double percent = (proteinLength > 0) ? (static_cast<double>(count) / proteinLength * 100.0) : 0.0;
                csvFile << "," << count << "," << percent;
            }
            csvFile << "\n";
        }
        csvFile.close();
        std::cout << "Analysis complete. Results saved to " << csvFilePath << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } 
    return 0;
}