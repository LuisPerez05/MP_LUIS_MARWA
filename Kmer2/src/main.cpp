/*
 * Metodología de la Programación: Kmer2
 * Curso 2023/2024
 */

/**
 * @file main.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 *
 * Created on 31 October de 2023, 14:30
 */

#include <iostream>

#include "Profile.h"

using namespace std;


/**
 * Shows help about the use of this program in the given output stream
 * @param outputStream The output stream where the help will be shown (for example,
 * cout, cerr, etc) 
 */
void showEnglishHelp(ostream& outputStream) {
    outputStream << "ERROR in Kmer2 parameters" << endl;
    outputStream << "Run with the following parameters:" << endl;
    outputStream << "kmer2 <outputFile.prf> <file1.prf> [<file2.prf> ... <filen.prf>] " << endl;
    outputStream << endl;
    outputStream << "Parameters:" << endl;
    outputStream << "<outputFile.prf>: output file where the join profile"
            << " will be saved" << endl;
    outputStream << "<file1.prf> [<file2.prf> ... <filen.prf>]: input profile files" << endl;  
    outputStream << endl;
    outputStream << "This program obtains a join profile (<outputFile.prf>)" 
            << " using the provided input profiles <file1.prf> [<file2.prf> ... <filen.prf>]" << endl;
    outputStream << endl;
}

/**
 * This program reads and normalizes an undefined number of Profile objects 
 * from the input files passed as parameters to main(). These normalized Profile 
 * objects are used to obtain the union of them, the join Profile. The join 
 * Profile is then zipped (pairs with frequency less or equals to zero are 
 * deleted) and sorted by decreasing order of frequency. In that order,
 * if there is any tie in frequencies, then alphabetical order of kmers is
 * applied. Finally, the resulting sorted Profile is saved in the 
 * output file. 
 * The program must receive at least an input file and an output file. 
 * The output Profile will have as profile identifier, the one of the first 
 * input file (<file1.prf>). 
 * If an input file <file*.prf> has a profile identifier different from the one
 * of the first file (<file1.prf>), then it will not be included in the union. 
 * 
 * Running sintax:
 * > kmer2 <outputFile.prf> <file1.prf> [<file2.prf> ... <filen.prf>]  
 * 
 * Running example:
 * > kmer2 /tmp/outputFile.prf ../Genomes/5pairsRNA.prf ../Genomes/6pairsRNA.prf
 * 
 * > cat /tmp/outputFile.prf
MP-KMER-T-1.0
red bug
7
C_ 6
AG 4
CC 4
CG 4
GC 2
G_ 2
_G 2
 */
int main(int argc, char* argv[]) {
    // This string contains the list of nucleotides that are considered as
    // valid within a genetic sequence. The rest of characters are considered as
    // unknown nucleotides 
    const string VALID_NUCLEOTIDES = "ACGT";
    

    // Check if the number of running arguments is correct, otherwise call to
    // showEnglishHelp(cerr) and end main()

    // Load and normalize the first input Profile file

    // Use a loop to load, normalize and obtain the join with the rest of Profiles

    // Zip the resulting Profile

    // Sort the zipped Profile

    // Save the final Profile to the output file
    
    #include <iostream>
#include <fstream>
#include <vector>
#include "Profile.h"

using namespace std;

/**
 * Shows help about the use of this program in the given output stream
 * @param outputStream The output stream where the help will be shown (for example,
 * cout, cerr, etc) 
 */
void showEnglishHelp(ostream& outputStream) {
    outputStream << "ERROR in Kmer2 parameters" << endl;
    outputStream << "Run with the following parameters:" << endl;
    outputStream << "kmer2 <outputFile.prf> <file1.prf> [<file2.prf> ... <filen.prf>] " << endl;
    outputStream << endl;
    outputStream << "Parameters:" << endl;
    outputStream << "<outputFile.prf>: output file where the join profile"
            << " will be saved" << endl;
    outputStream << "<file1.prf> [<file2.prf> ... <filen.prf>]: input profile files" << endl;  
    outputStream << endl;
    outputStream << "This program obtains a join profile (<outputFile.prf>)" 
            << " using the provided input profiles <file1.prf> [<file2.prf> ... <filen.prf>]" << endl;
    outputStream << endl;
}

int main(int argc, char* argv[]) {
    // Check if the number of running arguments is correct, otherwise call to showEnglishHelp(cerr) and end main()
    if (argc < 3) {
        showEnglishHelp(cerr);
        return 1;
    }

    // Load and normalize the first input Profile file
    const string VALID_NUCLEOTIDES = "ACGT";
    ifstream firstInputFile(argv[2]);
    if (!firstInputFile) {
        cerr << "Error: Unable to open the first input Profile file: " << argv[2] << endl;
        return 1;
    }

    Profile joinProfile;
    firstInputFile >> joinProfile; // Assuming operator>> is overloaded in Profile class for reading from file
    firstInputFile.close();

    // Normalize the first Profile
    // Assuming you have a normalize() method in the Profile class
    joinProfile.normalize(VALID_NUCLEOTIDES);

    // Use a loop to load, normalize, and obtain the join with the rest of Profiles
    for (int i = 3; i < argc; ++i) {
        ifstream inputFile(argv[i]);
        if (!inputFile) {
            cerr << "Warning: Unable to open input file: " << argv[i] << ". Skipping." << endl;
            continue;
        }

        Profile tempProfile;
        inputFile >> tempProfile; // Assuming operator>> is overloaded in Profile class for reading from file
        inputFile.close();

        // Normalize the loaded Profile
        tempProfile.normalize(VALID_NUCLEOTIDES);

        // Obtain the join with the loaded Profile
        joinProfile.join(tempProfile);
    }

    // Zip the resulting Profile (delete pairs with frequency less or equals to zero)
    joinProfile.zip(true, 0);

    // Sort the zipped Profile by decreasing order of frequency and then alphabetical order of kmers (if there is a tie)
    joinProfile.sort();

    // Save the final Profile to the output file
    ofstream outputFile(argv[1]);
    if (!outputFile) {
        cerr << "Error: Unable to open the output file: " << argv[1] << endl;
        return 1;
    }
    outputFile << joinProfile; // Assuming operator<< is overloaded in Profile class for writing to file
    outputFile.close();

    return 0;
}
