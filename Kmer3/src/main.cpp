/*
 * Metodología de la Programación: Kmer3
 * Curso 2023/2024
 */

/**
 * @file main.cpp
 * @author Luis Pérez Velasco <luispv05@correo.ugr.es>
 * @author Marwa Dris Azhir <marwadrisazhir@correo.ugr.es>
 *
 * Created on 16 November 2023, 14:15
 */

#include <iostream>
#include <cstring>
#include "Profile.h"

using namespace std;



/**
 * Shows help about the use of this program in the given output stream
 * @param outputStream The output stream where the help will be shown (for example,
 * cout, cerr, etc) 
 */
void showEnglishHelp(ostream& outputStream) {
    outputStream << "ERROR in Kmer3 parameters" << endl;
    outputStream << "Run with the following parameters:" << endl;
    outputStream << "kmer3 [-t min|max] <file1.prf> <file2.prf> [ ... <filen.prf>]" << endl;
    outputStream << endl;
    outputStream << "Parameters:" << endl;
    outputStream << "-t min | -t max: search for minimun distances or maximum distances (-t min by default)" << endl;
    outputStream << "<file1.prf>: source profile file for computing distances" << endl;
    outputStream << "<file2.prf> [ ... <filen.prf>]: target profile files for computing distances" << endl;  
    outputStream << endl;
    outputStream << "This program computes the distance from profile <file1.prf> to the rest" << endl;
    outputStream << endl;
}


/**
 * This program reads an undefined number of Profile objects from the set of 
 * files passed as parameters to main(). All the Profiles object, except the 
 * first one, must be stored in a dynamic array of Profile objects. Then, 
 * for each Profile in the dynamic array, this program prints to the 
 * standard output the name of the file of that Profile and the distance from 
 * the first Profile to the current Profile. 
 * Finally, the program should print in the standard output, the name of 
 * the file with the Profile with the minimum|maximum  distance to the Profile 
 * of the first file and its profile identifier.
 * 
 * At least, two Profile files are required to run this program.
 * 
 * This program assumes that the profile files are already normalized and 
 * sorted by frequency. This is not checked in this program. Unexpected results
 * will be obtained if those conditions are not met.
 * 
 * Running sintax:
 * > kmer3 [-t min|max] <file1.prf> <file2.prf> [  ... <filen.prf>] 
 * 
 * Running example:
 * > kmer3 ../Genomes/human1.prf ../Genomes/worm1.prf ../Genomes/mouse1.prf 
Distance to ../Genomes/worm1.prf: 0.330618
Distance to ../Genomes/mouse1.prf: 0.224901
Nearest profile file: ../Genomes/mouse1.prf
Identifier of the nearest profile: mus musculus
 * 
 * Running example:
 * > kmer3 -t max ../Genomes/human1.prf ../Genomes/worm1.prf ../Genomes/mouse1.prf 
Distance to ../Genomes/worm1.prf: 0.330618
Distance to ../Genomes/mouse1.prf: 0.224901
Farthest profile file: ../Genomes/worm1.prf
Identifier of the farthest profile: worm
 */
int main(int argc, char* argv[]) {
    
    int first_profile;
    int pos_max = 0;
    int pos_min = 0;
    bool min = true;
    
    // Process the main() arguments
    if(argc < 3){
        showEnglishHelp(cerr);
        return 1;
    }
    
    if(argv[1][0] == '-'){
        if(strcmp(argv[1], "-t") == 0){
            if(strcmp(argv[2], "max") == 0){
                min = false;
                first_profile = 3;
            }else{
                if(strcmp(argv[2], "min") != 0){
                    showEnglishHelp(cerr);
                    return 1;
                }else first_profile = 3;
            }
        }else{
            showEnglishHelp(cerr);
            return 1;

        }
        
    }else first_profile = 1;
    
    // Allocate a dynamic array of Profiles
    Profile *p;
    p = new Profile[argc - first_profile];
    
    
    // Load the input Profiles
    for (int pos = 0; pos < (argc - first_profile); pos++){
        p[pos].load(argv[pos + first_profile]);
    }
    
    // Calculate and print the distance from the first Profile to the rest
    double *dist;
    dist = new double[argc - first_profile];
    
    for(int pos = first_profile; pos < argc - 1; pos++){
        dist[pos - first_profile] = p[0].getDistance(p[pos - first_profile + 1]);
        
        if(min){
            if(dist[pos - first_profile] < dist[pos_min]){
                pos_min = pos - first_profile;
            }
        }else{
            if(dist[pos - first_profile] > dist[pos_max]){
                pos_max = pos - first_profile;
            }
        }
        
        cout << "Distance to " << argv[pos + 1] << ": " << dist[pos - first_profile] << endl;
    }
    
    // Print name of the file and identifier that takes min|max distance to the first one
    if(min){
        cout << "Nearest profile file: " << argv[pos_min + first_profile + 1] << endl;
        cout << "Identifier of the nearest profile: " << p[pos_min + 1].getProfileId() << endl;
    }else{
        cout << "Farthest profile file: " << argv[pos_max + first_profile + 1] << endl;
        cout << "Identifier of the farthest profile: " << p[pos_max + 1].getProfileId() << endl;
    }
    
    // Deallocate the dynamic array of Profile
    delete[] p;
    p = nullptr;
    
    delete[] dist;
    dist = nullptr;
    

    return 0;
}
