/*
 * Metodología de la Programación: Kmer0
 * Curso 2023/2024
 */

/* 
 * File:   main.cpp
 * @author Luis Pérez Velasco <luispv05@correo.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * 
 * Created on 24 October 2023, 13:58
 */

#include <iostream>
#include <string>

#include "Kmer.h"

using namespace std;

/**
 * This program first reads from the standard input an integer k (length of Kmer)
 * and a string with a genetic sequence. Then, it obtains from the genetic 
 * sequence, the list of kmers (of length k) and saves them in the array kmers. 
 * Then, the kmers are normalized. After that, the complementary kmers, 
 * converted to lowercase, are saved in the array complementaryKmers. Finally 
 * the kmers in the arrays kmers and complementaryKmers are shown in the 
 * standard output.
 * See the next example:
 * 
 * Running example:
 * > kmer0 < data/easyDNA5_missing.k0in
6
GCGCC<-->cgcgg
CGCCC<-->gcggg
GCCC_<-->cggg_
CCC_G<-->ggg_c
CC_G_<-->gg_c_
C_G_G<-->g_c_c
 */
int main(int argc, char* argv[]) {
    // This string contains the list of nucleotides that are considered as
    // valid within a genetic sequence. The rest of characters are considered as
    // unknown nucleotides 
    const string VALID_NUCLEOTIDES = "ACGT";
    
    // This string contains the list of complementary nucleotides for each
    // nucleotide in validNucleotides
    const string COMPLEMENTARY_NUCLEOTIDES = "TGCA";

    // This is a constant with the dimension of the array kmers
    const int DIM_ARRAY_KMERS = 100;
    
    // This is the array where the kmers of the input genetic sequence will be
    // saved
    Kmer kmers[DIM_ARRAY_KMERS];
    
    // This is the array where the complementary kmers will be
    // saved
    Kmer complementaryKmers[DIM_ARRAY_KMERS];
    
    // Read K (integer) and a string with the input nucleotides list
    
    // Obtain the kmers: find the kmers in the input string and put them in an array of Kmer
    
    // Normalize each Kmer in the array
    
    // Obtain the complementary kmers and turn them into lowercase

    // Show the list of kmers and complementary kmers as in the example
    
    int k, secuencia_genetica;
    
    cout <<"Introduzca k: ";
    cin >> k;
    cout <<"Introduzca la secuencia: ";
    cin >>secuencia_genetica;
    
    int num_fragmentos;
    
    if(secuencia_genetica.length()>=k){
        num_fragmentos = secuencia_genetica.length() - k + 1;
        
        if(num_fragmentos > DIM_ARRAY_KMERS){
            num_fragmentos = DIM_ARRAY_KMERS; 
        } else {
            num_fragmentos = 0;
        }
    }
    
     cout<<num_fragmentos<<endl;

    for (int i = 0; i < num_fragmentos  ; i++) {
       }
     
    // Obtener el Kmer desde la secuencia genética
    string secuencia_kmer = "";
    for (int j = 0; j < k; j++) {
            if (i + j < secuencia_genetica.length()) {
                char nucleotido = secuencia_genetica[i + j];
                if (esNucleotidoValido(nucleotido)) {
                    secuencia_kmer += nucleotido;
                } else {
                // Tratar con nucleótidos no válidos (si es necesario)
                }
            } else {
            // Tratar con el final de la secuencia genética (si es necesario)
        }
    }

    // Almacenar el Kmer en el arreglo kmers
    Kmer kmer(secuencia_kmer);
    kmers[i] = kmer;

    
    complementaryKmers[i] = kmers[i].complementary(VALID_NUCLEOTIDES,
        COMPLEMENTARY_NUCLEOTIDES);
        ToLower(complementaryKmers[i]);
    
    cout << kmers[i].toString() <<"<-->"<< complementaryKmers[i].toString();
    cout << endl;
    
        
    return 0;
}
