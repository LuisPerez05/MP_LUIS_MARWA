/*
 * Metodología de la Programación: Kmer1
 * Curso 2023/2024
 */

/** 
 * @file ArrayKmerFreqFunctions.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 27 October 2023, 12:00
 */


#include "ArrayKmerFreqFunctions.h"


void NormalizeArrayKmerFreq(KmerFreq array[], int nElements, std::string validNucleotides){ 
    
    // Loop to traverse and normalize each one of the kmers in array
    for(int i=0; i<nElements; i++){
        int totalFreq = array[i].getFrequency();
        // Normalize kmer i
        array[i].setFrequency(totalFreq);
    
    // Loop to traverse the kmers in array from position 1 to position nElements-1
        for(int j=i+1; j<nElements; j++){
        // index = Position of array[i].getKmer() in the subarray that begins
            int index = -1;
          //         at position 0 and ends at position i-1
            for(int k=0; k<i; k++){
                if(array[k].getKmer()==array[i].getKmer()){
                    index = k;
                    break;
                }
            }
          // If array[i].getKmer() was found in the the subarray from 0 to i-1
            if(index!= -1){
               // Accumulate the frequencies of the kmers at positions
                totalFreq += array[j].getFrequency();
               //    index and i in the kmer at position index
               // Delete from the array, the kmer at position i 
               
            }
        }
    }
}

void ReadArrayKmerFreq(KmerFreq array[], int dim, int nElements){
    
}

void PrintArrayKmerFreq(KmerFreq array[], int nElements){
    
}

void SwapElementsArrayKmerFreq(KmerFreq array[], int nElements, int first, int second){
    
}

int FindKmerInArrayKmerFreq(KmerFreq array[], Kmer kmer, int initialPos, int finalPos){
    
}

void SortArrayKmerFreq(KmerFreq array[], int nElements){
    
}

void DeletePosArrayKmerFreq(KmerFreq array[], int nElements, int pos){
    
}

void ZipArrayKmerFreq(KmerFreq array[], int nElements, bool deleteMissing=false, int lowerBound=0){
    
}