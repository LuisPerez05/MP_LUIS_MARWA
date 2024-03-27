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
    if(nElements>dim){
        nElements = dim;
    }
    else if(nElements < 0){
        nElements = 0;
    }
    
    for(int i=0; i<nElements; i++){
        array[i].getKmer();
    }
}

void PrintArrayKmerFreq(KmerFreq array[], int nElements){
    std::cout<<"Número de elementos: "<< nElements <<std::endl;
    for(int i=0; i<nElements; i++){
        array[i].toString();
    }
}

void SwapElementsArrayKmerFreq(KmerFreq array[], int nElements, int first, int second){
    if(first < 0 || first >=nElements || second < 0 || second >= nElements){
        throw std::out_of_range("Index fuera del rango");
    }
    std::swap(array[first], array[second]);
}

int FindKmerInArrayKmerFreq(KmerFreq array[], Kmer kmer, int initialPos, int finalPos){
    for(int i=initialPos; i<=finalPos; i++){
        if(array[i].getKmer()== kmer){
            return i;
        }
    }
    return -1;
}

void SortArrayKmerFreq(KmerFreq array[], int nElements){
    for (int i = 0; i < nElements - 1; i++) {
        for (int j = i + 1; j < nElements; j++) {
            if (array[i].getFrequency() < array[j].getFrequency()) {
                std::swap(array[i], array[j]);
            }
        }
    }
}

void DeletePosArrayKmerFreq(KmerFreq array[], int nElements, int pos){
    if (pos < 0 || pos >= nElements)
        throw std::out_of_range("Index fuera de rango");

    for (int i = pos; i < nElements - 1; i++) {
        array[i] = array[i + 1];
    }
    nElements--;
}


void ZipArrayKmerFreq(KmerFreq array[], int nElements, bool deleteMissing=false, int lowerBound=0){
    for (int i = 0; i < nElements; i++) {
        if () {
            DeletePosArrayKmerFreq(array, nElements, i);
            i--;
        } else if (array[i].getFrequency() <= ) {
            DeletePosArrayKmerFreq(array, nElements, i);
            i--;
        }
    }
}