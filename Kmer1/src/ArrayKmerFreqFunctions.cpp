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

// Función para normalizar las frecuencias de los kmers en el array
void NormalizeArrayKmerFreq(KmerFreq array[], int& nElements, std::string validNucleotides){ 
    
    // Loop to traverse and normalize each one of the kmers in array
    for(int i=0; i<nElements; i++){
    // Normalize kmer i
       array[i].getKmer().normalize(validNucleotides);
    }
    // Loop to traverse the kmers in array from position 1 to position nElements-1
    for(int i=1; i<nElements; i++){
        // index = Position of array[i].getKmer() in the subarray that begins    
        //         at position 0 and ends at position i-1
        int index = FindKmerInArrayKmerFreq(array, array[i].getKmer(), 0, i-1);    
        // If array[i].getKmer() was found in the the subarray from 0 to i-1
        if(index != -1){
            // Accumulate the frequencies of the kmers at positions
            array[index].setFrequency(array[index].getFrequency() + array[i].getFrequency());   
            //    index and i in the kmer at position index
            // Delete from the array, the kmer at position i 
            DeletePosArrayKmerFreq(array, nElements, i); 
            nElements--;
            i--;
        }
    }
}

// Función para leer los elementos del array de kmers y sus frecuencias
void ReadArrayKmerFreq(KmerFreq array[], int dim, int& nElements){
    // Ajustar el número de elementos si excede la dimensión máxima o es negativo
    if(nElements>dim){
        nElements = dim;
    }
    else if(nElements < 0){
        nElements = 0;
    }
     // Bucle para leer los kmers y sus frecuencias
    for(int i=0; i<nElements; i++){
        int frequency;
        std::string Tkmer;
        std::cin>>Tkmer>>frequency;
        array[i].setKmer(Kmer(Tkmer));
        array[i].setFrequency(frequency);
    }
}

// Función para imprimir los elementos del array de kmers y sus frecuencias
void PrintArrayKmerFreq(KmerFreq array[], int nElements){
    // Bucle para imprimir cada elemento del array
    for(int i=0; i<nElements; i++){
        std::cout<<array[i].toString()<<std::endl;
    }
}

// Función para intercambiar elementos en el array de kmers
void SwapElementsArrayKmerFreq(KmerFreq array[], int nElements, int first, int second){
     // Verificar que los índices estén dentro del rango
    if(first < 0 || first >=nElements || second < 0 || second >= nElements){
        throw std::out_of_range("Index fuera del rango");
    }
    // Intercambiar los elementos en las posiciones first y second
    std::swap(array[first], array[second]);
}

// Función para encontrar un kmer en el array de kmers
int FindKmerInArrayKmerFreq(KmerFreq array[], Kmer kmer, int initialPos, int finalPos){
    // Bucle para buscar el kmer en el rango especificado
    for(int i=initialPos; i<=finalPos; i++){
        // Si se encuentra el kmer, se devuelve su posición
        if(array[i].getKmer().toString()== kmer.toString()){
            return i;
        }
    }
    // Si no se encuentra, se devuelve -1
    return -1;
}

// Función para ordenar el array de kmers por frecuencia
void SortArrayKmerFreq(KmerFreq array[], int nElements){
    // Algoritmo de ordenamiento de burbuja
    for (int i = 0; i < nElements - 1; i++) {
        for (int j = i + 1; j < nElements; j++) {
            if (array[i].getFrequency() < array[j].getFrequency()) {
                std::swap(array[i], array[j]);
            }
        }
    }
}

// Función para eliminar un elemento en una posición del array de kmers
void DeletePosArrayKmerFreq(KmerFreq array[], int& nElements, int pos){
    // Verificar que la posición esté dentro del rango
    if (pos < 0 || pos >= nElements)
        throw std::out_of_range("Index fuera de rango");
    // Desplazar los elementos a la izquierda para llenar el espacio
    for (int i = pos; i < nElements - 1; i++) {
        array[i] = array[i + 1];
    }
    // Decrementar el número de elementos
    nElements--;
}

// Función para comprimir el array de kmers, eliminando los que tienen frecuencias menores o iguales a lowerBound
void ZipArrayKmerFreq(KmerFreq array[], int& nElements, bool deleteMissing=false, int lowerBound=0){
    int i=0;
    // Bucle para iterar sobre los elementos del array
    while(i<nElements){
        // Si se desea eliminar los kmers con nucleótidos faltantes y el kmer actual contiene uno, se elimina
        if(deleteMissing && array[i].getKmer().toString().at()){
            for(int j=i; j<nElements; j++){
                if(deleteMissing = true){
                    array[j] = array [j+1]; 
                }
            }
        } 
        // Si la frecuencia del kmer es menor o igual a lowerBound, se elimina
        else if (array[i].getFrequency() <= lowerBound) {
            for(int j = i; j < nElements - 1; j++) {
                array[j] = array[j+1];
            }
            
            nElements--;
        } 
        // Si no se cumple ninguna condición, se pasa al siguiente kmer
        else {
            i++;
        }
    }
}