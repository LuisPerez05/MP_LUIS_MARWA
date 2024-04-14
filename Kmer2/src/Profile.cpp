/*
 * Metodología de la Programación: Kmer2
 * Curso 2023/2024
 */

/** 
 * @file Profile.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 29 January 2023, 11:00
 */


#include "Profile.h"

using namespace std;

const string Profile::MAGIC_STRING_T="MP-KMER-T-1.0";

Profile::Profile() { 
    _profileId = "desconocido"; 
    _size = 0; 
} 

Profile::Profile(int size) { 
    _profileId = "desconocido"; 
    if (size < 0 || size > DIM_VECTOR_KMER_FREQ) { 
        throw std::out_of_range("Tamaño no válido por el vecotr de Kmers en el profile"); 
    } 
    _size = size; 
    for (int i = 0; i < _size; ++i) { 
        _vectorKmerFreq[i] = KmerFreq(Kmer::MISSING_NUCLEOTIDE, 0); 
    } 
} 

std::string Profile::getProfileId() { 
    return _profileId; 
} 

void Profile::setProfileId(std::string id) { 
    _profileId = id; 
} 

void Profile::deletePos(int pos) { 
    if (pos < 0 || pos >= _size) { 
        throw std::out_of_range("Posición inválida"); 
    } 
    for (int i = pos; i < _size - 1; ++i) { 
        _vectorKmerFreq[i] = _vectorKmerFreq[i + 1]; 
    } 
    _size--; 
} 

void Profile::zip(bool deleteMissing, int lowerBound) { 
    int i = 0; 
    while (i < _size) { 
        if ((deleteMissing && _vectorKmerFreq[i].getKmer().containsMissingNt()) || _vectorKmerFreq[i].getFreq() <= lowerBound) { 
            deletePos(i); 
        } else { 
            i++; 
        } 
    } 
} 

void Profile::join(Profile profile) { 
    int initialSize = _size; 
    int newSize = initialSize + profile._size; 
    if (newSize > DIM_VECTOR_KMER_FREQ) { 
        throw std::runtime_error("Posición máxima de capacidad excedida del array_vectorKmerFreq"); 
    } 
    for (int i = 0; i < profile._size; ++i) { 
        _vectorKmerFreq[initialSize + i] = profile._vectorKmerFreq[i]; 
    } 
    _size = newSize; 
} 