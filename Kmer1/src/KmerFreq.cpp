/*
 * Metodología de la Programación: Kmer1
 * Curso 2023/2024
 */

/** 
 * @file KmerFreq.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 27 de octubre de 2023, 11:03
 */

#include "KmerFreq.h"

// Constructor por defecto de KmerFreq
KmerFreq::KmerFreq() : _kmer(), _frequency(0){}

// Obtener el kmer asociado
Kmer KmerFreq::getKmer(){
    return _kmer;
}

// Obtener la frecuencia del kmer
int KmerFreq::getFrequency(){
    return _frequency;
}

// Establecer el kmer asociado
void KmerFreq::setKmer(Kmer kmer){
    _kmer = kmer;
}

// Establecer la frecuencia del kmer
void KmerFreq::setFrequency(int frequency){
    // Verificar si la frecuencia es negativa
    if(frequency<0){
        throw std::out_of_range("Frecuencia no puede ser negativa");
    }
    _frequency = frequency;
}

// Obtener una representación en cadena del KmerFreq
std::string KmerFreq::toString() const{
    return _kmer.toString() + " " + std::to_string(_frequency);
}