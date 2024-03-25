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

KmerFreq::KmerFreq() : _kmer(), _frequency(0){}


Kmer KmerFreq::getKmer(){
    return _kmer;
}

int KmerFreq::getFrequency(){
    return _frequency;
}

void KmerFreq::setKmer(Kmer kmer){
    _kmer = kmer;
}

void KmerFreq::setFrequency(int Frequency){
    if(frequency<0){
        throw std::out_of_range("Frecuencia no puede ser negativa");
    }
    _frequency = frequency;
}

std::string KmerFreq::toString() const{
    return (Kmer::_text);
}