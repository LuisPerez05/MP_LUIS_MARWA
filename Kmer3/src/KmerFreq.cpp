/*
 * Metodología de la Programación: Kmer3
 * Curso 2023/2024
 */

/** 
 * @file KmerFreq.cpp
 * @author Luis Pérez Velasco <luispv05@correo.ugr.es>
 * @author Marwa Dris Azhir <marwadrisazhir@correo.ugr.es>
 * 
 * Created on 16 November 2023, 14:15
 */


#include "KmerFreq.h"

using namespace std;

KmerFreq::KmerFreq(){
    _kmer = Kmer(1);
    _frequency = 0;
}

const Kmer& KmerFreq::getKmer() const{
    return _kmer;

}

int KmerFreq::getFrequency() const{
    return _frequency;
}

void KmerFreq::setKmer(const Kmer& kmer) {
    _kmer = kmer;
}

void KmerFreq::setFrequency(int frequency){
    if(frequency < 0){
        throw std::out_of_range("setFrequency(int frequency): frequency < 0");
    }
    
    _frequency = frequency;
}

string KmerFreq::toString() const{
    string aux = _kmer.toString() + " " + to_string(_frequency);
    return aux;
}

bool equal(const KmerFreq& first, const KmerFreq& second){

    return (first.getKmer().toString() == second.getKmer().toString());
}

bool greaterKmerFreq(const KmerFreq& first, const KmerFreq& second){
    bool sol = false;
         
    if(first.getFrequency() > second.getFrequency()){
        sol = true;
    }else if(first.getFrequency() == second.getFrequency()){
                sol = (first.getKmer().toString() < second.getKmer().toString());
            }
    
    return sol;
}

