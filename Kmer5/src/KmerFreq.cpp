/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file KmerFreq.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 22 December 2023, 10:00
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
        throw std::out_of_range("void KmerFreq::setFrequency(int frequency): frequency < 0");
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

void KmerFreq::write(std::ostream& outputStream) {
    _kmer.write(outputStream);
    outputStream.write(reinterpret_cast<const char*>(&_frequency), sizeof(_frequency));
}

void KmerFreq::read(std::istream& inputStream) {
    _kmer.read(inputStream);
    inputStream.read(reinterpret_cast<char*>(&_frequency), sizeof(_frequency));
}

std::ostream& operator<<(std::ostream& os, const KmerFreq& kmerFreq) {
    os << kmerFreq._kmer << " " << kmerFreq._frequency;
    return os;
}

std::istream& operator>>(std::istream& is, KmerFreq& kmerFreq) {
    string kmerStr;
    is >> kmerStr;
    Kmer kmer(kmerStr);
    kmerFreq.setKmer(kmer);
    is >> kmerFreq._frequency;
    return is;
}

bool operator>(KmerFreq kmerFreq1, KmerFreq kmerFreq2) {
    if (kmerFreq1._frequency > kmerFreq2._frequency) {
        return true;
    } else if (kmerFreq1._frequency == kmerFreq2._frequency) {
        return kmerFreq1._kmer > kmerFreq2._kmer;
    } else {
        return false;
    }
}

bool operator<(KmerFreq kmerFreq1, KmerFreq kmerFreq2) {
    if (kmerFreq1._frequency < kmerFreq2._frequency) {
        return true;
    } else if (kmerFreq1._frequency == kmerFreq2._frequency) {
        return kmerFreq1._kmer < kmerFreq2._kmer;
    } else {
        return false;
    }
}

bool operator==(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {
    return kmerFreq1 == kmerFreq2;
}
 
bool operator!=(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {
    return kmerFreq1 != kmerFreq2;
}

bool operator>=(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {
    return kmerFreq1 >= kmerFreq2;
}

bool operator<=(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {
    return kmerFreq1 <= kmerFreq2;
}