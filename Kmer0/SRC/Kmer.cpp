/*
 * Metodología de la Programación: Kmer0
 * Curso 2023/2024
 */

/** 
 * @file Kmer.cpp
 * @author Luis Pérez Velasco <luispv05@correo.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * 
 * Created on 24 October 2023, 14:00
 */

#include "Kmer.h"

using namespace std;

Kmer::Kmer(int k ){
    if (k<= 0){
        throw invalid_argument("k debe ser mayor a 0");
    }    
    _text = string(k, MISSING_NUCLEOTIDE);
       
}
   
   
Kmer::Kmer(const string& text){
    if (text.empty()) {
        throw invalid_argument("El texto no puede estar vacio");
    }
    _text = text;
}

   
int Kmer::getK() const{
    return _text.length();
}
   
   
int Kmer::size() const {
    return _text.size();
}
   
   
string Kmer::toString() const{
    return _text;
}

const char& Kmer::at(int index) const{
    if ((index < 0 )||(index > _text.size()-1)) {
            throw out_of_range("Indice fuera del rango del k-mero");
    }
    return _text[index];
}

   
char& Kmer::at(int index){
    if ((index < 0 )||(index > _text.size()-1)) {
        throw out_of_range("indice fuera del rango del k-mero");
    }
    return _text[index];
}
 
   
void Kmer::normalize(const string& validNucleotides){
    for(char& nucleotides : _text ){
        if(!IsValidNucleotide(nucleotides, validNucleotides)){
            nucleotides = MISSING_NUCLEOTIDE;
        }
        else {
            nucleotides = ToUpper(nucleotides);
        }
    }
}
   
   
Kmer Kmer::complementary(const string& nucleotides,const string& complementaryNucleotides) const {
    if (nucleotides.size() != complementaryNucleotides.size()) {
        throw invalid_argument("El tamaño de los nucleotidos y de los nucleotidos complementarios no coinciden");
    }
    Kmer kmerComplementario(_text);
    for (char& nucleotides : kmerComplementario._text) {
        auto pos = nucleotides.find(ToUpper(nucleotides));
        if (pos != string npos) {
            nucleotides = Complemetarynucleotides[pos];
        }
    }
    return kmerComplementario;
}
   

bool IsValidNucleotide(char nucleotide, const string& validNucleotides){
    bool es_valido = false;
   
    int i = 0;
    while ((i < (int)validNucleotides.length()) && !es_valido) {
       
        if (nucleotide == validNucleotides[i]){
            es_valido = true;
        }
        i++;
       
    }
   
    return es_valido;
}


void ToLower(Kmer& kmer){
   
   
    for (int i = 0; i < kmer.size(); i++) {
       
        if ( (isalpha(kmer.at(i))) && (kmer.at(i) >= 'A') && ( kmer.at(i) <= 'Z' ) ) {
           
            kmer.at(i) += 'a' -'A';
        }
    }
}


void ToUpper(Kmer& kmer){
   
   
    for (int i = 0; i < kmer.size(); i++) {
       
        if ( (isalpha(kmer.at(i))) && (kmer.at(i) >= 'a') && ( kmer.at(i) <= 'z' ) ) {
           
            kmer.at(i) += 'A' -'a';
        }
    }
}

