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
#include<string>
#include "Kmer.h"

Kmer::Kmer(int k){
    if (k<= 0){
        throw std::invalid_argument("k debe ser mayor a 0");
    }    
    std::string text(k, MISSING_NUCLEOTIDE);
    Kmer::_text = text;
       
}
   
   
Kmer::Kmer(const std::string& text){
    if (text.empty()) {
        throw std::invalid_argument("El texto no puede estar vacio");
    }
    Kmer::_text = text;
}

   
int Kmer::getK() const{
    return ((Kmer::_text).size());
}
   
   
int Kmer::size() const {
    return ((Kmer::_text).size());
}
   
   
std::string Kmer::toString() const{
    return (Kmer::_text);
}

const char& Kmer::at(int index) const{
    if(index>(Kmer::_text.size()-1)){
        throw std::out_of_range("Indice fuera del rango del k-mero");
    }
    return (Kmer::_text.at(index));
}

   
char& Kmer::at(int index){
    if(index>(Kmer::_text.size()-1)){
        throw std::out_of_range("indice fuera del rango del k-mero");
    }
    return (Kmer::_text.at(index));
}
 
   
void Kmer::normalize(const std::string& validNucleotides){
    int k = Kmer::_text.size();
    std::string solution;
    for(int n = 0; n < k; n++){
        char c_act = Kmer::_text.at(n);
        if(IsValidNucleotide(c_act, validNucleotides)){
        }else{
            bool minus = false;
            int pos = 0;
            int max = validNucleotides.size();
            while(!minus && (pos < max)){
                if(c_act == tolower(validNucleotides.at(pos))) minus = true;
                pos++;
            }
            if(minus){
                c_act = toupper(c_act);
            }else{
                c_act = MISSING_NUCLEOTIDE;
            }
        }
        solution = solution + c_act;
    }
    Kmer::_text = solution;
}
   
   
Kmer Kmer::complementary(const std::string& nucleotides,const std::string& complementaryNucleotides) const {
    if (nucleotides.size() != complementaryNucleotides.size()) {
        throw std::invalid_argument("El tamaño de los nucleotidos y de los nucleotidos complementarios no coinciden");
    }
    int k = Kmer::_text.size();
    std::string solution;
    
    for(int i = 0; i < k; i++){
        char c_act = Kmer::_text.at(i);
        int max = nucleotides.size();
        bool same = false;
        int pos = 0;
        while ((!same) && (pos < max)){
            if(c_act == nucleotides.at(pos)){
                same = true;
            }
            pos++;
        }
        if(same){
            pos--;
            c_act = complementaryNucleotides.at(pos);
        }
        solution = solution + c_act;
    }
    return(Kmer(solution));
}
   

bool IsValidNucleotide(char nucleotide, const std::string& validNucleotides){
    bool es_valido = false;
    int k = validNucleotides.size();
    int i = 0;
    while((i < k) && (!es_valido)) {
       if(nucleotide == validNucleotides.at(i)){
            es_valido = true;
        }
        i++;
    }
    return(es_valido);
}


void ToLower(Kmer& kmer){
    int k = kmer.size();
    for(int i = 0; i < k; i++){
        kmer.at(i) = tolower(kmer.at(i));
    }
}


void ToUpper(Kmer& kmer){
    int k = kmer.size();
    for(int i = 0; i < k; i++){
        kmer.at(i) = toupper(kmer.at(i));
    }   
}

