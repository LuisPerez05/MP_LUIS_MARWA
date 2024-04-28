/*
 * Metodología de la Programación: Kmer3
 * Curso 2023/2024
 */

/** 
 * @file Kmer.cpp
 * @author Luis Pérez Velasco <luispv05@correo.ugr.es>
 * @author Marwa Dris Azhir <marwadrisazhir@correo.ugr.es>
 
 * Created on 16 November 2023, 14:15
 */

#include "Kmer.h"

using namespace std;

Kmer::Kmer(int k){
    
    if(k <= 0){
        throw std::invalid_argument("Kmer(int k): k was <=0");
    }
    
    string text(k , MISSING_NUCLEOTIDE);    
    
    Kmer::_text = text;
}

Kmer::Kmer(const string& text){
    
    if(text.empty()){
        throw std::invalid_argument("Kmer(const std::string& text): text was empty");
    }
    
    Kmer::_text = text;

}

int Kmer::getK() const{
    return ((Kmer::_text).size());
}

int Kmer::size() const{
    return ((Kmer::_text).size());
}

string Kmer::toString() const{
    return(Kmer::_text);
}

const char& Kmer::at(int index) const{
    
    if(index > (Kmer::_text.size() - 1)){
        throw std::out_of_range("at(int index) const: index > (k-1)");
    }
    
    return(Kmer::_text.at(index));
}

char& Kmer::at(int index){
    
    if(index > (Kmer::_text.size() - 1)){
        throw std::out_of_range("at(int index): index > (k-1)");
    }
    
    return(Kmer::_text.at(index));
}

void Kmer::toLower(){
    
    int k = Kmer::_text.size();
    
    for(int n = 0; n < k; n++){
        
        Kmer::_text.at(n) = tolower(Kmer::_text.at(n));
            
    }
    
    
}

void Kmer::toUpper(){
    
    int k = Kmer::_text.size();
    
    for(int n = 0; n < k; n++){
        
        Kmer::_text.at(n) = toupper(Kmer::_text.at(n));
            
    }

}

void Kmer::normalize(const string& validNucleotides){
    int k = Kmer::_text.size();
    string solution;
    
    toUpper();
    
    for(int n = 0; n < k; n++){
        char c_act = Kmer::_text.at(n);
        
        
        if(!(IsValidNucleotide(c_act, validNucleotides))){
            c_act = MISSING_NUCLEOTIDE;
        }

        
        solution = solution + c_act;
    }
    
    Kmer::_text = solution;
}

Kmer Kmer::complementary(const string& nucleotides, 
               const string& complementaryNucleotides) const{
    
    if(nucleotides.size() != complementaryNucleotides.size()){
        throw std::invalid_argument("complementary(const std::string& nucleotides," 
               "const std::string& complementaryNucleotides) const: nucleotides.size()"
                "!= complementaryNucleotides.size()");
    }
    
    int k = Kmer::_text.size();
    string solution;
    
    for(int n = 0; n < k; n++){
        char c_act = Kmer::_text.at(n);
        
        int pos = nucleotides.find(c_act);
        
        if(pos != -1){
            c_act = complementaryNucleotides.at(pos);
        }
        

        
        solution = solution + c_act;
        
    }
    
    return(Kmer(solution));
    
}


bool IsValidNucleotide(char nucleotide, const string& validNucleotides){
    bool result = false;
    int k = validNucleotides.size();
    
    int n = 0;
    while( (n < k) && !result){
        
        if(nucleotide == validNucleotides.at(n)){
                result = true;
         }
        n++;
    }
        

    
    
    return(result);
}

void ToLower(Kmer& kmer){
    int k = kmer.size();
    
    for(int n = 0; n < k; n++){
        
        kmer.at(n) = tolower(kmer.at(n));
            
    }
    
}

void ToUpper(Kmer& kmer){
    int k = kmer.size();
    
    for(int n = 0; n < k; n++){
        
        kmer.at(n) = toupper(kmer.at(n));
            
    }
    
}
