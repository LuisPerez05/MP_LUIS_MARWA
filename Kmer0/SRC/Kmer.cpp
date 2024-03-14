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

Kmer::Kmer(int k ){
       
        if (k<= 0){
            throw invalid_argument("k debe ser mayor a 0");
        }
       
        // Construir el objeto kmer con k caracteres, todos esblecidos como
        // MISSING_NUCLEOTIDES
        _text = string(k, MISSING_NUCLEOTIDE);
       
}
   
   
Kmer::Kmer(const string& text){
       
        if (text.empty()) {
            throw invalid_argument("El texto no puede estar vacio");
        }
       
        // Copiar el valor del string en el objeto kmer tipo string
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
       
        if ((index < 0 )||(index > size()-1)) {
            throw out_of_range("indice fuera del rango del k-mero");
        }
       
        return _text[index];
    }

   
    char& Kmer::at(int index){
                               
        if ((index < 0 )||(index > size()-1)) {
            throw out_of_range("indice fuera del rango del k-mero");
        }
       
        return _text[index];
       
    }
 
   
    void Kmer::normalize(const string& validNucleotides){
        ToUpper(*this);
       
        for(int i = 0; i < size(); i++){
           
            if(!IsValidNucleotide(_text[i], validNucleotides)){
                _text[i] = MISSING_NUCLEOTIDE;
            }

            else if ((_text[i] >= 'a') && (_text[i] <= 'z')) {

                _text[i] = toupper(_text[i]);

            }
           
        }
       
    }
   
   
    Kmer Kmer::complementary(const string& nucleotides,
         const string& complementaryNucleotides) const {
           
        // Check if the sizes of nucleotides and complementaryNucleotides r
        // the same
        if (nucleotides.length() != complementaryNucleotides.length()) {
            throw std::invalid_argument("Sizes of nucleotides and "
                    " complementaryNucleotides must be the same");
        }

        // Create a new Kmer object to store the complementary sequence
        Kmer complementaryKmer(_text);

        // Iterate through each nucleotide in the Kmer and find its complement
        for (int i = 0; i < (int)_text.length(); ++i) {
            char& nucleotide = complementaryKmer._text[i];
            char currentNucleotide = nucleotide;
            bool found = false;

            // Iterate through nucleotides and their complements to find the  
            // corresponding complement
            int j = 0;
            while ( !found && (j < (int)nucleotides.length()) ) {
               
                if ( currentNucleotide == nucleotides[j] ) {
                    nucleotide = complementaryNucleotides[j];
                    found = true;

                }

                j++;
            }
            // If the nucleotide is not found in the nucleotides list, keep it
            // unchanged
            if (!found) {
                nucleotide = currentNucleotide;
            }
           
        }

        return complementaryKmer;
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

