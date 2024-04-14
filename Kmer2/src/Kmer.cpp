/*
 * Metodología de la Programación: Kmer2
 * Curso 2023/2024
 */

/** 
 * @file Kmer.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 31 October 2023, 14:15
 */

#include "Kmer.h"

using namespace std;

// Constructor para inicializar un objeto Kmer con una secuencia de longitud k
Kmer::Kmer(int k){
    if (k<= 0){ // Verifica si k es menor o igual a cero
        throw std::invalid_argument("k debe ser mayor a 0"); // Lanza un mensaje si k es menor o igual a cero
    }    
    std::string text(k, MISSING_NUCLEOTIDE);  // Crea una cadena de texto de longitud k con caracteres MISSING_NUCLEOTIDE
    Kmer::_text = text;   // Asigna la cadena de texto al atributo _text de la clase Kmer
       
}
   
// Constructor para inicializar un objeto Kmer con una secuencia proporcionada   
Kmer::Kmer(const std::string& text){
    if (text.empty()) { // Verifica si la cadena de texto está vacía
        throw std::invalid_argument("El texto no puede estar vacio"); // Lanza un mensaje si la cadena de texto está vacía
    }
    Kmer::_text = text; // Asigna la cadena de texto al atributo _text de la clase Kmer
}

// Devuelve la longitud de la secuencia de ADN   
int Kmer::getK() const{
    return ((Kmer::_text).size());
}
   
// Devuelve la longitud de la secuencia de ADN   
int Kmer::size() const {
    return ((Kmer::_text).size());
}
   
// Devuelve la secuencia de ADN como una cadena de texto   
std::string Kmer::toString() const{
    return (Kmer::_text);
}

// Devuelve una referencia al carácter en la posición index de la secuencia de ADN (const)
const char& Kmer::at(int index) const{
    if(index>(Kmer::_text.size()-1)){  // Verifica si el índice está dentro del rango de la secuencia
        throw std::out_of_range("Indice fuera del rango del k-mero"); // Lanza un mensaje si el índice está fuera del rango
    }
    return (Kmer::_text.at(index)); // Devuelve una referencia al carácter en la posición index
}


// Devuelve una referencia al carácter en la posición index de la secuencia de ADN
char& Kmer::at(int index){
    if(index>(Kmer::_text.size()-1)){ // Verifica si el índice está dentro del rango de la secuencia
        throw std::out_of_range("indice fuera del rango del k-mero"); // Lanza un mensaje si el índice está fuera del rango
    }
    return (Kmer::_text.at(index)); // Devuelve una referencia al carácter en la posición index
}
 

// Normaliza la secuencia de ADN utilizando la lista de nucleótidos válidos
void Kmer::normalize(const std::string& validNucleotides){
    int k = Kmer::_text.size(); // Obtiene la longitud de la secuencia de ADN
    std::string solution; // Cadena de texto para almacenar la secuencia normalizada
    for(int n = 0; n < k; n++){  // Interaciona sobre cada carácter de la secuencia
        char actual = Kmer::_text.at(n);  // Obtiene el carácter actual
        if(IsValidNucleotide(actual, validNucleotides)){ // Verifica si el carácter actual es un nucleótido válido
        }else{  // Si es válido, no hace nada
            bool minus = false; // Si no es válido, lo normaliza
            int pos = 0;
            int max = validNucleotides.size();
            while(!minus && (pos < max)){   // Busca el carácter en la lista de nucleótidos válidos
                if(actual == tolower(validNucleotides.at(pos))) minus = true;
                pos++;
            }
            if(minus){ // Si se encuentra en minúsculas, lo convierte a mayúsculas; de lo contrario, lo marca como MISSING_NUCLEOTIDE
                actual = toupper(actual);
            }else{
                actual = MISSING_NUCLEOTIDE;
            }
        }
        solution = solution + actual;  // Agrega el carácter normalizado a la solución
    }
    Kmer::_text = solution; // Actualiza la secuencia de ADN con la solución normalizada
}
   

// Devuelve un nuevo objeto Kmer que representa el complemento de la secuencia de ADN actual
Kmer Kmer::complementary(const std::string& nucleotides,const std::string& complementaryNucleotides) const {
    if (nucleotides.size() != complementaryNucleotides.size()) { // Verifica si los tamaños de las listas de nucleótidos y sus complementarios coinciden
        throw std::invalid_argument("El tamaño de los nucleotidos y de los nucleotidos complementarios no coinciden"); // Lanza un mensaje si los tamaños no coinciden
    }
    int k = Kmer::_text.size(); // Obtiene la longitud de la secuencia de ADN
    std::string solution; // Cadena de texto para almacenar la solución
    
    for(int i = 0; i < k; i++){ // Interacciona sobre cada carácter de la secuencia
        char actual = Kmer::_text.at(i); // Obtiene el carácter actual
        int max = nucleotides.size();
        bool same = false;
        int pos = 0;
        while ((!same) && (pos < max)){  // Busca el carácter en la lista de nucleótidos
            if(actual == nucleotides.at(pos)){
                same = true;
            }
            pos++;
        }
        if(same){ // Si se encuentra en la lista de nucleótidos, lo reemplaza por su complemento
            pos--;
            actual = complementaryNucleotides.at(pos);
        }
        solution = solution + actual; // Agrega el carácter a la solución
    }
    return(Kmer(solution)); // Devuelve un nuevo objeto Kmer con la solución
}
 
// Verifica si un carácter dado es un nucleótido válido
bool IsValidNucleotide(char nucleotide, const std::string& validNucleotides){
    bool es_valido = false; // Bandera para indicar si el carácter es válido
    int k = validNucleotides.size(); // Obtiene la longitud de la lista de nucleótidos válidos
    int i = 0;
    while((i < k) && (!es_valido)) {  // Interacciona sobre cada nucleótido válido
       if(nucleotide == validNucleotides.at(i)){ // Verifica si el carácter actual es un nucleótido válido
            es_valido = true; // Establece la bandera en verdadero si es válido
        }
        i++;
    }
    return(es_valido);  // Devuelve el resultado
}

// Convierte todos los caracteres de la secuencia de ADN a minúsculas
void ToLower(Kmer& kmer){
    int k = kmer.size(); // Obtiene la longitud de la secuencia de ADN
    for(int i = 0; i < k; i++){  // Interacciona sobre cada carácter de la secuencia
        kmer.at(i) = tolower(kmer.at(i)); // Convierte el carácter a minúsculas
    }
}


// Convierte todos los caracteres de la secuencia de ADN a mayúsculas
void ToUpper(Kmer& kmer){
    int k = kmer.size(); // Obtiene la longitud de la secuencia de ADN
    for(int i = 0; i < k; i++){ // Interacciona sobre cada carácter de la secuencia
        kmer.at(i) = toupper(kmer.at(i));  // Convierte el carácter a mayúsculas
    }   
}