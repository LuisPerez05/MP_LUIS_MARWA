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

string Profile::getProfileId() { 
    return _profileId; 
} 

void Profile::setProfileId(string id) { 
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

void Profile::readFromFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("No se pudo abrir el archivo");
    }
    file >> _profileId;
    int size;
    file >> size;
    if (size < 0 || size > DIM_VECTOR_KMER_FREQ) {
        throw out_of_range("Tamaño no válido por el vector de Kmers en el profile");
    }
    _size = size;
    for (int i = 0; i < _size; ++i) {
        char nucleotide;
        int frequency;
        file >> nucleotide >> frequency;
        _vectorKmerFreq[i] = KmerFreq(Kmer(nucleotide), frequency);
    }
    file.close();
}

KmerFreq& Profile::at(int index) {
    if (index < 0 || index >= _size) {
        throw out_of_range("Índice fuera de rango en el perfil");
    }
    return _vectorKmerFreq[index];
}

const KmerFreq& Profile::at(int index) const {
    if (index < 0 || index >= _size) {
        throw out_of_range("Índice fuera de rango en el perfil");
    }
    return _vectorKmerFreq[index];
}

int Profile::getSize() {
    return _size;
}

int Profile::getCapacity() {
    return DIM_VECTOR_KMER_FREQ;
}
int Profile::findKmer(Kmer kmer, int initialPos, int finalPos) {
    if (initialPos < 0 || finalPos >= _size || initialPos > finalPos) {
        throw out_of_range("Posiciones iniciales/finales inválidas para la búsqueda de kmer");
    }
    for (int i = initialPos; i <= finalPos; ++i) {
        if (_vectorKmerFreq[i].getKmer() == kmer) {
            return i;
        }
    }
    return -1;
}

int Profile::findKmer(Kmer kmer) {
    return findKmer(kmer, 0, _size - 1);
}
void Profile::sort() {
    std::sort(_vectorKmerFreq, _vectorKmerFreq + _size, [](const KmerFreq& a, const KmerFreq& b) {
        if (a.getFreq() != b.getFreq()) {
            return a.getFreq() > b.getFreq();
        } else {
            return a.getKmer().toString() < b.getKmer().toString();
        }
    });
}

void Profile::save(char fileName[]) {
    ofstream file(fileName);
    if (!file.is_open()) {
        throw ios_base::failure("No se pudo abrir el archivo para escritura");
    }
    file << MAGIC_STRING_T << endl;
    file << _profileId << endl;
    file << _size << endl;
    for (int i = 0; i < _size; ++i) {
        file << _vectorKmerFreq[i].getKmer().toString() << " " << _vectorKmerFreq[i].getFreq() << endl;
    }
    file.close();
}

void Profile::load(char fileName[]) {
    ifstream file(fileName);
    if (!file.is_open()) {
        throw ios_base::failure("No se pudo abrir el archivo para lectura");
    }
    string magicString;
    getline(file, magicString);
    if (magicString != MAGIC_STRING_T) {
        throw invalid_argument("Cadena mágica no válida");
    }
    getline(file, _profileId);
    int size;
    file >> size;
    if (size < 0 || size > DIM_VECTOR_KMER_FREQ) {
        throw out_of_range("Número de kmers no válido");
    }
    _size = size;
    for (int i = 0; i < _size; ++i) {
        char nucleotide;
        int frequency;
        file >> nucleotide >> frequency;
        _vectorKmerFreq[i] = KmerFreq(Kmer(nucleotide), frequency);
    }
    file.close();
}

void Profile::append(KmerFreq kmerFreq) {
    if (_size >= DIM_VECTOR_KMER_FREQ) {
        throw out_of_range("Perfil completo, no se pueden añadir más KmerFreq");
    }
    for (int i = 0; i < _size; ++i) {
        if (_vectorKmerFreq[i].getKmer() == kmerFreq.getKmer()) {
            _vectorKmerFreq[i].setFreq(_vectorKmerFreq[i].getFreq() + kmerFreq.getFreq());
            return;
        }
    }
    _vectorKmerFreq[_size++] = kmerFreq;
}
void Profile::normalize(string validNucleotides) {
    // Convertir todos los caracteres a mayúsculas y reemplazar los inválidos con MISSING_NUCLEOTIDE
    for (int i = 0; i < _size; ++i) {
        string kmerString = _vectorKmerFreq[i].getKmer().toString();
        for (char& c : kmerString) {
            if (validNucleotides.find(toupper(c)) == string::npos) {
                c = Kmer::MISSING_NUCLEOTIDE;
            } else {
                c = toupper(c);
            }
        }
        _vectorKmerFreq[i].setKmer(Kmer(kmerString));
    }

    // Fusionar kmers idénticos y actualizar las frecuencias
    for (int i = 0; i < _size; ++i) {
        for (int j = i + 1; j < _size; ++j) {
            if (_vectorKmerFreq[i].getKmer().toString() == _vectorKmerFreq[j].getKmer().toString()) {
                _vectorKmerFreq[i].setFreq(_vectorKmerFreq[i].getFreq() + _vectorKmerFreq[j].getFreq());
                // Eliminar el kmer repetido moviendo los kmers restantes hacia atrás
                for (int k = j; k < _size - 1; ++k) {
                    _vectorKmerFreq[k] = _vectorKmerFreq[k + 1];
                }
                _size--; // Decrementar el tamaño después de eliminar un kmer repetido
                j--; // Ajustar el índice para volver a comprobar si hay kmers repetidos
            }
        }
    }
}


