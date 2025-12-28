#include <iostream>
#include <fstream>
#include <string>
#include <codecvt>
#include <locale>
#include <cstdint>
#include <vector>
#include <iomanip>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>

using namespace std;

// Generer nøkler
vector<uint64_t> make_keys(uint64_t m) {
    vector<uint64_t> A(m);
    uint64_t x = 0;              
    mt19937 rng(987654);
    uniform_int_distribution<uint64_t> step(1, 1000);
    for (size_t i = 0; i < m; i++) {
        x += step(rng);     // unik økende
        A[i] = x;
    }
    shuffle(A.begin(), A.end(), rng);   // Stokke lista
    return A;
}

// Polynomial rolling 
static size_t hash_poly(uint64_t key, size_t tableSize) {
    return (key * 11400714819323198485ULL) % tableSize;     // Knuths multiplikativ
}

class HashTable {
public:  
    virtual void insert(uint64_t key) = 0;
    virtual bool contains(uint64_t key) const = 0;
    virtual void printCollisions() const = 0;
    virtual ~HashTable() {}
protected: 
    size_t tableSize;
    explicit HashTable(size_t sz) : tableSize(sz) {}
};

class HashTableLinear : public HashTable {
public: 
    static constexpr size_t fixedTableSize = 8388608;  

private:
    vector<uint64_t> table;
    vector<uint8_t> used;
    size_t inserted = 0;
    size_t totalCollisions = 0;
    size_t maxProbes = 0;

public:
    explicit HashTableLinear(size_t sz) : HashTable(sz), table(sz), used(sz, false) {}
    
    void reset() {
        fill(used.begin(), used.end(), uint8_t{0});
        inserted = 0;
        totalCollisions = 0;
        maxProbes = 0;
    }

    void insert(uint64_t key) override {
        size_t idx = hash_poly(key, tableSize);
        size_t probes = 1;
        size_t start = idx;
        while (used[idx]) {
            idx = (idx + 1) & (tableSize - 1);     // linear probing
            probes++;
            if (idx == start) {
                throw runtime_error("Lineær probing: full tabell.");
            }
        }
        table[idx] = key;
        used[idx] = true;    
        inserted++;    

        totalCollisions += (probes - 1);
        if (probes > maxProbes) maxProbes = probes;
    }

    bool contains(uint64_t key) const override {
        size_t idx = hash_poly(key, tableSize);  // finn "home" slot
        size_t start = idx;
        while (used[idx]) {
            if (table[idx] == key) return true;
            // Bitmasking
            idx = (idx + 1) & (tableSize - 1);     //Flytt til neste slot om false
            if (idx == start) break;
        }
        return false;
    }

    void printCollisions() const override {
        cout << "\nTotall elementer: " << inserted << "\n";
        cout << "Totall kollisjoner: " << totalCollisions << "\n";
        // unngå å dele på 0
        cout << fixed << setprecision(3)
             << "Kollisjoner per insetting: " << (inserted ? (double)totalCollisions / inserted : 0.0) << "\n";
    }
};

class HashTableDouble : public HashTable {
public: 
    static constexpr size_t fixedTableSize = 8'000'009;

private:
    vector<uint64_t> table;
    vector<uint8_t> used;
    size_t inserted = 0;
    size_t totalCollisions = 0;
    size_t maxProbes = 0;

public:
    explicit HashTableDouble(size_t sz) : HashTable(sz), table(sz), used(sz, false) {}
    
    void reset() {
        fill(used.begin(), used.end(), uint8_t{0});
        inserted = 0;
        totalCollisions = 0;
        maxProbes = 0;
    }


    static size_t h2(uint64_t key, size_t tableSize) {
        return 1 + ((key * 7046029254386353131ULL) % (tableSize - 1));
    }

    void insert(uint64_t key) override {
        size_t idx = hash_poly(key, tableSize);
        size_t probes = 1;
        if (!used[idx]) {
            table[idx] = key;
            used[idx] = true;
            inserted++;
            if (probes > maxProbes) maxProbes = probes;
            return;
        }

        size_t step = h2(key, tableSize);   // Regne h2 bare når det oppstår kollisjoner
        size_t start = idx;
        while (used[idx]) {
            idx = (idx + step) % tableSize;
            probes++;
            if (idx == start) {
                throw runtime_error("Dobbel hashing: fant ikke slot.");
            }
        }

        table[idx] = key;
        used[idx] = true;
        inserted++;

        totalCollisions += (probes - 1);
        if (probes > maxProbes) maxProbes = probes;
    }

    bool contains(uint64_t key) const override {
        size_t idx = hash_poly(key, tableSize);
        size_t step = h2(key, tableSize);
        size_t start = idx;

        while (used[idx]) {
            if (table[idx] == key) return true;
            idx = (idx + step) % tableSize;
            if (idx == start) break;    // full sirkel
        }
        return false;
    }

    void printCollisions() const override {
        cout << "\nTotall elementer: " << inserted << "\n";
        cout << "Totall kollisjoner: " << totalCollisions << "\n";
        // unngå å dele på 0
        cout << fixed << setprecision(3)
             << "Kollisjoner per insetting: " << (inserted ? (double)totalCollisions / inserted : 0.0) << "\n";
    }
};


int main() {

    const size_t ML = HashTableLinear::fixedTableSize;
    const size_t MD = HashTableDouble::fixedTableSize;
    vector<uint64_t> keysLin = make_keys(ML);
    vector<uint64_t> keysDub = make_keys(MD);
    HashTableLinear htLin(ML);
    HashTableDouble htDub(MD);

    mt19937_64 rng(123456789);
    uniform_int_distribution<size_t> pickStartL(0, ML-1);
    uniform_int_distribution<size_t> pickStepL(1, ML-1);
    uniform_int_distribution<size_t> pickStartD(0, MD-1);
    uniform_int_distribution<size_t> pickStepD(1, MD-1);

    vector<double> fills = {0.5, 0.8, 0.9, 0.99, 1.0};
    for (double f : fills) {
        const size_t kL = (f>= 1.0) ? ML : size_t(f * ML);

        size_t startL = pickStartL(rng);
        size_t stepL = pickStepL(rng);

        htLin.reset();
        auto t0 = chrono::steady_clock::now();
        for (size_t done = 0, i = startL; done < kL; done++, i = (i + stepL) & (ML - 1)) {
            htLin.insert(keysLin[i]);
        }
        auto t1 = chrono::steady_clock::now();
        chrono::duration<double> dur1 = t1 - t0;

        cout << "\n[Lineær probing] Fyllingsgrad = " << (f*100) << "% (" << kL << "/" << ML << ")\n";
        htLin.printCollisions();
        cout << "Tid: " << dur1.count() << " sek\n";

        size_t kD = (f >= 1.0) ? MD : size_t(f * MD);
        size_t startD = pickStartD(rng);
        size_t stepD = pickStepD(rng);

        htDub.reset();
        auto t2 = chrono::steady_clock::now();
        for (size_t done = 0, i = startD; done < kD; done++, i = (i + stepD) % MD) {
            htDub.insert(keysDub[i]);
        }
        auto t3 = chrono::steady_clock::now();
        chrono::duration<double> dur2 = t3 - t2;

        cout << "\n[Dobbel hashing] Fyllingsgrad = " << (f*100) << "% (" << kD << "/" << MD << ")\n";
        htDub.printCollisions();
        cout << "Tid: " << dur2.count() << " sek\n";
    }

    return 0;
}