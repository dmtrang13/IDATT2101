#include <iostream>
#include <fstream>
#include <string>
#include <codecvt>
#include <locale>
#include <cstdint>
#include <vector>
#include <iomanip>
#include <cstdlib>

using namespace std;

// Håndtere æ, ø, å 
static vector<uint32_t> utf8_codepoints(const string& s) {
    static wstring_convert<codecvt_utf8<char32_t>, char32_t> conv;
    u32string u32 = conv.from_bytes(s);
    vector<uint32_t> nameNum;
    nameNum.reserve(u32.size());
    for (char32_t cp : u32) {
        if (cp == U',' || cp == U' ') continue;
        nameNum.push_back(static_cast<uint32_t>(cp));
    }
    return nameNum;
}

// Polynomial rolling
static size_t hash_poly(const string& name, size_t tableSize, uint64_t B = 31) {
    uint64_t h = 0;
    for (uint32_t cp : utf8_codepoints(name)) {
        h = h * B + cp; // Horner
    }
    return static_cast<size_t>(h % tableSize);
}

// Valgbar komma eller mellomrom mellom navn og etternavn
static string canonical(const string& s) {
    string t; t.reserve(s.size());
    for (unsigned char c : s) if (c != ',' && c != ' ') t.push_back(c);
    return t;
}

// Lenket liste-node
struct Node {
    string key;
    Node* next;
    explicit Node(string k, Node* n=nullptr) : key(move(k)), next(n) {}
};

// Hashtabell med separat chaining
struct HashTable {
    vector<Node*> buckets;
    size_t tableSize;
    size_t count = 0;
    size_t insert_collisions = 0;
    size_t search_collisions = 0;

    explicit HashTable(size_t sz) : buckets(sz, nullptr), tableSize(sz) {}

    size_t indexOf(const string& name) const { return hash_poly(name, tableSize); }

    void insert(const string& name) {
        size_t idx = indexOf(name);
        if (buckets[idx] != nullptr) {
            cout << "Kollisjon (innsetting): [" << name << "] vs [" << buckets[idx]->key << "]\n";
            insert_collisions++; // <- viktig
        }
        buckets[idx] = new Node(name, buckets[idx]);
        count++;
    }

    bool contains(const string& name) {
        size_t idx = indexOf(name);
        Node* p = buckets[idx];
        if (!p) return false;

        const string qcan = canonical(name);
        if (canonical(p->key) == qcan) return true;

        for (Node* q = p->next; q; q = q->next) {
            cout << "Kollisjon (søk): [" << name << "] vs [" << p->key << "]\n";
            search_collisions++; // <- viktig
            if (canonical(q->key) == qcan) return true;
            p = q;
        }
        return false;
    }

    static size_t chainSize(Node* head) {
        size_t s = 0;
        for (Node* p = head; p; p = p->next) s++;
        return s;
    }
    void printCollisions() const {
        size_t totalCollision = 0;
        for (size_t i = 0; i < tableSize; i++) {
            size_t s = chainSize(buckets[i]);
            if (s > 1) {
                cout << "Indeks " << i << " (listestørrelse " << s << "): ";
                for (Node* p = buckets[i]; p; p = p->next) cout << "[" << p->key << "] ";
                cout << "\n";
                totalCollision += (s - 1);
            }
        }
        cout << "\nTotalt elementer: " << count << "\n";
    }

    void printLoadFactor() const {
        double load = tableSize ? (double)count / tableSize : 0.0;
        cout << fixed << setprecision(3)
             << "Lastefaktor: " << load << " (n = " << count << ", m = " << tableSize << ")\n";
    }

    //Totalt antall kollisjoner under innsetting og ratio
    void printStats() const {
        cout << "Totalt kollisjoner (innsetting): " << insert_collisions << "\n";
        cout << "Totalt kollisjoner (søk): " << search_collisions << "\n";
        cout << fixed << setprecision(3)
             << "Kollisjoner per person (innsetting): "
             << (count ? (double)insert_collisions / count : 0.0) << "\n";
    }

    ~HashTable() {
        for (Node* head : buckets) {
            while (head) { Node* next = head->next; delete head; head = next; }
        }
    }
};

int main() {
    ifstream nameFile("navn.txt");
    if (!nameFile) {
        cout << "Feil ved åpning av navn.txt" << endl;
        return EXIT_FAILURE;
    }

    const size_t tableSize = 128;   // Nærmest toerpotens
    HashTable ht(tableSize);

    string line;
    size_t lineCount = 0;
    while (getline(nameFile, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back(); // CRLF
        if (line.empty()) continue;
        ht.insert(line);
        lineCount++;
    }

    cout << "Sette inn " << lineCount << " navn i tabell med størrelse " << tableSize << ".\n";
    ht.printCollisions();  

    string q = "Trang Minh Duong"; 
    cout << "\nSøke opp \"" << q << "\": " << (ht.contains(q) ? "Funnet" : "Ikke funnet") << "\n\n";

    ht.printLoadFactor();
    ht.printStats(); 
}