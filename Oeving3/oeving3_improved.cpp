#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <iomanip>

using namespace std;


// Generer data
vector<int> make_dups(int n) {
    vector<int> A(n);
    mt19937 rng(123456);                         // Mersenne Twister for random tall
    uniform_int_distribution<int> dist(0, 1'000'000);
    // Hver partall indeks = random tall, hver oddetall indeks = 42
    for (int i = 0; i < n; i++) A[i] = (i%2 == 0) ? dist(rng) : 42;
    return A;
}

vector<int> make_unique(int n) {
    vector<int> A(n);
    iota(A.begin(), A.end(), 1);                // Fyll inn lista med unike tall
    mt19937 rng(987654);
    shuffle(A.begin(), A.end(), rng);           // Stokke lista
    return A;
}

long long checksum(vector<int>& A) {
    return std::accumulate(A.begin(), A.end(), 0LL);
}

// Benchmark for å måle programmets ytelse (for både sortert lista og resortert lista)
enum class BenchMode { Copy, InPlace };
template <class Sort>
long long bench_ns(vector<int> &a, Sort sorter, BenchMode mode = BenchMode::Copy) {

    vector<int> tmp;                             // Lagre for copy-modus
    vector<int> *A = &a;                         // Sorter original liste (InPlace)
    if (mode == BenchMode::Copy) {
        tmp = a;                                 // Lage kopien
        A = &tmp;                                // Sorterer kopien
    }

    long long before = checksum(*A);             

    chrono::steady_clock::time_point t0 = chrono::steady_clock::now();
    sorter(*A);                                  // Sorter en kopi av lista
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
    if(!is_sorted(A->begin(), A->end())) cerr << "Lista er ikke sortert.\n";
    long long after = checksum(*A);
    if (before != after)
        std::cerr << "Sjekksum mismatch!\n";

    long long ns = chrono::duration_cast<chrono::nanoseconds>(t1-t0).count();
    return ns;
}

// Hoare oppdeling
int partition(vector<int>& A, int l, int h) {
    int pivot = A[(l + h) / 2];
    int i = l-1; int j = h+1;
    
    while (true) {
        do { i++; } while (A[i] < pivot);
        do { j--; } while (A[j] > pivot);
        if (i >= j) return j;
        swap(A[i],A[j]);
    }
}

// Sett opp sentinels A[0] og A[n-1]
void set_sentinels(vector<int>& A) {
    int n = (int)A.size();
    if (n <= 1) return;         
    int min_i = 0, max_i = 0;
    for (int i = 1; i < n; ++i) {
        if (A[i] < A[min_i]) min_i = i;
        if (A[i] > A[max_i]) max_i = i;         // Holde styr over minst i og størst i
    }
    if (A[min_i] == A[max_i]) return;      
    swap(A[0], A[min_i]);                       // Alle elementene er like => putt minste i i begynnelsen
    if (max_i == 0) max_i = min_i;              // Størst i i begynnelsen => bytt posisjon med minst i
    swap(A[n-1], A[max_i]);                     // Flytt størst i mot slutten av lista
}

// Innsettingsort for små deltabeller
void insertion_sort(vector<int> &A, int l, int h) {
    for (int i = l+1; i <= h; i++) {
        int x = A[i], j = i - 1;                // Sett x som nøkkelverdi, start sammenligning med resten 
        while (j >= l && A[j] > x) {            // Sjekk grensen og om elementet er større enn x
            A[j+1] = A[j];                      // Flytt det større elementet til høyre, åpning ved j
            j--;                                // Fortsett å flytte j til venstre til x er på riktig plass
        }
        A[j+1] = x;                             // Putt x i dens riktig plass
    }
}

// Sorterer deltabellen (indre område [l...h], l >= 1, h <= n-2) 
void quicksort_core(vector<int>& A, int l, int h) {
    int cutOff = 24;                            // Innsettingsort cutoff

    while (l < h) {
        // Sikkerhetssjek slik at den leser ikke utenfor grensen 
        assert(l >= 1 && h + 1 < (int)A.size());  
        if (A[l-1] == A[h+1]) return;           // Om alle elementene er like
        // Bruk quicksort for stor liste (indekser > cutOff), insettingsort for deltabeller (indekser < cutOff)
        if (h - l + 1 <= cutOff) {              
            insertion_sort(A, l, h);
            return; 
        }

        int j = partition(A, l, h);  
        if (j - l < h - j) {                    // Om venstre side er mindre
            quicksort_core(A, l, j);            // Rekursere den mindre siden [l...j]
            l = j + 1;                          // Halerekursere på den større siden via å løkke gjennom [j+1...h]
        } else {                                // Høyre side mindre eller like stor
            quicksort_core(A, j + 1, h);        // Rekursere den mindre siden [j+1...h]
            h = j;                              // Løkka fortsetter på større siden
        }
    }
}

void quicksort_improved(vector<int>& A) {
    int n = (int)A.size();
    if (n <= 2) {
        if (n == 2 && A[1] < A[0]) swap(A[0], A[1]); 
        return;
    }

    set_sentinels(A);
    if (A.front() == A.back()) return;          // Alle elementene er like
    quicksort_core(A, 1, n-2);                  // Sorterer bare det indre området
}


int main() {
    // Øker hastigheten på I/O
    ios::sync_with_stdio(false);                // Slår av synkronisering med C std
    cin.tie(nullptr);                           // Løsner cin fra cout slik lesninger ikke tvinger fram flush
    int n = 1'000'000;
    vector<int> unique = make_unique(n);
    
    cout << "Sum før:  " << checksum(unique) << "\n";
    // Sorter kopien av lista (Copy)
    long long t_unique1 = bench_ns(unique, [&] (vector<int> &A){quicksort_improved(A); }, BenchMode::Copy);
    cout << "Improvisert quicksort, med unike tall: " << t_unique1 << "ns\n";
    quicksort_improved(unique);
    cout << "Sum etter: " << checksum(unique) << "\n";

    // Resorter faktisk lista (InPlace)
    long long t_unique2 = bench_ns(unique, [&] (vector<int> &A){quicksort_improved(A); }, BenchMode::InPlace);
    // Skriver ut med 4 desimaler
    cout << "Sortere samme lista igjen: " << t_unique2 << "ns (forhold "
        << fixed << setprecision(4) << (double)t_unique2 / max(1LL, t_unique1) << "x)\n";

    vector<int> dups = make_dups(n);
    cout << "Sum før:  " << checksum(dups) << "\n";
    long long t_dups = bench_ns(dups, [&] (vector<int> &A){quicksort_improved(A); });
    cout << "Improvisert quicksort, med mange duplikater: " << t_dups << "ns\n";
    quicksort_improved(dups);
    cout << "Sum etter: " << checksum(dups) << "\n";
    // Resorter faktisk lista (InPlace)
    long long t_dups2 = bench_ns(dups, [&] (vector<int> &A){quicksort_improved(A); }, BenchMode::InPlace);
    // Skriver ut med 4 desimaler
    cout << "Sortere samme lista igjen: " << t_dups2 << "ns (forhold "
        << fixed << setprecision(4) << (double)t_dups2 / max(1LL, t_dups) << "x)\n";

    return 0;
}