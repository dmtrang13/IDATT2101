#include <iostream>
#include <vector>
#include <random>
#include <chrono>
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

void quickSort(vector<int>& A, int l,int h) {
    if (l < h) {
        int j = partition(A,l,h);
        quickSort(A,l,j);
        quickSort(A,j+1,h);
    }
}

int main() {
    // Øker hastigheten på I/O
    ios::sync_with_stdio(false);                // Slår av synkronisering med C std
    cin.tie(nullptr);                           // Løsner cin fra cout slik lesninger ikke tvinger fram flush
    int n = 1'000'000;
    vector<int> unique = make_unique(n);
    
    cout << "Sum før:  " << checksum(unique) << "\n";
    // Sorter kopien av lista (Copy)
    long long t_unique1 = bench_ns(unique, [&] (vector<int> &A){quickSort(A, 0, (int)A.size() - 1); }, BenchMode::Copy);
    cout << "Vanlig quicksort, med unike tall: " << t_unique1 << "ns\n";
    quickSort(unique, 0, (int)unique.size() - 1);
    if (!is_sorted(unique.begin(), unique.end())) std::cerr << "Ikke sortert!\n";
    cout << "Sum etter: " << checksum(unique) << "\n";

    // Resorter faktisk lista (InPlace)
    long long t_unique2 = bench_ns(unique, [&] (vector<int> &A){quickSort(A, 0, (int)A.size() - 1); }, BenchMode::InPlace);
    // Skriver ut med 4 desimaler
    cout << "Sortere samme lista igjen: " << t_unique2 << "ns (forhold "
        << fixed << setprecision(4) << (double)t_unique2 / max(1LL, t_unique1) << "x)\n";

    vector<int> dups = make_dups(n);
    cout << "Sum før:  " << checksum(dups) << "\n";
    long long t_dups = bench_ns(dups, [&] (vector<int> &A){quickSort(A, 0, (int)A.size() - 1); });
    cout << "Vanlig quicksort, med mange duplikater: " << t_dups << "ns\n";
    quickSort(dups, 0, (int)dups.size() - 1);
    if (!is_sorted(dups.begin(), dups.end())) std::cerr << "Ikke sortert!\n";
    cout << "Sum etter: " << checksum(dups) << "\n";
    // Resorter faktisk lista (InPlace)
    long long t_dups2 = bench_ns(dups, [&] (vector<int> &A){quickSort(A, 0, (int)A.size() - 1); }, BenchMode::InPlace);
    // Skriver ut med 4 desimaler
    cout << "Sortere samme lista igjen: " << t_dups2 << "ns (forhold "
        << fixed << setprecision(4) << (double)t_dups2 / max(1LL, t_dups) << "x)\n";

    return 0;
}