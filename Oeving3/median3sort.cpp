#include <iostream>
#include <vector>
#include <chrono>

using namespace std;

// Partially sorts a[v], a[m], a[h] so a[v] <= a[m] <= a[h]; returns m
int median3sort(vector<int>& a, int v, int h) {
    int m = (v + h) / 2;
    if (a[v] > a[m]) swap(a[v], a[m]);
    if (a[m] > a[h]) {
        swap(a[m], a[h]);
        if (a[v] > a[m]) swap(a[v], a[m]);
    }
    return m;
}

// Sedgewick/Hoare-style partition with pivot stored at h-1
int partition_median3(vector<int>& a, int v, int h) {
    int m = median3sort(a, v, h);
    int pivot = a[m];
    swap(a[m], a[h - 1]);         // move pivot out of the way
    int iv = v;                   // left scanner
    int ih = h - 1;               // right scanner (pivot is at h-1)

    for (;;) {
        while (a[++iv] < pivot) ; // safe due to sentinels from median3sort
        while (a[--ih] > pivot) ;
        if (iv >= ih) break;
        swap(a[iv], a[ih]);
    }

    swap(a[iv], a[h - 1]);        // put pivot into its final place
    return iv;                    // return pivot index
}

void quicksort(vector<int>& a, int v, int h) {
    if (h - v > 2) {
        int p = partition_median3(a, v, h);
        quicksort(a, v, p - 1);
        quicksort(a, p + 1, h);
    } else {
        // Sort small segment (size 2 or 3) in place
        median3sort(a, v, h);
    }
}

int main() {
    vector<int> a = {10, 7, 8, 9, 1, 5};
    quicksort(a, 0, (int)a.size() - 1);
    for (int x : a) cout << x << ' ';
    return 0;
}
