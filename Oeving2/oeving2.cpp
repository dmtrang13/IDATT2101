#include <iostream>
#include <chrono>
using namespace std;

int recursion1 (int n, int x) {
    if (n == 1) {
        return n, x;
    }
    return x + recursion1(n-1,x);
}

int recursion2(int n, int x) {
    if (n == 1) {
        return n, x;
    }
    if (n%2 == 0) {
        return recursion2(n/2, x+x);
    }
    return x + recursion2((n-1)/2, x+x);
}

int main() {
    int n = 43339;
    int x = 2;
    chrono::steady_clock::time_point begin1 = chrono::steady_clock::now();
    recursion1(n,x);
    chrono::steady_clock::time_point end1 = chrono::steady_clock::now();

    cout << "Tidsmåling: " << chrono::duration_cast<chrono::nanoseconds> (end1 - begin1).count() << "[ns]\n";

    chrono::steady_clock::time_point begin2 = chrono::steady_clock::now();
    recursion2(n,x);
    chrono::steady_clock::time_point end2 = chrono::steady_clock::now();

    cout << "Tidsmåling: " << chrono::duration_cast<chrono::nanoseconds> (end2 - begin2).count() << "[ns]\n";
}
