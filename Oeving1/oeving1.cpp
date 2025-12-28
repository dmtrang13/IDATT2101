#include <iostream>
#include <vector>
#include <numeric>
#include <random>
#include <algorithm>
using namespace std;

vector<int> brute_force(vector<int> deltas) {
    vector<int> best_transaction = {0, 0, 0};
    int operations = 0;
    for (int i = 0; i < deltas.size(); i++) {
        for (int j = i + 1; j < deltas.size(); j++) {
            int profit = reduce(deltas.begin() + i + 1, deltas.begin() + j);
            operations++;
            if (profit > best_transaction[2]) {
                best_transaction[0] = i + 1;
                best_transaction[1] = j;
                best_transaction[2] = profit;
            }
        }
    }
    cout << operations << "\n";
    return best_transaction;
}

int main() {
    //vector<int> deltas = {-1, 3, -9, 2, 2, -1, 2, -1, -5}; //Opprinnelig datasett
    for (int i = 0; i < 10; i++) {
        vector<int> deltas(i);
        generate(deltas.begin(), deltas.end(), [&](){ return rand() % 21 - 10; }); //Random tall for asymptotisk analyse
        vector<int> best = brute_force(deltas);


        cout << best[0] << " " << best[1] << " " << best[2] << "\n";
    }

    return 0;
}