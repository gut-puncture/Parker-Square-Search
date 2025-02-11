#include <iostream>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

// Check if x is a perfect square.
inline bool isPerfectSquare(long long x) {
    if(x < 0) return false;
    long long r = (long long)std::llround(std::sqrt((long double)x));
    return r*r == x;
}

int main(){
    // Set maximum for a, b, c (their values are the square roots of our entries).
    const int maxVal = 100000; // adjust as needed
    int count = 0;
    // Iterate over possible a, b, c.
    // Note: a, b, c must be >=1 so that a^2, b^2, c^2 > 0.
    #pragma omp parallel for collapse(3) schedule(dynamic)
    for (int a = 1; a <= maxVal; a++){
        for (int b = 1; b <= maxVal; b++){
            for (int c = 1; c <= maxVal; c++){
                // Our base squares:
                long long A = (long long)a * a;    // top-left = a²
                long long E = (long long)b * b;    // center = b²
                long long C = (long long)c * c;    // top-right = c²

                // Derived entries:
                long long B = 3 * E - A - C;       // top-middle
                long long D = E + C - A;           // middle-left
                long long F = E + A - C;           // middle-right
                long long G = 2 * E - C;           // bottom-left
                long long I = 2 * E - A;           // bottom-right
                long long H = A + C - E;           // bottom-middle

                // Check none of the derived numbers is negative.
                if(B < 0 || D < 0 || F < 0 || G < 0 || H < 0 || I < 0)
                    continue;
                // Check each derived value is a perfect square.
                if(!isPerfectSquare(B)) continue;
                if(!isPerfectSquare(D)) continue;
                if(!isPerfectSquare(F)) continue;
                if(!isPerfectSquare(G)) continue;
                if(!isPerfectSquare(H)) continue;
                if(!isPerfectSquare(I)) continue;

                // Check all 9 numbers are distinct.
                long long sq[9] = { A, B, C, D, E, F, G, H, I };
                bool distinct = true;
                for(int i = 0; i < 9 && distinct; i++){
                    for(int j = i + 1; j < 9; j++){
                        if(sq[i] == sq[j]) { distinct = false; break; }
                    }
                }
                if(!distinct) continue;

                // Found a valid magic square of perfect squares.
                #pragma omp critical
                {
                    std::cout << "Magic square found (a=" << a << ", b=" << b << ", c=" << c << "):\n";
                    std::cout << A << "\t" << B << "\t" << C << "\n";
                    std::cout << D << "\t" << E << "\t" << F << "\n";
                    std::cout << G << "\t" << H << "\t" << I << "\n";
                    std::cout << "--------------------------\n";
                    count++;
                }
            }
        }
    }
    std::cout << "Total magic squares found: " << count << "\n";
    return 0;
}
