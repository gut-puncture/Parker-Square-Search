#include <stdio.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h> // Include OpenMP header

// ... (PellSolution struct, is_perfect_square, generate_pell_solutions - same as before) ...
// Structure to hold a solution to the Pell equation using GMP integers.
typedef struct {
    mpz_t N;
    mpz_t M;
} PellSolution;

// Function to check if a GMP integer is a perfect square.
bool is_perfect_square(mpz_t n) {
    if (mpz_sgn(n) < 0) {
        return false;
    }
    mpz_t root;
    mpz_init(root);
    mpz_sqrt(root, n);
    mpz_mul(root, root, root);
    bool result = mpz_cmp(root, n) == 0;
    mpz_clear(root);
    return result;
}

// Function to generate solutions to Pell equation using GMP
int generate_pell_solutions(mpz_t D, unsigned long long limit, PellSolution *sols) {
    int count = 0;
    mpz_t M, N2, temp;
    mpz_init(M);
    mpz_init(N2);
    mpz_init(temp);

    for (unsigned long long m = 1; m <= limit; ++m) {
        mpz_set_ui(M, m);
        mpz_mul(temp, M, M);
        mpz_mul_ui(temp, temp, 2);
        mpz_add(N2, D, temp);

        if (mpz_sgn(N2) < 0) {
            continue;
        }
        if (is_perfect_square(N2)) {
            mpz_t N;
            mpz_init(N);
            mpz_sqrt(N, N2);
            mpz_init(sols[count].N);
            mpz_init(sols[count].M);
            mpz_set(sols[count].N, N);
            mpz_set(sols[count].M, M);
            count++;

            if (mpz_sgn(N) != 0) {
                mpz_init(sols[count].N);
                mpz_init(sols[count].M);
                mpz_neg(sols[count].N, N);
                mpz_set(sols[count].M, M);
                count++;
            }
            if (mpz_sgn(M) != 0) {
                mpz_init(sols[count].N);
                mpz_init(sols[count].M);
                mpz_set(sols[count].N, N);
                mpz_neg(sols[count].M, M);
                count++;
            }
            if (mpz_sgn(N) != 0 && mpz_sgn(M) != 0) {
                mpz_init(sols[count].N);
                mpz_init(sols[count].M);
                mpz_neg(sols[count].N, N);
                mpz_neg(sols[count].M, M);
                count++;
            }
            mpz_clear(N);
        }
    }

    mpz_clear(M);
    mpz_clear(N2);
    mpz_clear(temp);
    return count;
}
// --- Checkpointing Functions ---

// Loads the checkpoint from the file.
void load_checkpoint(unsigned long long *D_val, int *total_candidates) {
    FILE *fp = fopen("/content/drive/MyDrive/magic_square_checkpoint.txt", "r");
    if (fp == NULL) {
        // File doesn't exist (first run), or error opening.
        *D_val = 0; // Will be set to D_min later
        *total_candidates = 0;
        printf("No checkpoint file found. Starting from beginning.\n");
        return;
    }

    if (fscanf(fp, "%llu %d", D_val, total_candidates) != 2) {
        // File is corrupted or in the wrong format.
        fprintf(stderr, "Error reading checkpoint file. Starting from beginning.\n");
        *D_val = 0; // Will be set to D_min later
        *total_candidates = 0;
    } else {
        printf("Resuming from D_val = %llu, total_candidates = %d\n", *D_val, *total_candidates);
    }

    fclose(fp);
}

// Saves the checkpoint to the file.
void save_checkpoint(unsigned long long D_val, int total_candidates) {
    FILE *fp = fopen("/content/drive/MyDrive/magic_square_checkpoint.txt", "w");
    if (fp == NULL) {
        fprintf(stderr, "Error opening checkpoint file for writing!\n");
        return; // Don't exit, just try again next time.
    }

    fprintf(fp, "%llu %d\n", D_val, total_candidates);
    printf("Checkpoint saved: D_val = %llu, total_candidates = %d\n", D_val, total_candidates); // Added logging
    fclose(fp);
}



int search_magic_square(unsigned long long pell_limit, unsigned long long D_min, unsigned long long D_max) {
    int total_candidates = 0;
    unsigned long long current_D_val = D_min;
    unsigned long long log_interval = 1000; // Log progress every 1000 D values

    // Load checkpoint (if any).
    load_checkpoint(&current_D_val, &total_candidates);
    if (current_D_val < D_min || current_D_val>D_max) {
        current_D_val = D_min; // Ensure we start within the correct range.
    }


    // Allocate enough space.
    PellSolution *sols = malloc(4ULL * pell_limit * 2 * sizeof(PellSolution));

    if (sols == NULL) {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    mpz_t D_val, Y, X, V, U, e1, e2, e, a_val, c_val;
    mpz_t b2, d2, f2, h2, i2, g2;
    mpz_t A, B, C, D_cell, E, F, G, H, I;
    // Initialize mpz_t variables...
    mpz_init(D_val);
    mpz_init(Y);
    mpz_init(X);
    mpz_init(V);
    mpz_init(U);
    mpz_init(e1);
    mpz_init(e2);
    mpz_init(e);
    mpz_init(a_val);
    mpz_init(c_val);
    mpz_init(b2);
    mpz_init(d2);
    mpz_init(f2);
    mpz_init(h2);
    mpz_init(i2);
    mpz_init(g2);
    mpz_init(A);
    mpz_init(B);
    mpz_init(C);
    mpz_init(D_cell);
    mpz_init(E);
    mpz_init(F);
    mpz_init(G);
    mpz_init(H);
    mpz_init(I);


    // --- OpenMP Parallel Region ---
    #pragma omp parallel for reduction(+:total_candidates) shared(sols) firstprivate(pell_limit, D_max, D_val, Y, X, V, U, e1, e2, e, a_val, c_val, b2, d2, f2, h2, i2, g2, A, B, C, D_cell, E, F, G, H, I)
    for (unsigned long long d = current_D_val; d <= D_max; ++d) {
          //Re-init local vars
        mpz_set_ui(D_val, d);
        int num_sols = generate_pell_solutions(D_val, pell_limit, sols);
        if (num_sols == 0) {
              continue;
        }
        // Iterate through all pairs of solutions.
        for (int i = 0; i < num_sols; ++i) {
            for (int j = 0; j < num_sols; ++j) {

                mpz_set(Y, sols[i].N);
                mpz_set(X, sols[i].M);
                mpz_set(V, sols[j].N);
                mpz_set(U, sols[j].M);

                if (mpz_cmp(Y,V) == 0 && mpz_cmp(X,U)==0) {
                    continue;
                }

                // Calculate e, a, and c.
                mpz_mul(e1, Y, Y);
                mpz_mul(e, X, X);
                mpz_mul_ui(e,e,2);
                mpz_add(e1,e1,e);

                mpz_mul(e2,V,V);
                mpz_mul(e,U,U);
                mpz_mul_ui(e,e,2);
                mpz_add(e2,e2,e);

                if(mpz_cmp(e1,e2) != 0) continue;
                mpz_div_ui(e,e1,2);

                mpz_mul(a_val,Y,Y);
                mpz_mul(e,X,X);
                mpz_mul_ui(e,e,2);
                mpz_sub(a_val,a_val,e);
                mpz_div_ui(a_val, a_val, 2);

                mpz_mul(c_val,V,V);
                mpz_mul(e,U,U);
                mpz_mul_ui(e,e,2);
                mpz_sub(c_val, c_val, e);
                mpz_div_ui(c_val,c_val,2);

                // Calculate the other entries.
                mpz_mul(b2, e, e);
                mpz_mul_ui(b2, b2, 3);
                mpz_mul(e, a_val, a_val);
                mpz_sub(b2, b2, e);
                mpz_mul(e, c_val, c_val);
                mpz_sub(b2, b2, e);

                mpz_mul(d2,e,e);
                mpz_mul(f2, a_val, a_val);
                mpz_sub(d2, d2, f2);
                mpz_mul(f2, c_val, c_val);
                mpz_add(d2, d2, f2);

                mpz_mul(f2, e,e);
                mpz_mul(h2, a_val,a_val);
                mpz_add(f2,f2,h2);
                mpz_mul(h2,c_val,c_val);
                mpz_sub(f2,f2,h2);

                mpz_mul(h2, a_val, a_val);
                mpz_mul(i2, c_val, c_val);
                mpz_add(h2,h2,i2);
                mpz_mul(i2,e,e);
                mpz_sub(h2,h2,i2);

                mpz_mul(i2, e, e);
                mpz_mul_ui(i2, i2, 2);
                mpz_mul(g2, a_val, a_val);
                mpz_sub(i2,i2,g2);

                mpz_mul(g2, e, e);
                mpz_mul_ui(g2, g2, 2);
                mpz_mul(e, c_val, c_val);
                mpz_sub(g2, g2, e);

                // Check if all entries are perfect squares.
                if (!is_perfect_square(b2) || !is_perfect_square(d2) ||
                    !is_perfect_square(f2) || !is_perfect_square(h2) ||
                    !is_perfect_square(i2) || !is_perfect_square(g2)) {
                    continue;
                }

                // Assemble the square.
                mpz_mul(A, a_val, a_val);
                mpz_set(B, b2);
                mpz_mul(C, c_val, c_val);
                mpz_set(D_cell, d2);
                mpz_mul(E, e, e);
                mpz_set(F, f2);
                mpz_set(G, g2);
                mpz_set(H, h2);
                mpz_set(I, i2);

                mpz_t square[9];
                for(int k = 0; k < 9; k++) mpz_init(square[k]); //Should be initialized here to be thread-safe.

                mpz_set(square[0], A);
                mpz_set(square[1], B);
                mpz_set(square[2], C);
                mpz_set(square[3], D_cell);
                mpz_set(square[4], E);
                mpz_set(square[5], F);
                mpz_set(square[6], G);
                mpz_set(square[7], H);
                mpz_set(square[8], I);

                // Check for uniqueness of elements.
                bool unique = true;
                for (int k = 0; k < 9; ++k) {
                    for (int l = k + 1; l < 9; ++l) {
                        if (mpz_cmp(square[k], square[l]) == 0) {
                            unique = false;
                            break;
                        }
                    }
                    if (!unique) break;
                }

                if (unique) {
                    total_candidates++;  // Atomic increment for thread safety
                    gmp_printf("Candidate: a=%Zd, e=%Zd, c=%Zd, square=[%Zd, %Zd, %Zd, %Zd, %Zd, %Zd, %Zd, %Zd, %Zd]\n",
                           a_val, e, c_val, A, B, C, D_cell, E, F, G, H, I);
                }
                for(int k = 0; k < 9; k++) mpz_clear(square[k]);//Should be cleared in thread-safe area.

            }
        }
        // --- Checkpointing and Logging (INSIDE OpenMP critical section) ---
        #pragma omp critical
        {
            if ((d - current_D_val + 1) % log_interval == 0) {
                printf("Processed D_val: %llu\n", d); // Log current D_val
                save_checkpoint(d + 1, total_candidates);  //Save progress after each interval.
            }

        }
    } // End of the OpenMP parallel for loop.

        // --- Final Checkpointing ---
        #pragma omp critical
        {
            save_checkpoint(D_max + 1, total_candidates); //final save!
        }

    // Clean up allocated memory
     for (int i = 0; i < 4ULL * pell_limit * 2; i++) {
        mpz_clear(sols[i].N);
        mpz_clear(sols[i].M);
    }
    free(sols);

    // Clear all mpz_t variables...
    mpz_clear(D_val);
    mpz_clear(Y);
    mpz_clear(X);
    mpz_clear(V);
    mpz_clear(U);
    mpz_clear(e1);
    mpz_clear(e2);
    mpz_clear(e);
    mpz_clear(a_val);
    mpz_clear(c_val);
    mpz_clear(b2);
    mpz_clear(d2);
    mpz_clear(f2);
    mpz_clear(h2);
    mpz_clear(i2);
    mpz_clear(g2);
    mpz_clear(A);
    mpz_clear(B);
    mpz_clear(C);
    mpz_clear(D_cell);
    mpz_clear(E);
    mpz_clear(F);
    mpz_clear(G);
    mpz_clear(H);
    mpz_clear(I);

    return total_candidates;
}

int main() {
    // Initial search parameters (adjust as needed)
    unsigned long long pell_limit = 100000000;  // 100 million
    unsigned long long D_min = 199999996;     // ~200 million
    unsigned long long D_max = 3000000000;   // 2 billion

    int results = search_magic_square(pell_limit, D_min, D_max);
    printf("Total candidates: %d\n", results);
    return 0;
}
