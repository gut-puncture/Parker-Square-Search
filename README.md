# Finding a Parker Square

This project aims to find a 3x3 magic square where every entry is a distinct integer perfect square. Such a square is called a Parker Square (context: https://youtu.be/aOT_bG-vWyg?si=EMINyWgpb-p7gHZC).
This is a classic, unsolved problem in recreational mathematics.  Instead of brute-forcing combinations of numbers, we use a mathematical approach based on the Pell equation, significantly reducing the search space.

## The Problem

A 3x3 magic square is an arrangement of nine distinct numbers in a square grid such that the sum of the numbers in each row, column, and diagonal is the same.  The "magic constant" is the value of this sum.  In a magic square of *squares*, all nine numbers are also perfect squares.
A B C
D E F
G H I
We are looking for integer solutions to this problem where A, B, C, D, E, F, G, H, and I are all perfect squares, and all the row, column, and diagonal sums are equal.  Furthermore, all nine numbers must be *distinct*.

## The Methodology - Short Version
The 9 numbers of the square can be reduced to three variables: a,e and c. Then these variables are again written in terms of X,Y,U and V. This is done so that we can write these two equations:
Y² - 2X² = D
V² - 2U² = D
This is done for two reasons:
1. Both these equations are a version of the standard "Pell Equation" and the solutions to all Pell Equations can be algorithmically generated.
2. The first point means that for a random value of D, I can find X,Y,U and V and with these 4 variables, I can find a,e and c and with these 3 variables I can find the all 9 numbers of the Parker square.

This mathematical reduction means that I only have to iterate over 1 variable and can generate all others, with the constraints of the Parker Square still encoded.


## The Methodology - Long Version
The key insight is to relate the elements of the magic square to solutions of a specific Diophantine equation – the Pell equation. This approach drastically reduces the search space compared to a brute-force method.

1.  **Defining `a`, `e`, and `c`:**

    We define:

    *   `e² = E` (the center element of the magic square)
    *   `a² = A` (the top-left element)
    *   `c² = C` (the top-right element)

    Crucially, `a`, `e`, and `c` are integers (since A, E, and C are perfect squares).

2.  **Deriving Relationships:**

    It can be shown (see references below for the full derivation) that the other elements of the magic square can be expressed in terms of `a`, `e`, and `c`:

    *   `B = b² = 3e² - a² - c²`
    *   `D = d² = e² - a² + c²`
    *   `F = f² = e² + a² - c²`
    *   `G = g² = 2e² - c²`
    *   `H = h² = a² + c² - e²`
    *   `I = i² = 2e² - a²`

    This is a significant reduction. Instead of searching for nine independent numbers, we only need to find three: `a`, `e`, and `c`.

3.  **The Magic Constant and `e`:**

    The magic constant (the sum of each row, column, and diagonal) can be shown to be `3e²`.  This further reinforces the central role of `e`.

4.  **Introducing `X`, `Y`, `U`, and `V`:**
    We can express `a`, `c`, and `e` in relation to the solutions of two Pell-like equations.

    Let's consider the following two equations:
     Y² - 2X² = D
     V² - 2U² = D

    *   If we find solutions (Y, X) and (V, U) to these equations, we can define:
        *  `e = (Y² + 2X²)/2 = (V² + 2U²)/2`
        *   `a = (Y² - 2X²)/2`
        *  `c = (V² - 2U²)/2`

    * The value 'D' becomes a *parameter* of our search. It is the common value of the two Pell-like equations.

5. **Why the Pell Equation?**
   The standard Pell equation is of the form x² - Dy² = 1. The equation we have above is called a Pell-like equation. It can be proven that solutions to Pell-like equations can be methodically generated. The reason to reduce it to this form is that methods exist to find solutions to these type of equations.

6.  **The Search Algorithm:**

    Our algorithm works as follows:

    *   **Loop over `D`:**  We iterate through a range of possible values for `D` (from `D_min` to `D_max`).
    *   **Solve Pell-like Equations:** For each `D`, we find solutions (Y, X) and (V, U) to the equations Y² - 2X² = D and V² - 2U² = D, up to a specified limit (`pell_limit`). We use a modified, brute-force search for the solutions.  For each solution (N,M), we add to a solution set (-N,M), (N,-M) and (-N,-M).
    *   **Generate Candidate Magic Squares:**  For each *pair* of solutions (Y, X) and (V, U), we calculate `a`, `e`, and `c`.  We then calculate the other elements of the potential magic square (b², d², f², g², h², i²) using the formulas in step 2.
    *   **Check for Validity:** We check two critical conditions:
        *   All the calculated values (b², d², f², g², h², i²) must be perfect squares.
        *   All nine elements of the magic square (A, B, C, D, E, F, G, H, I) must be distinct.
    *   **Output:** If a candidate passes both checks, we've found a valid magic square of squares, and we print the results.

7. **Checkpointing**
    The search space is very large. Because the program may be interrupted, we implement a checkpoint system that saves the progress in a file named `magic_square_checkpoint.txt` which is loaded when the search restarts.

## Code Implementation

The code is implemented in C using the GMP library for arbitrary-precision arithmetic. This is necessary because the numbers involved can become extremely large.  The code uses OpenMP for CPU parallelization to speed up the search.

## Running the Code

1.  **Install Dependencies:**  You'll need the GMP library and an OpenMP-enabled C compiler (typically GCC). On a Debian/Ubuntu system (like Colab):

    ```bash
    sudo apt-get update
    sudo apt-get install -y libgmp-dev libomp-dev
    ```

2.  **Compile:** Compile the code using GCC with the `-lgmp` and `-fopenmp` flags:

    ```bash
    gcc -o magic_square_gmp_omp magic_square_gmp_omp.c -lgmp -fopenmp
    ```

3.  **Run:** Execute the compiled program:

    ```bash
    ./magic_square_gmp_omp
    ```

4.  **Parameters:**  The program takes three parameters (although they are currently hardcoded within `main()`):

    *   `pell_limit`:  The maximum value for `X`, `Y`, `U`, and `V` in the Pell equation solutions. This is the most crucial parameter for expanding the search space.
    *   `D_min`:  The minimum value for `D` to consider.  This is related to the lower bound on the size of the magic square numbers.  For squares with elements > 10^16, `D_min` should be around 2 * 10^8.
    *   `D_max`:  The maximum value for `D` to consider.

## Optimization

The code is parallelized using OpenMP, distributing the search over multiple CPU cores.  This significantly reduces the runtime. The checkpoint system described above ensures that the progress is saved regularly.

## Limitations

*   The code uses a brute-force search for Pell equation solutions. While much more efficient than a brute-force search for the magic square itself, this is still computationally expensive.
*   The code cannot directly utilize TPUs.  A complete rewrite in a framework like TensorFlow or JAX would be necessary, but is likely impractical due to the nature of the problem (arbitrary-precision integer arithmetic and branching).

## Future Improvements

*   Explore more efficient algorithms for finding Pell equation solutions (e.g., using continued fractions).  This could significantly speed up the search.
*   Investigate further mathematical properties of magic squares of squares to potentially prune the search space even more.
