#include<stdlib.h>

void basic_gemm(int, int, int,
                const double *, int,
                const double *, int,
                double *, int);
void pack_a(int start, const double *a, const int lda);
void pack_b(int start, const int n, const double *b, const int ldb);
void micro_kernel(int loop_2, int loop_3, int loop_4, int ldc, double *c);

// Original values
// const int m_r = 4;
// const int n_r = 8;
// const int k_c = 256;
// const int m_c = 512;

// Experimental values
// int m_r = 64;
// int n_r = 8;
// int k_c = 256;
// int m_c = 512;

// Experimental 2 values - from parameter testing results
int m_r = 128;
int n_r = 8;
int k_c = 128;
int m_c = 512;


// Global variable to store the packed subarray of A
double *A_packed;
// Global variable to store the spliced subarray of A_packed
double *A_splice;
// Global variable to store the packed subarray of B
double *B_packed;
// Global variable to store the spliced subarray of B_packed
double *B_splice;

/* Compute C = C + A*B
 *
 * C has rank m x n
 * A has rank m x k
 * B has rank k x n
 * ldX is the leading dimension of the respective matrix.
 * 
 * Only works for values m, n, k such that:
 *  m % m_c == m % m_r == 0
 *  n % n_r == 0
 *  k % k_c == 0
 *
 * All matrices are stored in column major format.
 */
void optimised_gemm(int m, int n, int k,
                    const double *a, int lda,
                    const double *b, int ldb,
                    double *c, int ldc)
{
    // Allocate memory to A_packed and B_packed
    // A_packed will be k_c*m_c, B_packed will be k_c*n
    A_packed = calloc(m_c*k_c, sizeof(double));
    B_packed = calloc(k_c*n, sizeof(double));
    // Allocate memory to A_splice and B_splice
    // A_splice will be m_r*k_c doubles, B_splice will be k_c*n_r doubles
    A_splice = calloc(m_r*k_c, sizeof(double));
    B_splice = calloc(k_c*n_r, sizeof(double));

    // Split A into columns k_c wide and B into rows k_c tall. Access these simultaneously as the i'th index of column/row
    for(int loop_1 = 0; loop_1 < k; loop_1 += k_c) {
        // Pack the row from B
        int b_start = loop_1; // The start is loop_1 values down in the left-most column of B
        pack_b(b_start, n, b, ldb);
        // Split the column from A into blocks m_c tall
        for(int loop_2 = 0; loop_2 < m; loop_2 += m_c) {
            // Pack the block from A
            int a_start = loop_2 + loop_1*lda; // The start is loop_2 value down in the loop_1'th column of A
            pack_a(a_start, a, lda);
            // Split the row from B into columns n_r wide
            for(int loop_3 = 0; loop_3 < n; loop_3 += n_r) {
                // Extract this column into a new data structure
                // B_splice = memcpy(B_splice, &B_packed[loop_3*k_c], n_r*k_c*sizeof(double));
                // Only store the pointer to the start of B_splice to avoid doing memcpy
                B_splice = (B_packed + loop_3*k_c);

                // Split the block from the column from A into rows m_r tall
                for(int loop_4 = 0; loop_4 < m_c; loop_4 += m_r) {
                    // Extract this row into a new data structure
                    // A_splice = memcpy(A_splice, &A_packed[loop_4*k_c], m_r*k_c*sizeof(double));
                    // Do as aboe with B_splice
                    A_splice = (A_packed + loop_4*k_c);

                    // Multiply the row from A with the column from B store it in temp_output
                    micro_kernel(loop_2, loop_3, loop_4, ldc, c);
                }
            }
        }
    }
    // Free the allocated memory
    free(A_packed);
    free(B_packed);
    // free(A_splice);
    // free(B_splice);
}

/* Pack A from the given start index as needed for BLIS.
 * 
 * Output is returned using the global A_packed variable
 */
void pack_a(int start, const double *a, const int lda) {
    // Output is stored in A_packed

    // store which index you're inserting into
    int output_index = 0;

    // Loop over each row in the output
    for(int row = 0; row < m_c/m_r; row++) {
        // Loop over each column in this row
        for(int column = 0; column < k_c; column++) {
            // Loop over each value in this column
            for(int value = 0; value < m_r; value++) {
                // Select the appropriate value from A
                // This can be found at start + value + column*lda + row*m_r
                // start shifts everything down the A matrix, hence just added
                // value shifts downwards in each column, hence just added
                // column shifts rightwards within each row, so we add column*lda
                // row shifts m_r steps downwards, as each row is of height m_r, hence row*m_r is added
                A_packed[output_index] = a[start + value + column*lda + row*m_r];
                // increment the index
                output_index++;
            }
        }
    }
}

/* Pack B from the given start index as needed for BLIS.
 * 
 * Output is returned using the global B_packed variable
 */
void pack_b(int start, const int n, const double *b, const int ldb) {
    // Output is stored in B_packed

    // store which index you're inserting into
    int output_index = 0;

    // Loop over each column in the output
    for(int column = 0; column < n/n_r; column++) {
        // Loop over each line in this column
        for(int line = 0; line < k_c; line++) {
            // Loop over each value in this line
            for(int value = 0; value < n_r; value++) {
                // Select the appropriate value from B
                // This can be found at start + value*ldb + line + column*ldb*n_r
                // start shifts everything down the B matrix, hence just added
                // value shifts rightwards on each row, so we add value*ldb
                // line shifts downwards, hence just added
                // column shifts n_r steps rightwards, as each column is of width n_r, hence column*ldb*n_r is added
                B_packed[output_index] = b[start + value*ldb + line + column*ldb*n_r];
                // increment the index
                output_index++;
            }
        }
    }
}

/* Apply the micro kernel, i.e. matrix multiplication of A_splice * B_splice
 * 
 * All required variables for multiplication are globally defined.
 * The input variables define where in C the output should be stored:
 *      loop_2 + loop_4 tells us how far down the column we are
 *      loop_3 tells us how many columns right we are, i.e. loop_3*ldc steps in the 1D matrix
 * If called out of order may result in segmentation fault (if A_splice or B_splice are not correctly allocated memory/initialized)
 */
void micro_kernel(int loop_2, int loop_3, int loop_4, int ldc, double *c) {
    // Compute the start index
    int start = loop_2 + loop_4 + ldc*loop_3;
    // Loop over each column in A_splice and row in B_splice as pairs (i.e. column 0 row 0, column 1 row 1, etc.)
    for(int index = 0; index < k_c; index++) {
        // Loop over each value in the row of B_splice
        for(int B_val = 0; B_val < n_r; B_val++) {
            // Loop over each value in the column of A_splice
            for(int A_val = 0; A_val < m_r; A_val++) {
                // Multiply the values together, and store the result in kernel_output
                // As we loop over the row first and the column second, the output will automatically be in column-major format
                // The correct value of A_splice is found at [A_val + index*m_r], as index tracks which column we're in, and each column is m_r long. The value in B_splice is similarly found at [B_val + index*n_r], as each row is n_r long.
                c[start + A_val + B_val*ldc] += A_splice[A_val + index*m_r] * B_splice[B_val + index*n_r];
                // Start = loop_2 + loop_4 + loop_3 * ldc
                // Adjustment = A_val + B_val*ldc
            }
        }
    }
}

/* Main function used for testing as required.
 *
 * Currently set up to test A_packed (make sure the correct m_r, n_r, k_c, m_c values are set)
 */
// int main(int argc, char const *argv[])
// {
//     // Initialize variables
//     int start = 0; // Start index for packing
//     int lda = 4; // Leading dimension of A
//     int ldb = 4; // Leading dimension of B
//     int n = 4; // Width of the B matrix

//     A_packed = calloc(k_c*m_c, sizeof(double));
//     B_packed = calloc(k_c*n, sizeof(double));
//     kernel_output = calloc(m_r*n_r, sizeof(double));
//     double *b = calloc(16, sizeof(double));
//     double *a = calloc(16, sizeof(double));

//     a[0] = 0;
//     a[1] = 4;
//     a[2] = 8;
//     a[3] = 12;
//     a[4] = 1;
//     a[5] = 5;
//     a[6] = 9;
//     a[7] = 13;
//     a[8] = 2;
//     a[9] = 6;
//     a[10] = 10;
//     a[11] = 14;
//     a[12] = 3;
//     a[13] = 7;
//     a[14] = 11;
//     a[15] = 15;
//     b[0] = 0;
//     b[1] = 4;
//     b[2] = 8;
//     b[3] = 12;
//     b[4] = 1;
//     b[5] = 5;
//     b[6] = 9;
//     b[7] = 13;
//     b[8] = 2;
//     b[9] = 6;
//     b[10] = 10;
//     b[11] = 14;
//     b[12] = 3;
//     b[13] = 7;
//     b[14] = 11;
//     b[15] = 15;

//     pack_a(start, a, lda);
//     pack_b(start, n, b, ldb);
//     micro_kernel();

//     for(int i = 0; i < m_r*n_r; i++) {
//         printf("%d: %f\n", i, kernel_output[i]);
//     }
    
//     return 0;
// }
