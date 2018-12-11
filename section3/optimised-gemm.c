void basic_gemm(int, int, int,
                const double *, int,
                const double *, int,
                double *, int);

const int m_r = 4;
const int n_r = 8;
const int k_c = 256;
const int m_c = 512;

// Global variable to store the packed subarray of A
double *A_packed;
// Global variable to store the packed subarray of B
double *B_packed;
// Global variable to store the output of the micro kernel (matrix multiplication)
double *kernel_output;

/* TODO
 *
 * 1. Implement B packing function
 * 2. Implement A packing function
 * 3. Rewrite basic-gemm for mix of column-major and row-major formats
 * 4. Write the correct indexing logic for each loop:
 *      - Loop 1
 *      - Loop 2
 *      - Loop 3
 *      - Loop 4
 *      - Loop 5
 *      - Loop 6
 */


/* Compute C = C + A*B
 *
 * C has rank m x n
 * A has rank m x k
 * B has rank k x n
 * ldX is the leading dimension of the respective matrix.
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
    A_packed = calloc(k_c*m_c, sizeof(double));
    B_packed = calloc(k_c*n, sizeof(double));

    // Split A into columns k_c wide. Access these simultaneously as the i'th index of column/row
    for(int loop_1 = 0; loop_1 < k; loop_1 += k_c) {
        // Pack the row from B
        //TODO pack B
        // Split the column from A into blocks m_c tall
        for(int loop_2 = 0; loop_2 < k; loop_2 += m_c) {
            // Pack the block from A
            //TODO pack A
            // Split the row from B into columns n_r wide
            for(int loop_3 = 0; loop_3 < n; n_r) {
                // Extract this column into a new data structure

                // Split the block from the column from A into rows m_r tall
                for(int loop_4 = 0; loop_4 < m_c; loop_4 += m_r) {
                    // Extract this row into a new data structure

                    // Multiply the row from A with the column from B store it in temp_output

                    // Loop over the columns of temp_output
                    for(int loop_5 = 0; loop_5 < n_r; loop_5++) {
                        // Loop over each value in these columns
                        for(int loop_6 = 0; loop_6 < m_r; loop_6++) {
                            // Move this value into the correct index of C
                            //TODO insert into C
                        }
                    }
                }
            }
        }
    }
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

/* Apply the micro kernel, i.e. matrix multiplication of A_packed * B_packed
 * 
 * All required variables are globally defined, so no inputs are required.
 * If called out of order may result in segmentation fault (if kernel_output, A_packed or B_packed are not correctly allocated memory/initialized)
 * 
 * Overwrites any previously stored values in kernel_output when called.
 */
void micro_kernel() {
    // Create a counter variable to keep track of which index of kernel_output we're inserting into
    int counter = 0;
    // Loop over each column in A_packed and row in B_packed as pairs (i.e. column 0 row 0, column 1 row 1, etc.)
    for(int index = 0; index < k_c; index++) {
        // Loop over each value in the row of B_packed
        for(int B_val = 0; B_val < n_r; B_val++) {
            // Loop over each value in the column of A_packed
            for(int A_val = 0; A_val < m_r; A_val++) {
                // Multiply the values together, and store the result in kernel_output
                // As we loop over the row first and the column second, the output will automatically be in column-major format
                // The correct value of A_packed is found at [A_val + index*m_r], as index tracks which column we're in, and each column is m_r long. The value in B_packed is similarly found at [B_val + index*n_r], as each row is n_r long.
                kernel_output[counter] = A_packed[A_val + index*m_r] * B_packed[B_val + index*n_r];
                // Increment the counter variable
                counter++;
            }
        }
    }
}

/* Main function used for testing as required.
 *
 * Currently set up to test A_packed (make sure the correct m_r, n_r, k_c, m_c values are set)
 */
int main(int argc, char const *argv[])
{
    double *b = calloc(16, sizeof(double));
    A_packed = calloc(16, sizeof(double));
    b[0] = 0;
    b[1] = 4;
    b[2] = 8;
    b[3] = 12;
    b[4] = 1;
    b[5] = 5;
    b[6] = 9;
    b[7] = 13;
    b[8] = 2;
    b[9] = 6;
    b[10] = 10;
    b[11] = 14;
    b[12] = 3;
    b[13] = 7;
    b[14] = 11;
    b[15] = 15;
    int start = 0;
    int lda = 4;

    pack_a(start, b, lda);

    for(int i = 0; i < 16; i++) {
        printf("%d: %f\n", i, A_packed[i]);
    }
    
    return 0;
}
