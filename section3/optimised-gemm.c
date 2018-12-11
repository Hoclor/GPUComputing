void basic_gemm(int, int, int,
                const double *, int,
                const double *, int,
                double *, int);

const int m_r = 4;
const int n_r = 8;
const int k_c = 256;
const int m_c = 512;
const double *A_packed;

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
