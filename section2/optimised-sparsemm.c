#include "utils.h"

void basic_sparsemm(const COO, const COO, COO *);
void basic_sparsemm_sum(const COO, const COO, const COO,
                        const COO, const COO, const COO,
                        COO *);

/* Computes C = A*B.
 * C should be allocated by this routine.
 */
void optimised_sparsemm(const COO A, const COO B, COO *C)
{
    // Initialize C
    // Get m and n from A and B: m = A.m, n = B.n
    // Take NZ as (A.NZ + b.NZ) * 5
    alloc_sparse(A->m, B->n, (A.NZ + B.NZ)*5, C);
    // To access e.g. NZ of C, do (*C)->NZ

    // Loop over all (a,b) coordinates in A
    // For each non-zero of these, look at the coordinates in B
    // If a coordinate (b, c) is non-zero, add (a,b)*(b,c) to (a,c) in output

    // If we run out of space, use realloc_sparse(NZ, C) to increase allocated memory
    return basic_sparsemm(A, B, C);
}

/* Computes O = (A + B + C) (D + E + F).
 * O should be allocated by this routine.
 */
void optimised_sparsemm_sum(const COO A, const COO B, const COO C,
                            const COO D, const COO E, const COO F,
                            COO *O)
{
    return basic_sparsemm_sum(A, B, C, D, E, F, O);
}
