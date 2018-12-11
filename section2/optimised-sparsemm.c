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
    // Create a temporary _p_COO object
    *C = NULL;
    COO Ctemp;
    // Create a counter for Ctemp usedNZ
    int usedNZ = 0;
    // Get m and n from A and B: m = A.m, n = B.n
    // Take NZ as (A.NZ + b.NZ) * 5
    alloc_sparse(A->m, B->n, (A->NZ + B->NZ)*5, &Ctemp);

    // Loop over all (a,b) coordinates in A
    // For each of these, look at the coordinates in B
    // For each coordinate (b, c), add (a,b)*(b,c) to (a,c) in output
    for(int a_coord_index = 0; a_coord_index < A->NZ; a_coord_index++) {
        // Get the coordinate in A
        int A_i = A->coords[a_coord_index].i;
        int A_j = A->coords[a_coord_index].j;
        for(int b_coord_index = 0; b_coord_index < B->NZ; b_coord_index++) {
            // Check if the B coordinate matches the A coordinate (i.e. A = (a,b), B = (b,c))
            int B_i = B->coords[a_coord_index].i;
            int B_j = B->coords[a_coord_index].j;
            if(A_j == B_i) {
                // Compute (a,b)*(b,c)
                double resultData = A->data[a_coord_index]*B->data[b_coord_index];
                // Loop over all output coords so far to check if we have already computed a value at (a, c)
                int coordExists = 0;
                for(int c_coord_index = 0; c_coord_index < usedNZ; c_coord_index++) {
                    // Check if this coordinate matches (a, c)
                    int C_i = Ctemp->coords[c_coord_index].i;
                    int C_j = Ctemp->coords[c_coord_index].j;
                    if(C_i == A_i && C_j == B_j) {
                        // Add the resultData to  (C_i, C_j)
                        Ctemp->data[c_coord_index] += resultData;
                        coordExists = 1;
                        break;
                    }
                }
                // Check if we found the coordinate
                if(!coordExists) {
                    // If we did not, first check if we've used up the available memory
                    if(usedNZ == Ctemp->NZ) {
                        // Reallocate Ctemp to NZ + A.NZ + B.NZ
                        realloc_sparse(Ctemp->NZ + A->NZ + B->NZ, &Ctemp);
                    }
                    // If we did not, add a new coord to C and increment usedNZ
                    struct coord newCoord;
                    newCoord.i = A_i;
                    newCoord.j = A_j;
                    Ctemp->coords[usedNZ] = newCoord;
                    Ctemp->data[usedNZ] = resultData;
                    usedNZ++;
                }
            }
        }
    }
    // Multiplication is done, return the result in C
    *C = Ctemp;
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
