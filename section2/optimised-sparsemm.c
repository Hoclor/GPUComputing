#include "utils.h"
#include<stdlib.h>

struct sorting_coords {
    int i, j;
    double dat;
};
typedef struct sorting_coords triplet;

void basic_sparsemm(const COO, const COO, COO *);
void basic_sparsemm_sum(const COO, const COO, const COO,
                        const COO, const COO, const COO,
                        COO *);
int cmpfunc(void *a, void *b);
void sort_COO(const COO A, COO *O);
void add_sparse_matrices(const COO A, const COO B, COO *O);


int cmpfunc(void *a, void *b) {
    triplet *tripletA = (triplet*) a;
    triplet *tripletB = (triplet*) b;

    // -1 == don't switch ('lower' to the left)

    // if A's row < B's row, return -1
    if(tripletA->i < tripletB->i) {
        return -1;
    }
    // if they're the same and A's column <= B's column, return -1
    else if(tripletA->i == tripletB->i && tripletA->j <= tripletA->j) {
        return -1;
    }
    // Otherwise, A's row > B's row, or A's row = B's row but A's column > B's column, so swap them
    return 1;
}

/* Sorts the coordinates and data of the input COO by row, then column
 * The result is returned in-place.
 */
void sort_COO(const COO A, COO *O) {
    // Create a COO to store the output in
    COO sp;
    alloc_sparse(A->m, A->n, A->NZ, &sp);

    // Create the list of sorting_coords
    struct sorting_coords *list = calloc(A->NZ, sizeof(struct sorting_coords));
    for(int i = 0; i < A->NZ; i++) {
        list[i].i = A->coords[i].i;
        list[i].j = A->coords[i].j;
        list[i].dat = A->data[i];
    }
    // Sort the list by row
    qsort(list, A->NZ, sizeof(struct sorting_coords), cmpfunc);

    // Assign the coords and data of sp from the sorted list
    for(int i = 0; i < sp->NZ; i++) {
        sp->coords[i].i = list[i].i;
        sp->coords[i].j = list[i].j;
        sp->data[i] = list[i].dat;
    }
    // Point O to sp
    *O = sp;
}

/* Computes C = A*B.
 * C should be allocated by this routine.
 */
void optimised_sparsemm(const COO unsortedA, const COO unsortedB, COO *C)
{
    // Create a temporary _p_COO object
    *C = NULL;
    COO Ctemp;
    // Create a counter for Ctemp usedNZ
    int usedNZ = 0;
    // Get m and n from A and B: m = A.m, n = B.n
    // Take NZ as (A.NZ + b.NZ) * 5
    alloc_sparse(unsortedA->m, unsortedB->n, (unsortedA->NZ + unsortedB->NZ)*5, &Ctemp);

    // sort the A and B matrices so that it is sorted by row first, then by column
    COO A;
    COO B;
    sort_COO(unsortedA, &A);
    sort_COO(unsortedB, &B);

    // Loop over all (a,b) coordinates in A
    // For each of these, look at the coordinates in B
    // For each coordinate (b, c), add (a,b)*(b,c) to (a,c) in output
    for(int a_coord_index = 0; a_coord_index < A->NZ; a_coord_index++) {
        // Get the coordinate in A
        int b_row_found = 0; 
        int A_i = A->coords[a_coord_index].i;
        int A_j = A->coords[a_coord_index].j;
        for(int b_coord_index = 0; b_coord_index < B->NZ; b_coord_index++) {
            // Check if the B.i coordinate matches the A.j coordinate (i.e. A = (a,b), B = (b,c))
            int B_i = B->coords[b_coord_index].i;
            int B_j = B->coords[b_coord_index].j;
            if(A_j == B_i) {
                b_row_found = 1;
                // Compute (a,b)*(b,c)
                double resultData = A->data[a_coord_index] * B->data[b_coord_index];
                // Loop over all output coords so far to check if we have already computed a value at (a, c)
                int coordExists = 0;
                int c_row_found = 0;
                for(int c_coord_index = 0; c_coord_index < usedNZ; c_coord_index++) {
                    // Check if this coordinate matches (a, c)
                    int C_i = Ctemp->coords[c_coord_index].i;
                    int C_j = Ctemp->coords[c_coord_index].j;
                    if(C_i == A_i) {
                        c_row_found = 1;
                        if(C_j == B_j) {
                            // Add the resultData to  (C_i, C_j)
                            Ctemp->data[c_coord_index] += resultData;
                            coordExists = 1;
                            break;
                        }
                    } else if(c_row_found) {
                        // We won't find the desired coordinates, as C is sorted by row then by column,
                        // and we have now searched all coords with the correct row
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
                    newCoord.j = B_j;
                    Ctemp->coords[usedNZ] = newCoord;
                    Ctemp->data[usedNZ] = resultData;
                    usedNZ++;
                }
            } else if(b_row_found) {
                // We won't find any other matches for this A, so skip to the next coordinate
                break;
            }
        }
    }
    // Update Ctemp's NZ to the actual used NZ
    Ctemp->NZ = usedNZ;
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

void add_sparse_matrices(const COO A, const COO B, COO *O) {
    // Create a COO for the output
    COO sp;
    // Guess that we will need A->NZ + B->NZ entries (this is the cap)
    alloc_sparse(A->m, A->n, A->NZ + B->NZ, &sp);

    // Perform a merge-sort-style merge of A and B coordinates - requires they are sorted, so sort them first
    COO sortedA;
    COO sortedB;

    sort_COO(A, &sortedA);
    sort_COO(B, &sortedB);

    // Keep a counter of the current element of A, and same for B
    int currentA = 0;
    int currentB = 0;
    // Keep a counter of which item we're inserting into sp
    int output_counter = 0;

    while(currentA < A->NZ && currentB < B->NZ) {
        // Compare the coordinates of the current A item and current B item
        // If itemA < itemB, append itemA to sp and increment currentA
        // If itemB < itemA, do the same but for B
        // If itemA == itemB, add them together, append the result, and increment both currentA and currentB
        int itemAi = sortedA->coords[currentA].i;
        int itemAj = sortedA->coords[currentA].j;
        int itemBi = sortedB->coords[currentB].i;
        int itemBj = sortedB->coords[currentB].j;

        if(itemAi < itemBi || (itemAi == itemBi && itemAj < itemBj)) {
            // Append itemA and increment currentA
            sp->coords[output_counter].i = itemAi;
            sp->coords[output_counter].j = itemAj;
            sp->data[output_counter] = sortedA->data[currentA];
            output_counter++;
            currentA++;
        // Check if they're the same
        } else if(itemBi == itemAi && itemBj == itemAj) {
            // Append the sum of their data entries
            sp->coords[output_counter].i = itemAi;
            sp->coords[output_counter].j = itemAj;
            sp->data[output_counter] = sortedA->data[currentA] + sortedB->data[currentB];
            output_counter++;
            currentA++;
            currentB++;
        } else {
            // The final case is that itemB < itemA coordinates wise, so append itemB and increment currentB
            sp->coords[output_counter].i = itemBi;
            sp->coords[output_counter].j = itemBj;
            sp->data[output_counter] = sortedB->data[currentB];
            output_counter++;
            currentB++;
        }
    }
    // Keep appending from sortedA or sortedB, whichever isn't empty
    while(currentA < A->NZ) {
        // Get the item from A
        int itemAi = sortedA->coords[currentA].i;
        int itemAj = sortedA->coords[currentA].j;
        // Input it into sp
        sp->coords[output_counter].i = itemAi;
        sp->coords[output_counter].j = itemAj;
        sp->data[output_counter] = sortedA->data[currentA];
        // increment currentA and output_counter
        currentA++;
        output_counter++;
    }
    while(currentB < B->NZ) {
        int itemBi = sortedB->coords[currentB].i;
        int itemBj = sortedB->coords[currentB].j;
        // Input it into sp
        sp->coords[output_counter].i = itemBi;
        sp->coords[output_counter].j = itemBj;
        sp->data[output_counter] = sortedB->data[currentB];
        // increment currentB and output_counter
        currentB++;
        output_counter++;
    }
    // Set sp->NZ to the actual NZ (output_counter)
    sp->NZ = output_counter;
    // return the output in the pointer *O
    *O = sp;
}