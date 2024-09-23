#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "three-phase-header.h"

int beta_max;
int N, K;  // node number and group number
double** Distances;   // distance matrix between elements
double** DistancesT;
int* LB; // Lower bound of the number of elements in the group i 
int* UB; // Upper bound of the number of elements in the group i
double theta, theta_max, theta_min; 
double alpha;
int beta_min; 
int eta_max;

// main processure
int maxNumberIterations;
double start_time, Time_limit;
Solution S_b; //best solution
Solution CS;
Solution *S; //S_i
Solution *O; //O_i in crossover
Neighborhood *Neighbors;

//double neighboorhood local search
int *s; // partition array for each v
double objective;

// Matrix M
double **Delta_Matrix;  // incremental matrix 
double **Delta_Matrix_p1;
double **Delta_Matrix_p2;
double *groupDiversity;
double *groupDiversity_p1;
double *groupDiversity_p2;
int *SelectEle;
int *SelectEleTemp;
int *SelectGroup;
int *p1;
int *p2;

// for crossover
int *vectorElement;
int *groupElement;
int *LBGroup;
int *UBGroup;
int *BigThanLB;
int *ub;

//directed pertubation
double** Avg;
int* Rd, * UnderLB; //Rd=R
int *SizeG; //c_g

/* Exchange Method for Anticlustering Based on a Distance matrix
 * 
 * param *distances: vector of data points (in R, this is a distance matrix,
 *         the matrix structure must be restored in C)
 * param *N_in: The number of elements (i.e., number of "rows" in *data)
 * param *K_in: The number of clusters
 * param *C: The number of categories
 * param *lower_bound: Minimum number of elements in each anticluster.
 * param *upper_bound: Maximum number of elements in each anticluster.
 * param *Beta_max: The algorithm begins with a pool of random initial solutions of size beta_max. 
 *                   Over time, the size of the solution pool decreases linearly until it reaches beta_min.
 * param *time_limit:  Maximum execution time of the algorithm (in seconds)
 * param *Theta_max: Parameter for the strength of undirected perturbation, which decreases linearly over time from theta_max to theta_min..
 * param *Theta_min: Parameter for the strength of undirected perturbation, which decreases linearly over time from theta_max to theta_min..
 * param *Beta_min: The minimum solution pool size the algorithm should reach before making a determination.
 * param *Eta_max: .
 * param *Alpha: Parameter for weitghing the discrimitation of a slighlty worse local optiomal child solution
 *               in Yang et al. set to 0.05 (might differ due to different implemetnation of calculation).
 * param *result: Calculated assignment of elements to clusters. Emptz vector.
 * param *cost: Value of objective function.
 * param *mem_error: This is passed with value 0 and only receives the value 1 
 *       if a memory error occurs when executing this function. The caller needs
 *       to test if this value is 1 after execution.
 * 
 * 
 * The return value is assigned to the argument `result`, via pointer
*/
void three_phase_search_disperion(
                      double *distances, 
                      int *N_in,
                      int *K_in, 
                      int *upper_bound, 
                      int *lower_bound, 
					  int *number_of_iterations,
                      int *clusters,
											int *Beta_max, 
											int *time_limit,
											double *Theta_max,
											double *Theta_min,
											int *Beta_min,
											int *Eta_max,
                                            double *Alpha,
											int *result,
											double *score,
											int *mem_error) {

  N = *N_in;
  K = *K_in;
  beta_max = *Beta_max;  
  theta = *Theta_max;
  theta_max = *Theta_max;
  beta_min = *Beta_min;
  eta_max = *Eta_max;
  Time_limit =  *time_limit;
  alpha = *Alpha;
  maxNumberIterations = *number_of_iterations;
  
  // Allocate memory for Distances and DistancesT arrays
  Distances = (double**)malloc(N * sizeof(double*));
  if (Distances == NULL) { *mem_error = 1; return; }
  DistancesT = (double**)malloc(N * sizeof(double*));
  if (DistancesT == NULL) { *mem_error = 1; return; }
  for (int i = 0; i < N; i++) {
    Distances[i] = (double*)malloc(N * sizeof(double));
    if (Distances[i] == NULL) { *mem_error = 1; return; }
    DistancesT[i] = (double*)malloc(N * sizeof(double));
    if (DistancesT[i] == NULL) { *mem_error = 1; return; }
  }
  
  // Fill Distances and DistancesT with values from input
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      Distances[i][j] = distances[i * N + j];
      DistancesT[i][j] = 2 * distances[i * N + j];
    }
  }
  
  if (clusters[0] == -1 ) {
    // Allocate memory for LB and UB arrays
    LB = (int*)malloc(K * sizeof(int));
    if (LB == NULL) { *mem_error = 1; return; }
    UB = (int*)malloc(K * sizeof(int));
    if (UB == NULL) { *mem_error = 1; return; }
    for (int i = 0; i < K; i++) {
        LB[i] = *lower_bound;  // Assuming lower_bound is a pointer to an int
        UB[i] = *upper_bound;  // Assuming upper_bound is a pointer to an int
    }
  } else {
    LB = (int*)malloc(K * sizeof(int));
    if (LB == NULL) { *mem_error = 1; return; }
    UB = (int*)malloc(K * sizeof(int));
    if (UB == NULL) { *mem_error = 1; return; }
    for (int i = 0; i < K; i++) {
        LB[i] = clusters[i];  // Assuming lower_bound is a pointer to an int
        UB[i] = clusters[i];  // Assuming upper_bound is a pointer to an int
    }
  }
  
  AssignMemory();
  if (*mem_error == 1) {
    return;
  }
  
  BuildNeighbors();
  
  SearchAlgorithm();
  
  //save S_b -> solution with result
  for (int i = 0; i < N; i++){
    result[i] = S_b.s[i];
  }
  *score = S_b.cost;
  
  // Remember to free the allocated memory after use
  for (int i = 0; i < N; i++) {
    free(Distances[i]); Distances[i] = NULL;
    free(DistancesT[i]); DistancesT[i] = NULL;
  }
  free(Distances); Distances = NULL;
  free(DistancesT); DistancesT = NULL;
  free(LB); LB = NULL;
  free(UB); UB = NULL;
  
  ReleaseMemory();
}



void SearchAlgorithm() {
    /* Algorithm 1: The main procedure of TPSDP. */

    int eta;
    int pickedSolution;

    //important! for windows and linux there is a differnt definition of this time
    //on windows its the wall time, on linux the CPU time
    clock_t start_time = clock();
    S_b.cost = -INFINITY;
    
    // Initial population generation
    int i, j, k;
    for (i = 0; i < beta_max; i++) {
        InitialSol(&CS);
        for (j = 0; j < N; j++) S[i].s[j] = CS.s[j];
        for (k = 0; k < K; k++) S[i].SizeG[k] = CS.SizeG[k];
        S[i].cost = CS.cost;
        if (S[i].cost > S_b.cost) {
            for (j = 0; j < N; j++) S_b.s[j] = S[i].s[j];
            for (k = 0; k < K; k++) S_b.SizeG[k] = S[i].SizeG[k];
            S_b.cost = S[i].cost;
        }
    }
 
    int counter = 1;
    while (counter <= maxNumberIterations) {
        
        eta = (int)(theta * N / K);
        for (i = 0; i < beta_max; i++) {
            for (j = 0; j < N; j++) O[i].s[j] = S[i].s[j];
            for (k = 0; k < K; k++) O[i].SizeG[k] = S[i].SizeG[k];
            O[i].cost = S[i].cost;
        }
        // Strong Perturbation and Local Search
        for (i = 0; i < beta_max; i++) {
            UndirectedPerturbation(eta, S[i].s, S[i].SizeG);
            DoubleNeighborhoodLocalSearch(S[i].s, S[i].SizeG, &S[i].cost);
            if (S[i].cost > S_b.cost) {
                for (j = 0; j < N; j++) S_b.s[j] = S[i].s[j];
                for (k = 0; k < K; k++) S_b.SizeG[k] = S[i].SizeG[k];
                S_b.cost = S[i].cost;
            }
        }

        // Crossover and Local Search
        if (beta_max > 1) {
            for (i = 0; i < beta_max; i++) {
                pickedSolution = random_int(beta_max);
                do {
                    pickedSolution = (pickedSolution + 1) % beta_max;
                } while (pickedSolution == i);
                Crossover(S[i].s, S[pickedSolution].s, O[i].s, O[i].SizeG);
                DoubleNeighborhoodLocalSearch(O[i].s, O[i].SizeG, &O[i].cost);
            }
            for (i = 0; i < beta_max; i++) {
                if (O[i].cost >= S[i].cost) {
                    for (j = 0; j < N; j++) S[i].s[j] = O[i].s[j];
                    for (k = 0; k < K; k++) S[i].SizeG[k] = O[i].SizeG[k];
                    S[i].cost = O[i].cost;
                } else if (LocalSearchCriterionCalcutlation(O[i].s, S[i].s, O[i].cost, S[i].cost) > 1) {
                    for (j = 0; j < N; j++) S[i].s[j] = O[i].s[j];
                    for (k = 0; k < K; k++) S[i].SizeG[k] = O[i].SizeG[k];
                    S[i].cost = O[i].cost;
                }
                if (S[i].cost > S_b.cost) {
                    for (j = 0; j < N; j++) S_b.s[j] = S[i].s[j];
                    for (k = 0; k < K; k++) S_b.SizeG[k] = S[i].SizeG[k];
                    S_b.cost = S[i].cost;
                }
            }
        }

        // Direct Perturbation and Local Search
        for (i = 0; i < beta_max; i++) {
            DirectPerturbation(eta_max, S[i].s, S[i].SizeG);
            DoubleNeighborhoodLocalSearch(S[i].s, S[i].SizeG, &S[i].cost);
            if (S[i].cost > S_b.cost) {
                for (j = 0; j < N; j++) S_b.s[j] = S[i].s[j];
                for (k = 0; k < K; k++) S_b.SizeG[k] = S[i].SizeG[k];
                S_b.cost = S[i].cost;
            }
        }

        // Linearly decrease population size
        // Note: Implement sort function based on the comparison function `Cmpare`
        qsort(S, beta_max, sizeof(Solution), Cmpare);
        beta_max = (int)(beta_min - beta_max) * counter / maxNumberIterations + beta_max;
        theta = theta_max - (theta_max - theta_min) * counter / maxNumberIterations;
        counter++;
    }

    // Stop measuring time and calculate the elapsed time
    clock_t end_time = clock();
    double elapsed_time = (double) (end_time - start_time)/CLOCKS_PER_SEC;
    Rprintf("The run time of the distance_clustering algortihm in seconds is: %f\n", elapsed_time);
}

/* Algorithm 2: initializes a solution S */
void InitialSol(Solution *S) {
    /* Algorithm 2: initializes a solution S */
    RandomInitialSol(S->s, S->SizeG);
    DoubleNeighborhoodLocalSearch(S->s, S->SizeG, &(S->cost));
}

void DoubleNeighborhoodLocalSearch(int partition[], int SizeGroup[], double* cost) {
    const double DELTA_THRESHOLD = 0.0001;  // Define a constant for comparison threshold
    int i, v, g, u;
    int oldGroup, oldGroup1, t;
    int imp;

    // Initialize the partition array
    for (i = 0; i < N; i++) s[i] = partition[i];

    // Build the delta_f matrix for cost changes
    BuildDeltaMatrix();

    // Initialize the delta_f value
    double delta_f = -99999.0;

    do {
        imp = 0;  // Reset improvement flag

        // First loop: Move individual elements to improve partition
        for (v = 0; v < N; v++) {
            for (g = 0; g < K; g++) {
                // Check if moving `v` to `group` is valid and beneficial
                if ((s[v] != g) && (SizeGroup[s[v]] > LB[s[v]]) && (SizeGroup[g] < UB[g])) {
                    delta_f = Delta_Matrix[v][g] - Delta_Matrix[v][s[v]];

                    if (delta_f > DELTA_THRESHOLD) {
                        oldGroup = s[v];

                        // Update delta_f matrix for the move
                        OneMoveUpdateDeltaMatrix(v, oldGroup, g);

                        // Update group sizes
                        SizeGroup[oldGroup] -= 1;
                        SizeGroup[g] += 1;

                        // Assign v to new group
                        s[v] = g;

                        // Update total cost
                        objective += delta_f;

                        // Mark as improved
                        imp = 1;
                    }
                }
            }
        }

        // Second loop: Swap pairs of elements between groups
        for (v = 0; v < N; v++) {
            for (u = v + 1; u < N; u++) {
                // Only swap if nodes are in different groups
                if (s[v] != s[u]) {
                    delta_f = (Delta_Matrix[v][s[u]] - Delta_Matrix[v][s[v]])
                          + (Delta_Matrix[u][s[v]] - Delta_Matrix[u][s[u]])
                          - DistancesT[v][u];

                    if (delta_f > DELTA_THRESHOLD) {
                        oldGroup = s[v];
                        oldGroup1 = s[u];

                        // Update delta_f matrix M for the swap
                        OneMoveUpdateDeltaMatrix(v, oldGroup, oldGroup1);
                        OneMoveUpdateDeltaMatrix(u, oldGroup1, oldGroup);

                        // Swap the two nodes between groups
                        t = s[v];
                        s[v] = s[u];
                        s[u] = t;

                        // Update total cost
                        objective += delta_f;

                        // Mark as improved
                        imp = 1;
                    }
                }
            }
        }
    } while (imp == 1);  // Continue until no improvement is made

    // Update the partition array with the final assignments
    for (i = 0; i < N; i++) partition[i] = s[i];
    *cost = objective;
}

void UndirectedPerturbation(int L, int partition[], int SizeGroup[]) {
    /* Algorithm 4: Undirected Perturbation. Applies a strong perturbation to the partition */

    int current_index;
    int v, g, x, y;
    int oldGroup, oldGroup1, swap;


    for (int i = 0; i < N; i++) {
        s[i] = partition[i];
    }

    theta = L; 
    int count = 0;
    int NumberNeighbors = N * (N - 1) / 2 + N * K;

    // Perturbation loop
    while (count < theta) {
        current_index = random_int(NumberNeighbors);

        if (Neighbors[current_index].type == 1) { // Type 1 neighbor: (element, group)
            v = Neighbors[current_index].v;
            g = Neighbors[current_index].g;

            // Apply perturbation if constraints are met
            if (s[v] != g && SizeGroup[s[v]] > LB[s[v]] && SizeGroup[g] < UB[g]) {
                oldGroup = s[v];
                SizeGroup[oldGroup]--;
                SizeGroup[g]++;
                s[v] = g;
                count++;
            }
        } else if (Neighbors[current_index].type == 2) { // Type 2 neighbor: (element x, element y)
            x = Neighbors[current_index].x;
            y = Neighbors[current_index].y;

            // Apply perturbation if elements are in different groups
            if (s[x] != s[y]) {
                oldGroup = s[x];
                oldGroup1 = s[y];
                swap = s[x];
                s[x] = s[y];
                s[y] = swap;
                count++;
            }
        }
    }

    // Copy the perturbed partition back to the original partition
    for (int i = 0; i < N; i++) {
        partition[i] = s[i];
    }
}

void DirectPerturbation(int eta_max, int partition[], int SizeGroup[]) {
    /* Algorithm 6: Directed Perturbation. 
	Iteratively refines partitions to balance group sizes and minimize costs */

    int i, j, k, L, number, minDeltaValue, minElement;

    // Initialize the partition and size groups
    for (i = 0; i < N; i++) s[i] = partition[i];
    for (j = 0; j < K; j++) SizeG[j] = SizeGroup[j];
    BuildDeltaMatrix();

    // Main loop for perturbation iterations
    for (L = 0; L < eta_max; L++) {

        // Reset tracking variables
        number = 0;
        for (i = 0; i < K; i++) {
            UnderLB[i] = 0;
            Rd[i] = -1;
            for (j = 0; j < K; j++) {
                Avg[i][j] = 0.0;
            }
        }

        // Find the minimum scoring element for each group
        for (k = 0; k < K; k++) {
            minDeltaValue = 99999999;
            minElement = 0;
            for (i = 0; i < N; i++) {
                if (s[i] == k) {
                    if (Delta_Matrix[i][k] < minDeltaValue) {
                        minDeltaValue = Delta_Matrix[i][k];
                        minElement = i;
                    }
                }
            }

            // Record the minimum element for removal
            Rd[k] = minElement;
            SizeG[k] -= 1;

            // If the group size falls below the lower bound, mark it
            if (SizeG[k] < LB[k]) {
                UnderLB[k] = 1;
                number += 1;
            }
        }

        // Rebuild the Delta matrix and average connections after removal
        for (i = 0; i < K; i++) {
            for (j = 0; j < K; j++) {
                Delta_Matrix[Rd[i]][s[Rd[j]]] = Delta_Matrix[Rd[i]][s[Rd[j]]] - Distances[Rd[i]][Rd[j]];
                Avg[s[Rd[i]]][s[Rd[j]]] = Delta_Matrix[Rd[i]][s[Rd[j]]] / SizeG[s[Rd[j]]];
            }
        }        
		
        // Handle groups that are under the lower bound (LB)
        int selectedGroup;
        int maxAvgCon;
        int nn = 0;
        int i;
        while (nn < number) {
            maxAvgCon = -9999;
            i = random_int(K);

            // Find the element with the highest average connection to the group
            do {
                i = (i + 1) % K;
            } while (UnderLB[i] == 0);
            for (j = 0; j < K; j++) {
                if (Avg[j][i] > maxAvgCon) {
                    maxAvgCon = Avg[j][i];
                    selectedGroup = j;
                }
            }

            // Move the selected element to the group
            SizeG[i] += 1;
            for (k = 0; k < K; k++) {
                if (Rd[k] != -1) {
                    Delta_Matrix[Rd[k]][i] += Distances[Rd[k]][Rd[selectedGroup]];
                    Avg[s[Rd[k]]][i] = Delta_Matrix[Rd[k]][i] / SizeG[i];
                }
            }

            // Clear old connections for the moved element and finalize the move
            for (k = 0; k < K; k++) {
                Avg[s[Rd[selectedGroup]]][k] = 0.0;
            }
            s[Rd[selectedGroup]] = i;
            UnderLB[i] = 0;
            Rd[selectedGroup] = -1;
            nn++;
        }

        // Handle remaining elements for groups that are above LB
        int groupWithMaxAvgCon;
        nn = 0;
        while (nn < K - number) {
            selectedGroup = rand() % K;
            do {
                selectedGroup = (selectedGroup + 1) % K;
            } while (Rd[selectedGroup] == -1);
            maxAvgCon = -9999;
            for (j = 0; j < K; j++) {
                if (Avg[selectedGroup][j] > maxAvgCon) {
                    maxAvgCon = Avg[selectedGroup][j];
                    groupWithMaxAvgCon = j;
                }
            }
            // Move the selected element to the group
            if (SizeG[groupWithMaxAvgCon] < UB[groupWithMaxAvgCon]) {
                SizeG[groupWithMaxAvgCon] += 1;
                for (k = 0; k < K; k++) {
                    if (Rd[k] != -1) {
                        Delta_Matrix[Rd[k]][groupWithMaxAvgCon] += Distances[Rd[k]][Rd[selectedGroup]];
                        Avg[s[Rd[k]]][groupWithMaxAvgCon] = Delta_Matrix[Rd[k]][groupWithMaxAvgCon] / SizeG[groupWithMaxAvgCon];
                    }
                }
                for (k = 0; k < K; k++) {
                	Avg[s[Rd[selectedGroup]]][k] = 0.0;
                }
                s[Rd[selectedGroup]] = groupWithMaxAvgCon;
				Rd[selectedGroup] = -1;
				nn += 1;
			}
			else {
				for (k = 0; k < K; k++) {
					Avg[k][groupWithMaxAvgCon] = 0.0;
				}
            }
        }
        BuildDeltaMatrix();
    }

    for (i = 0; i < N; i++) partition[i] = s[i];
    for (j = 0; j < K; j++) SizeGroup[j] = SizeG[j];
}

void Crossover(int partition1[], int partition2[], int score[], int scSizeGroup[]) {
    /* Algorithm 5: combines partitions in a way that maintains group constraints */

    int i, j, maxGroupDiversity, selectedGroup;
    int elementCount, groupCount;
    int targetGroup = -1;
    int processedCount;
    int selectedElement;
    int totalLowerBound, totalBelowLowerBound;

    // Initialize s and p1 with partition1
    for (i = 0; i < N; i++) {
        s[i] = partition1[i];
        p1[i] = partition1[i];
    }
    BuildDeltaMatrix();
    for (i = 0; i < N; i++) {
        for (j = 0; j < K; j++) {
            Delta_Matrix_p1[i][j] = Delta_Matrix[i][j];
        }
    }
    BuildGroupDiversityForCrossover();
    for (i = 0; i < K; i++) {
        groupDiversity_p1[i] = groupDiversity[i];
    }

    // Initialize s and p2 with partition2
    for (i = 0; i < N; i++) {
        s[i] = partition2[i];
        p2[i] = partition2[i];
    }
    BuildDeltaMatrix();
    for (i = 0; i < N; i++) {
        for (j = 0; j < K; j++) {
            Delta_Matrix_p2[i][j] = Delta_Matrix[i][j];
        }
    }
    BuildGroupDiversityForCrossover();
    for (i = 0; i < K; i++) {
        groupDiversity_p2[i] = groupDiversity[i];
    }

    // Initialize arrays
    for (i = 0; i < N; i++) {
        vectorElement[i] = i;
        score[i] = -1;
    }
    for (i = 0; i < K; i++) {
        LBGroup[i] = 0;
        UBGroup[i] = 0;
        BigThanLB[i] = 0;
        groupElement[i] = i;
        ub[i] = UB[i];
        scSizeGroup[i] = 0;
    }

    // Main crossover process
    for (i = 0; i < K; i++) {
        if (uniform_rnd_number() < 0.5) {
            // Process partition1
            maxGroupDiversity = -9999;
            for (j = 0; j < K; j++) {
                if (groupDiversity_p1[j] > maxGroupDiversity) {
                    maxGroupDiversity = groupDiversity_p1[j];
                    selectedGroup = j;
                }
            }

            elementCount = 0;
            for (j = 0; j < N; j++) {
                if (p1[j] == selectedGroup) {
                    SelectEle[elementCount++] = j;
                }
            }

            groupCount = 0;
            for (j = 0; j < K; j++) {
                if (ub[j] != -1 && ub[j] >= elementCount) {
                    SelectGroup[groupCount++] = j;
                }
            }

            if (groupCount == 0) { // No valid group found
                int minDiff = 999999;
                for (j = 0; j < K; j++) {
                    if (ub[j] != -1 && elementCount - ub[j] < minDiff) {
                        minDiff = elementCount - ub[j];
                        targetGroup = j;
                    }
                }

                processedCount = 0;
                while (processedCount < elementCount - minDiff) {
                    selectedElement = random_int(elementCount);
                    do {
                        selectedElement = (selectedElement + 1) % elementCount;
                    } while (SelectEle[selectedElement] == -1);

                    score[SelectEle[selectedElement]] = targetGroup;
                    SelectEleTemp[processedCount++] = SelectEle[selectedElement];
                    vectorElement[SelectEle[selectedElement]] = -1;
                    SelectEle[selectedElement] = -1;
                }
                elementCount = processedCount;
            } else {
                targetGroup = SelectGroup[random_int(groupCount)];
                for (j = 0; j < elementCount; j++) {
                    score[SelectEle[j]] = targetGroup;
                    vectorElement[SelectEle[j]] = -1;
                    SelectEleTemp[j] = SelectEle[j];
                }
            }
        } else {
            // Process partition2 (similar to partition1 logic)
            maxGroupDiversity = -9999;
            for (j = 0; j < K; j++) {
                if (groupDiversity_p2[j] > maxGroupDiversity) {
                    maxGroupDiversity = groupDiversity_p2[j];
                    selectedGroup = j;
                }
            }

            elementCount = 0;
            for (j = 0; j < N; j++) {
                if (p2[j] == selectedGroup) {
                    SelectEle[elementCount++] = j;
                }
            }

            groupCount = 0;
            for (j = 0; j < K; j++) {
                if (ub[j] != -1 && ub[j] >= elementCount) {
                    SelectGroup[groupCount++] = j;
                }
            }

            if (groupCount == 0) { // No valid group found
                int minDiff = 999999;
                for (j = 0; j < K; j++) {
                    if (ub[j] != -1 && elementCount - ub[j] < minDiff) {
                        minDiff = elementCount - ub[j];
                        targetGroup = j;
                    }
                }

                processedCount = 0;
                while (processedCount < elementCount - minDiff) {
                    selectedElement = random_int(elementCount);
                    do {
                        selectedElement = (selectedElement + 1) % elementCount;
                    } while (SelectEle[selectedElement] == -1);

                    score[SelectEle[selectedElement]] = targetGroup;
                    SelectEleTemp[processedCount++] = SelectEle[selectedElement];
                    vectorElement[SelectEle[selectedElement]] = -1;
                    SelectEle[selectedElement] = -1;
                }
                elementCount = processedCount;
            } else {
                targetGroup = SelectGroup[random_int(groupCount)];
                for (j = 0; j < elementCount; j++) {
                    score[SelectEle[j]] = targetGroup;
                    vectorElement[SelectEle[j]] = -1;
                    SelectEleTemp[j] = SelectEle[j];
                }
            }
        }

        // Update group diversity
        for (j = 0; j < elementCount; j++) {
            groupDiversity_p1[p1[SelectEleTemp[j]]] -= Delta_Matrix_p1[SelectEleTemp[j]][p1[SelectEleTemp[j]]];
            groupDiversity_p2[p2[SelectEleTemp[j]]] -= Delta_Matrix_p2[SelectEleTemp[j]][p2[SelectEleTemp[j]]];
            p1[SelectEleTemp[j]] = -1;
            p2[SelectEleTemp[j]] = -1;
        }

        ub[targetGroup] = -1;
        scSizeGroup[targetGroup] = elementCount;
    }

    // Adjust assignments to maintain group size constraints
    processedCount = 0;
    totalLowerBound = 0;
    totalBelowLowerBound = 0;
    for (i = 0; i < K; i++) {
        totalLowerBound += LB[i];
        if (scSizeGroup[i] < LB[i]) {
            processedCount += scSizeGroup[i];
            totalBelowLowerBound += scSizeGroup[i];
            LBGroup[i] = 1;
        } else {
            processedCount += LB[i];
        }
        if (scSizeGroup[i] > LB[i]) {
            BigThanLB[i] = 1;
        }
    }

    // Assign unprocessed elements
    for (i = 0; i < N; i++) {
        if (vectorElement[i] != -1) {
            processedCount++;
        }
    }

    while (processedCount < totalLowerBound) {
        targetGroup = random_int(K);
        do {
            targetGroup = (targetGroup + 1) % K;
        } while (BigThanLB[targetGroup] == 0);

        elementCount = 0;
        for (j = 0; j < N; j++) {
            if (score[j] == targetGroup) {
                SelectEle[elementCount++] = j;
            }
        }

        selectedElement = random_int(elementCount);
        score[SelectEle[selectedElement]] = -1;
        vectorElement[SelectEle[selectedElement]] = SelectEle[selectedElement];
        scSizeGroup[targetGroup]--;
        if (scSizeGroup[targetGroup] == LB[targetGroup]) {
            BigThanLB[targetGroup] = 0;
        }
        processedCount++;
    }

    int sumLB = 0;
    for (i = 0; i < K; i++) {
        if (LBGroup[i] == 1) {
            sumLB += LB[i];
        }
    }

    while (totalBelowLowerBound < sumLB) {
        targetGroup = random_int(K);
        do {
            targetGroup = (targetGroup + 1) % K;
        } while (LBGroup[targetGroup] == 0);

        elementCount = 0;
        for (i = 0; i < N; i++) {
            if (vectorElement[i] != -1) {
                SelectEle[elementCount++] = i;
            }
        }

        selectedElement = random_int(elementCount);
        score[SelectEle[selectedElement]] = targetGroup;
        vectorElement[SelectEle[selectedElement]] = -1;
        scSizeGroup[targetGroup]++;
        if (scSizeGroup[targetGroup] == LB[targetGroup]) {
            LBGroup[targetGroup] = 0;
        }
        totalBelowLowerBound++;
    }

    int totalSize = 0;
    for (i = 0; i < K; i++) {
        totalSize += scSizeGroup[i];
        if (scSizeGroup[i] < UB[i]) {
            UBGroup[i] = 1;
        }
    }

    while (totalSize < N) {
        targetGroup = random_int(K);
        do {
            targetGroup = (targetGroup + 1) % K;
        } while (UBGroup[targetGroup] == 0);

        elementCount = 0;
        for (i = 0; i < N; i++) {
            if (vectorElement[i] != -1) {
                SelectEle[elementCount++] = i;
            }
        }

        selectedElement = random_int(elementCount);
        score[SelectEle[selectedElement]] = targetGroup;
        vectorElement[SelectEle[selectedElement]] = -1;
        scSizeGroup[targetGroup]++;
        if (scSizeGroup[targetGroup] == UB[targetGroup]) {
            UBGroup[targetGroup] = 0;
        }
        totalSize++;
    }
}

double LocalSearchCriterionCalcutlation(int partition1[], int partition2[], double cost1, double cost2) {
    /* 
     * Evaluates the quality and dissimilarity of partitions.
     * It calculates the value that combines the ratio of costs and a
     * dissimilarity factor between partition1 and partition2.
     */
    
    // Handle potential division by zero
    if (cost2 == 0.0) {
        fprintf(stderr, "Error: Division by zero (cost2 is zero).\n");
        return -1;  
    }

    int i, j;
    int count = 0;
    int totalPairs = (N * (N - 1)) / 2;  // Number of unique pairs (i, j) where i < j

    // Loop over all pairs of elements to count dissimilar pairs
    for (i = 0; i < N - 1; i++) {
        for (j = i + 1; j < N; j++) {
            // Count cases where elements are grouped differently in the two partitions
            if ((partition1[i] == partition1[j]) !=  (partition2[i] == partition2[j])) {
                count++;
            }
        }
    }

    // Calculate ratio of costs + dissimilarity factor
    // should alpha not be alpha/2 ?
    double dissimilarityFactor = ((double)count / totalPairs) * K;
    return  cost1 / cost2 +  alpha *  dissimilarityFactor;
}

void BuildNeighbors() {
	/* Initializes the neighbor structure for optimization */
	int i, j;
	int count = 0;
	
    // Type 1 neighbors: (i, j) where each element i can be in group j
	for (i = 0; i < N; i++)
		for (j = 0; j < K; j++)
		{
			Neighbors[count].type = 1;
			Neighbors[count].v = i;
			Neighbors[count].g = j;
			count++;
		}

    // Type 2 neighbors: (i, j) where each pair of elements (i, j) are neighbors    
	for (i = 0; i < N; i++)
		for (j = i + 1; j < N; j++)
		{
			Neighbors[count].type = 2;
			Neighbors[count].x = i;
			Neighbors[count].y = j;
			count++;
		}
		
}

void AssignMemory() {
    /*  Allocates memory dynamically for various arrays and matrices necessary 
	for the algorithm's execution. This includes structures for population management, 
	distance matrices, diversity measures, and neighborhood exploration.
	*/
    
    s = (int*)malloc(N * sizeof(int));
    SizeG = (int*)malloc(K * sizeof(int));
    
    S = (Solution*)malloc(beta_max * sizeof(Solution));
    O = (Solution*)malloc(beta_max * sizeof(Solution));
    
    int i;
    groupDiversity = (double*)malloc(K * sizeof(double));
    groupDiversity_p1 = (double*)malloc(K * sizeof(double));
    groupDiversity_p2 = (double*)malloc(K * sizeof(double));
    
    for (i = 0; i < beta_max; i++) {
        S[i].s = (int*)malloc(N * sizeof(int));
        O[i].s = (int*)malloc(N * sizeof(int));
    }
    
    for (i = 0; i < beta_max; i++) {
        S[i].SizeG = (int*)malloc(K * sizeof(int));
        O[i].SizeG = (int*)malloc(K * sizeof(int));
    }
    
    CS.s = (int*)malloc(N * sizeof(int));
    S_b.s = (int*)malloc(N * sizeof(int));
    
    CS.SizeG = (int*)malloc(K * sizeof(int));
    S_b.SizeG = (int*)malloc(K * sizeof(int));
    
    Neighbors = (Neighborhood*)malloc((N * (N - 1) / 2 + N * K) * sizeof(Neighborhood));
    
    Avg = (double**)malloc(K * sizeof(double*));
    for (i = 0; i < K; i++) Avg[i] = (double*)malloc(K * sizeof(double));
    Rd = (int*)malloc(K * sizeof(int));
    for (i = 0; i < K; i++) Rd[i] = 0;
    UnderLB = (int*)malloc(K * sizeof(int));
    
    ub = (int*)malloc(K * sizeof(int));
    LBGroup = (int*)malloc(K * sizeof(int));
    UBGroup = (int*)malloc(K * sizeof(int));
    BigThanLB = (int*)malloc(K * sizeof(int));
    vectorElement = (int*)malloc(N * sizeof(int));
    groupElement = (int*)malloc(K * sizeof(int));
    SelectEle = (int*)malloc(N * sizeof(int));
    SelectGroup = (int*)malloc(K * sizeof(int));
    SelectEleTemp = (int*)malloc(N * sizeof(int));
    p1 = (int*)malloc(N * sizeof(int));
    p2 = (int*)malloc(N * sizeof(int));
}

void ReleaseMemory() {
    /* responsible for reading the input file, 
    initializing matrices, and setting constraints on group sizes. */ 
    
    free(s); s = NULL;
    free(SizeG); SizeG = NULL;
    
    free(CS.s); CS.s = NULL;
    free(CS.SizeG); CS.SizeG = NULL;
    free(S_b.s); S_b.s = NULL;
    free(S_b.SizeG); S_b.SizeG = NULL;
    
    free(LB); LB = NULL;
    free(UB); UB = NULL;
    free(Neighbors); Neighbors = NULL;
    
    free(groupDiversity); groupDiversity = NULL;
    free(groupDiversity_p1); groupDiversity_p1 = NULL;
    free(groupDiversity_p2); groupDiversity_p2 = NULL;
    free(Avg); Avg = NULL;
    free(Rd); Rd = NULL;
    free(UnderLB); UnderLB = NULL;
    free(ub); ub = NULL;
    free(LBGroup); LBGroup = NULL;
    free(UBGroup); UBGroup = NULL;
    free(BigThanLB); BigThanLB = NULL;
    free(vectorElement); vectorElement = NULL;
    free(groupElement); groupElement = NULL;
    free(SelectEle); SelectEle = NULL;
    free(SelectGroup); SelectGroup = NULL;
    free(SelectEleTemp); SelectEleTemp = NULL;
    free(p1); p1 = NULL;
    free(p2); p2 = NULL;
}
