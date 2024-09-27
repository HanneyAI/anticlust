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

/* TPSPD for Anticlustering Based on a Distance matrix
 * 
 * param *distannces: vector of data points (in R, this is a distance matrix,
 *         the matrix structure must be restored in C)
 * param *N_in: The number of elements (i.e., number of "rows" in *data)
 * param *K_in: The number of clusters. When lower_bound and upper_boound are set to the number of K, 
 *              the clusters will be equally sized. 
 * param *number_of_iterations: A number that defines how many times the steps in the search algorithm are repeated.
 * param *clusters: A predefined vector of length K specifies the number of elements in each cluster.
 *               If a default vector [-1] is provided, cluster sizes will be determined based 
 *               on the lower and upper bounds. When a cluster size array is provided, 
 *               the lower and upper bounds are ignored as they become redundant.
 * param *lower_bound: Minimum number of elements in each anticluster. 
 * param *upper_bound: Maximum number of elements in each anticluster.
 * param *Beta_max: The algorithm begins with a pool of random initial solutions of size beta_max. 
 *                   Over time, the size of the solution pool decreases linearly until it reaches beta_min.
 * param *elapsed_time: Measures the runtime of the algotihm (in seconds)
 * param *Theta_max: Parameter for the strength of undirected perturbation,
 *                   which decreases linearly over time from theta_max to theta_min..
 * param *Theta_min: Parameter for the strength of undirected perturbation, 
 *                   which decreases linearly over time from theta_max to theta_min..
 * param *Beta_min: The minimum solution pool size the algorithm should reach before making a determination.
 * param *Eta_max: A parameter that specifies how many times the steps in the direct perturbation are executed.
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
void three_phase_search_dynamic_population_size(
                      double *distances, 
                      int *N_in,
                      int *K_in, 
                      int *number_of_iterations,
                      int *clusters,
                      int *upper_bound, 
                      int *lower_bound, 
											int *Beta_max, 
											int *elapsed_time,
											double *Theta_max,
											double *Theta_min,
											int *Beta_min,
											int *Eta_max,
                                            double *Alpha,
											int *result,
											double *score,
											int *mem_error) {
  Rprintf("C is in diverity");
  N = *N_in;
  K = *K_in;
  beta_max = *Beta_max;  
  theta = *Theta_max;
  theta_max = *Theta_max;
  beta_min = *Beta_min;
  eta_max = *Eta_max;
  Time_limit =  *elapsed_time;
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
  *elapsed_time = Time_limit;
  
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
    Time_limit = elapsed_time;
}

/* Algorithm 2: initializes a solution S */
void InitialSol(Solution *S) {
    /* Algorithm 2: initializes a solution S */
    RandomInitialSol(S->s, S->SizeG);
    DoubleNeighborhoodLocalSearch(S->s, S->SizeG, &(S->cost));
}

int Cmpare(const void *a, const void *b) {
    /* Compares two solutions based on their cost */
    Solution *solA = (Solution *)a;
    Solution *solB = (Solution *)b;
    return (solB->cost - solA->cost); // Return positive if b is greater, negative if a is greater
}

// Function to swap two elements
void swap_elements(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}


void fisher_yates_shuffle(int arr[], int n) {
    // Fisher-Yates shuffle to generate a random permutation 
    for (int i = n - 1; i > 0; i--) {
        // Generate a random index j where 0 <= j <= i
        int j = random_int(i + 1);

        // Swap arr[i] and arr[j]
        swap_elements(&arr[i], &arr[j]);
    }
}

void RandomInitialSol(int s[], int SizeG[]) {
	/* Algorithm 2: initializes a random solution that respects the group size constraints */

    // Allocates memory
	int* isAssigned = (int *)malloc(N * sizeof(int));  // Tracks if an element is assigned
	int* groupSize =(int *)malloc(K * sizeof(int)); // Stores the size of each group
    int* permutedIndexList =(int *)malloc(N * sizeof(int)); // Stores the size of each group
    int* permutedGroupList =(int *)malloc(K * sizeof(int)); // Stores the size of each group
	
    int i;
	for (i = 0; i < K; i++) groupSize[i] = 0;
	for (i = 0; i < N; i++) isAssigned[i] = 0;
    for (i = 0; i < N; i++) permutedIndexList[i] = i;
    for (i = 0; i < K; i++) permutedGroupList[i] = i;

    fisher_yates_shuffle(permutedIndexList, N);
    fisher_yates_shuffle(permutedGroupList, K);

    // Calculate the total number of elements that need to satisfy the lower bounds
    int total_assigned = 0;
    int total_LB = 0;
    for (int i = 0; i < K; i++) {
        total_LB += LB[i];
    }

    // First phase: Assign elements to satisfy lower bound constraints (LB)
    int selected_element = 0;
    while (total_assigned < total_LB) {

        for (int group = 0; group < K; group++) {
            if (groupSize[group] < LB[group]) {
                s[permutedIndexList[selected_element]] = group;
                isAssigned[permutedIndexList[selected_element]] = 1;
                groupSize[group]++;
                total_assigned++;
                break;  // Move to the next element once assigned
            }
        }
        selected_element++;
    }

	// Second phase: Assign remaining elements, respecting the upper bound (UB)
    while (total_assigned < N) {
        for (int group = 0; group < K; group++) {
            if (groupSize[permutedGroupList[group]] < UB[permutedGroupList[group]]) {
                s[permutedIndexList[selected_element]] = permutedGroupList[group];
                isAssigned[permutedIndexList[selected_element]] = 1;
                groupSize[permutedGroupList[group]]++;
                total_assigned++;
                break;
            }
        }
        fisher_yates_shuffle(permutedGroupList, K);
        selected_element++;
    }

    // Copy the final group sizes into the output array SizeG
   	for (i = 0; i < K; i++)  SizeG[i] = groupSize[i];

    // Free allocated memory
	free(groupSize); groupSize = NULL;
	free(isAssigned); isAssigned = NULL;
    free(permutedIndexList); permutedIndexList = NULL; 
    free(permutedGroupList); permutedGroupList = NULL;
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
    int oldGroup, swap;


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
                if (Avg[j][i] > maxAvgCon && Rd[j]!=-1) {
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

// Function to process a partition 
void process_partition(double* groupDiversity, int* partition,  int* ub, int* childSolution,
     int* vectorElement, int K, int N, int *element_count, int *target_group) {
    int i, selectedGroup, processedCount, selectedElement;

    int elementCount = *element_count; 
    int targetGroup = *target_group;
    int maxGroupDiversity = -1;
    for (i = 0; i < K; i++) {
        if (groupDiversity_p1[i] > maxGroupDiversity) {
            maxGroupDiversity = groupDiversity_p1[i];
            selectedGroup = i;
        }   
    }

    for (i = 0; i < N; i++) {
        if (p1[i] == selectedGroup) {
            SelectEle[elementCount++] = i;
        }
    }

    int groupCount = 0;
    for (i = 0; i < K; i++) {
        if (ub[i] != -1 && ub[i] >= elementCount) {
            SelectGroup[groupCount++] = i;
        }
    }

    if (groupCount == 0) { // No valid group found
        int minDiff = 999999;
        for (i = 0; i < K; i++) {
            if (ub[i] != -1 && elementCount - ub[i] < minDiff) {
                minDiff = elementCount - ub[i];
                targetGroup = i;
            }
        }

        processedCount = 0;
        while (processedCount < elementCount - minDiff) {
            selectedElement = random_int(elementCount);
            do {
                selectedElement = (selectedElement + 1) % elementCount;
            } while (SelectEle[selectedElement] == -1);

            childSolution[SelectEle[selectedElement]] = targetGroup;
            SelectEleTemp[processedCount++] = SelectEle[selectedElement];
            vectorElement[SelectEle[selectedElement]] = -1;
            SelectEle[selectedElement] = -1;
        }
        elementCount = processedCount;
    } else {
        targetGroup = SelectGroup[random_int(groupCount)];
        for (i = 0; i < elementCount; i++) {
            childSolution[SelectEle[i]] = targetGroup;
            vectorElement[SelectEle[i]] = -1;
            SelectEleTemp[i] = SelectEle[i];
        }
    }
    *element_count = elementCount;
    *target_group = targetGroup;
}

void Crossover(int partition1[], int partition2[], int childSolution[], int scSizeGroup[]) {
    /* Algorithm 5: combines partitions in a way that maintains group constraints */

    int i, j;
    int elementCount, processedCount, selectedElement;
    int totalLowerBound, totalBelowLowerBound;

    // Initialize arrays
    for (i = 0; i < N; i++) {
        vectorElement[i] = i;
        childSolution[i] = -1;
    }
    for (i = 0; i < K; i++) {
        LBGroup[i] = 0;
        UBGroup[i] = 0;
        BigThanLB[i] = 0;
        groupElement[i] = i;
        ub[i] = UB[i];
        scSizeGroup[i] = 0;
    }

    // Initialize p1 with partition1
    for (i = 0; i < N; i++) {
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

    // Initialize p2 with partition2
    for (i = 0; i < N; i++) {
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

    int targetGroup = -1;
    // Main crossover process
    for (i = 0; i < K; i++) {
        elementCount = 0;
        if (uniform_rnd_number() < 0.5) {
            // Process partition 1
            process_partition(groupDiversity_p1, p1, ub, childSolution, vectorElement, K, N, &elementCount, &targetGroup);  
        } else {
            // Process partition2 (similar to partition1 logic)
            process_partition(groupDiversity_p2, p2, ub, childSolution, vectorElement, K, N, &elementCount, &targetGroup);  
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

    // Assign unprocessed elements to meet lower bounds Pseudo code line 22
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
            if (childSolution[j] == targetGroup) {
                SelectEle[elementCount++] = j;
            }
        }

        selectedElement = random_int(elementCount);
        childSolution[SelectEle[selectedElement]] = -1;
        vectorElement[SelectEle[selectedElement]] = SelectEle[selectedElement];
        scSizeGroup[targetGroup]--;
        if (scSizeGroup[targetGroup] == LB[targetGroup]) {
            BigThanLB[targetGroup] = 0;
        }
        processedCount++;
    }

    // Assign elements to meet lower bounds Pseudo code line 22
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
        childSolution[SelectEle[selectedElement]] = targetGroup;
        vectorElement[SelectEle[selectedElement]] = -1;
        scSizeGroup[targetGroup]++;
        if (scSizeGroup[targetGroup] == LB[targetGroup]) {
            LBGroup[targetGroup] = 0;
        }
        totalBelowLowerBound++;
    }

    // Assign elements to meet upper bounds Pseudo code line 23
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
        childSolution[SelectEle[selectedElement]] = targetGroup;
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

    Rprintf("Start lcoal serach criteriom");
    
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
		for (j = 0; j < K; j++) {
			Neighbors[count].type = 1;
			Neighbors[count].v = i;
			Neighbors[count].g = j;
			count++;
		}

    // Type 2 neighbors: (i, j) where each pair of elements (i, j) are neighbors    
	for (i = 0; i < N; i++)
		for (j = i + 1; j < N; j++) {
			Neighbors[count].type = 2;
			Neighbors[count].x = i;
			Neighbors[count].y = j;
			count++;
		}
		
}

void ClearDeltaMatrix() {
	/* Resets the delta_f matrix */
    for (int i = 0; i < N; ++i) {
        // should this not be over N ?!
        for (int j = 0; j < K; j++) {
            Delta_Matrix[i][j] = 0.0;
        }
    }
}

void ClearDeltaMatrixDispersion() {
	/* Resets the delta_f matrix */
    for (int i = 0; i < N; ++i) {
        // should this not be over N ?!
        for (int j = 0; j < K; j++) {
            Delta_Matrix[i][j] = INFINITY;
        }
    }
}

void BuildDeltaMatrix() {
	/*  Builds the delta_f matrix and calculates the objective function value */

	ClearDeltaMatrix();

    int i, j;
	// Update Delta_Matrix based on distances
    for (int i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            Delta_Matrix[i][s[j]] += Distances[i][j];
        }
    }

    // Calculate the objective function value
    objective = 0.0;
    for (i = 0; i < N; i++) {
        objective += Delta_Matrix[i][s[i]];
    }
    objective /= 2.0;
}

void BuildDeltaMatrixDispersion() {
	/*  Builds the delta_f matrix and calculates the objective function value */

	ClearDeltaMatrixDispersion();

    int i, j;
	// Update Delta_Matrix based on distances
    for (int i = 0; i < N-1; i++) {
        for (j = i+1; j < N; j++) {
            Delta_Matrix[i][s[j]] = fmin(Distances[i][j], Delta_Matrix[i][s[j]]);
        }
    } 

    // Calculate the objective function value
    objective = INFINITY;
    for (i = 0; i < N; i++) {
        objective = fmin(Delta_Matrix[i][s[i]], objective);
    }
}

void BuildGroupDiversityForCrossover() {
	/*  Builds group diversity values for crossover */
	
    // Initialize group diversity values to zero
	for (int i = 0; i < K; i++) groupDiversity[i] = 0.0;
	
	// Compute group diversity based on distances
    for (int i = 0; i < N; i++) {
        int group_i = s[i];
        for (int j = 0; j < N; j++) {
            if (group_i == s[j]) {
                groupDiversity[group_i] += Distances[i][j];
            }
        }
    }
}

void OneMoveUpdateDeltaMatrix(int i, int oldGroup, int newGroup) {
    // Update Delta_Matrix for all elements affected by the move
	for (int j = 0; j < N; j++) {
		if (j != i) {
			Delta_Matrix[j][oldGroup] -= Distances[i][j];
			Delta_Matrix[j][newGroup] += Distances[i][j];
		}
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
    for (i = 0; i < beta_max; i++) {
        S[i].s = (int*)malloc(N * sizeof(int));
        O[i].s = (int*)malloc(N * sizeof(int));
        S[i].SizeG = (int*)malloc(K * sizeof(int));
        O[i].SizeG = (int*)malloc(K * sizeof(int));
    }    
    
    Delta_Matrix = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta_Matrix[i] = (double*)malloc(K * sizeof(double));
    Delta_Matrix_p1 = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta_Matrix_p1[i] = (double*)malloc(K * sizeof(double));
    Delta_Matrix_p2 = (double**)malloc(N * sizeof(double*));
    for (i = 0; i < N; i++) Delta_Matrix_p2[i] = (double*)malloc(K * sizeof(double));
    groupDiversity = (double*)malloc(K * sizeof(double));
    groupDiversity_p1 = (double*)malloc(K * sizeof(double));
    groupDiversity_p2 = (double*)malloc(K * sizeof(double));
    
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

    int i;
    for (i = 0; i < beta_max; i++) {
        free(S[i].s); S[i].s = NULL;
        free(S[i].SizeG); S[i].SizeG = NULL;
        free(O[i].s); O[i].s = NULL;
        free(O[i].SizeG); O[i].SizeG = NULL;
    }
    free(S); S = NULL;
    free(O); O = NULL;
    
    free(LB); LB = NULL;
    free(UB); UB = NULL;
    free(Neighbors); Neighbors = NULL;
    
    for (i = 0; i < N; i++) {
        free(Delta_Matrix[i]); Delta_Matrix[i] = NULL;
        free(Delta_Matrix_p1[i]); Delta_Matrix_p1[i] = NULL;
        free(Delta_Matrix_p2[i]); Delta_Matrix_p2[i] = NULL;
    }
    free(Delta_Matrix); Delta_Matrix = NULL;
    free(Delta_Matrix_p1); Delta_Matrix_p1 = NULL;
    free(Delta_Matrix_p2); Delta_Matrix_p2 = NULL;
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

// Generate a random integer from zero to max-1
int random_int(int max) {
  GetRNGstate();
  double my_number = unif_rand();
  PutRNGstate();
  return (int) floor(my_number * max);
}
