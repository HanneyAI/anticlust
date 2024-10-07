#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <assert.h>
#include "three-phase-header.h"

extern int beta_max;
extern int N, K;  // node number and group number
extern double** Distances;   // distance matrix between elements
extern double** DistancesT;
extern int* LB; // Lower bound of the number of elements in the group i 
extern int* UB; // Upper bound of the number of elements in the group i
extern double theta, theta_max, theta_min; 
extern double alpha;
extern int beta_min; 
extern int eta_max;

// main processure
extern int maxNumberIterations;
extern double start_time, Time_limit;
extern Solution S_b; //best solution
extern Solution CS;
extern Solution *S; //S_i
extern Solution *O; //O_i in crossover
extern Neighborhood *Neighbors;

//double neighboorhood local search
extern int *s; // partition array for each v
extern double objective;
double *min_distance_per_cluster;
int **min_distance_tuple;

// Crosssover
double *min_distance_per_cluster_p1;
int **min_distance_tuple_p1;
double *min_distance_per_cluster_p2;
int **min_distance_tuple_p2;
extern int *SelectEle;
extern int *SelectEleTemp;
extern int *SelectGroup;
extern int *p1;
extern int *p2;

// for ;crossover
extern int *vectorElement;
extern int *groupElement;
extern int *LBGroup;
extern int *UBGroup;
extern int *BigThanLB;
extern int *ub;

//directed pertubation
extern int* Rd, * UnderLB; //Rd=R
extern int *SizeG; //c_g

void adding(int new_ind, int g, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster);
void removing(int removed_ind, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster);
void swapping(int ind1, int ind2, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster);
double evaluate_objective(double *s_min_distance_per_cluster);
void fill_arrays(int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster);
void initialize_arrays(int **s_min_distance_tuple, double *s_min_distance_per_cluster);
void DoubleNeighborhoodLocalSearchDispersion(int partition[], int SizeGroup[], double* cost);
void SearchAlgorithmDisperion(void);
void InitialSolDispersion(Solution *S);
void UndirectedPerturbationDispersion(int L, int partition[], int SizeGroup[]);
void DoubleNeighborhoodLocalSearchDispersion(int partition[], int SizeGroup[], double* cost);
void CrossoverDispersion(int partition1[], int partition2[], int score[], int scSizeGroup[]);
void DirectPerturbationDispersion(int eta_max, int partition[], int SizeGroup[]);
void AssignMemoryDispersion();
void ReleaseMemoryDispersion();

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
void three_phase_search_dispersion(
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


  Rprintf("C is in dispersion\n");

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
 // Rprintf("Distance matrix D =\n");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      Distances[i][j] = distances[i * N + j];
 //     Rprintf("%f ",Distances[i][j]);
      DistancesT[i][j] = 2 * distances[i * N + j];
    }
  //  printf("\n");
  }
  
  min_distance_per_cluster = (double *)malloc(K * sizeof(double));
  min_distance_tuple = malloc(K * sizeof(int *));
  for (int i = 0; i < K; i++) {
    min_distance_tuple[i] = malloc(2 * sizeof(int));
  }
  min_distance_per_cluster_p1 = (double *)malloc(K * sizeof(double));
  min_distance_tuple_p1 = malloc(K * sizeof(int *));
  for (int i = 0; i < K; i++) {
    min_distance_tuple_p1[i] = malloc(2 * sizeof(int)); // IMPORTANT: typo: missing p1 / p2 below
  }
  min_distance_per_cluster_p2 = (double *)malloc(K * sizeof(double));
  min_distance_tuple_p2 = malloc(K * sizeof(int *));
  for (int i = 0; i < K; i++) {
    min_distance_tuple_p2[i] = malloc(2 * sizeof(int));
    
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

  AssignMemoryDispersion();
  if (*mem_error == 1) {
    return;
  }
  

  BuildNeighbors();
  
  SearchAlgorithmDisperion();
  
  //save S_b -> solution with result
  for (int i = 0; i < N; i++){
    result[i] = S_b.s[i];
    Rprintf("result [%d] %d", i, S_b.s[i]);
  }
  *score = S_b.cost;
  Rprintf("cost %f", S_b.cost);
  
  // Remember to free the allocated memory after use
  for (int i = 0; i < K; i++){
        free(min_distance_tuple[i]); min_distance_tuple[i] = NULL;
        free(min_distance_tuple_p1[i]); min_distance_tuple_p1[i] = NULL;
        free(min_distance_tuple_p2[i]); min_distance_tuple_p2[i] = NULL;
  }
  free(min_distance_per_cluster); min_distance_per_cluster = NULL;
  free(min_distance_tuple); min_distance_tuple = NULL;
  free(min_distance_per_cluster_p1); min_distance_per_cluster_p1 = NULL;
  free(min_distance_tuple_p1); min_distance_tuple_p1 = NULL;
  free(min_distance_per_cluster_p2); min_distance_per_cluster_p2 = NULL;
  free(min_distance_tuple_p2); min_distance_tuple_p2 = NULL;
  for (int i = 0; i < N; i++) {
    free(Distances[i]); Distances[i] = NULL;
    free(DistancesT[i]); DistancesT[i] = NULL;
  }
  free(Distances); Distances = NULL;
  free(DistancesT); DistancesT = NULL;
  free(LB); LB = NULL;
  free(UB); UB = NULL;

  Rprintf("Until now runs trhough.");
  
  ReleaseMemoryDispersion();
}

void initialize_arrays(int **s_min_distance_tuple, double *s_min_distance_per_cluster) {
    for (int k = 0; k < K; k++) {
        s_min_distance_per_cluster[k] = INFINITY;
        s_min_distance_tuple[k][0] = 0;
        s_min_distance_tuple[k][1] = 0;
    }
}

void fill_arrays(int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster) {
    initialize_arrays(s_min_distance_tuple, s_min_distance_per_cluster);
    // Function to fill the arrays based on the distance matrix and cluster assignments
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            if (Distances[i][j] < s_min_distance_per_cluster[partition[i]] && partition[i] == partition[j]) { // IMPORTANT: add "&& partition[i] = ..." to check that both data points are in the same cluster!
                s_min_distance_per_cluster[partition[i]] = Distances[i][j];
                s_min_distance_tuple[partition[i]][0] = i;
                s_min_distance_tuple[partition[i]][1] = j;
            }
        }
    }
}


// Function to evaluate the objective function
double evaluate_objective(double *s_min_distance_per_cluster) {
    double f = s_min_distance_per_cluster[0];
    for (int k = 1; k < K; k++) {
        f = fmin(f, s_min_distance_per_cluster[k]);
    }
    return f;
}

// Function to add an element to a cluster
void adding(int new_ind, int cluster, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster) {
    for (int i = 0; i < N; i++) {
        if (partition[i] == cluster && i != new_ind) { // IMPORTANT: add "&& i != new_ind" to make sure to not include d[i][i] = 0!
            // printf("Calculating difference between %d and %d\n",i,new_ind);
            if (Distances[i][new_ind] < s_min_distance_per_cluster[cluster]) {
                s_min_distance_per_cluster[cluster] = Distances[i][new_ind];
                s_min_distance_tuple[cluster][0] = i;
                s_min_distance_tuple[cluster][1] = new_ind;
                // printf("New tuple (%d,%d) with score %f for cluster %d recorded\n",i,new_ind,s_min_distance_per_cluster[cluster],cluster);
            }
        }
    }
    partition[new_ind] = cluster;
}

// Function to remove an element from its cluster
void removing(int removed_ind, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster)  {
    int g = partition[removed_ind];
    partition[removed_ind] = -1;  // Temporarily hide the index
    if (s_min_distance_tuple[g][0] == removed_ind || s_min_distance_tuple[g][1] == removed_ind) {
		// only requires update of objective function if remove_ind appears in min_distance_tuple
        s_min_distance_per_cluster[g] = INFINITY;
        for (int i = 0; i < N - 1; i++) {
            if (partition[i] == g) {
                for (int j = i + 1; j < N; j++) {
                    if (partition[j] == g) {
                        if (Distances[i][j] < s_min_distance_per_cluster[g]) {
                            s_min_distance_per_cluster[g] = Distances[i][j];
                            s_min_distance_tuple[g][0] = i;
                            s_min_distance_tuple[g][1] = j;
                        }
                    }
                }
            }
        }
    }
}

void swapping(int ind1, int ind2, int *partition, int **s_min_distance_tuple, double *s_min_distance_per_cluster) {
    int g1 = partition[ind1];
    int g2 = partition[ind2];

    // Temporarily hide ind1
    partition[ind1] = -1;
    // IMPORTANT: add the next four lines
    // if either ind1/ind2 belong to their respective min_distance_tuple, ensure to recalculate the corresponding s_min_distance_per_cluster when adding ind1/ind2
    if (ind1 == s_min_distance_tuple[g1][0] || ind1 == s_min_distance_tuple[g1][1]) s_min_distance_per_cluster[g1] = INFINITY;
    if (ind2 == s_min_distance_tuple[g2][0] || ind2 == s_min_distance_tuple[g2][1]) s_min_distance_per_cluster[g2] = INFINITY;
    // the correct recalculation of s_min_distance_tuple happens in the adding-functions below

    // Add ind2 to the cluster of g1
    adding(ind2, g1, partition, s_min_distance_tuple, s_min_distance_per_cluster);

    // Add ind1 to the cluster of g2
    adding(ind1, g2, partition, s_min_distance_tuple, s_min_distance_per_cluster);
}


void SearchAlgorithmDisperion() {
    /* Algorithm 1: The main procedure of TPSDP. */

    int eta;
    int pickedSolution;

    //important! for windows and linux there is a differnt definition of this time
    //on windows its the wall time, on linux the CPU time
    clock_t start_time = clock();
    S_b.cost = -1;  //dispersion should be maximized
    int i;
   // printf("Initial best solution (should have random entries):\n");
    //for(i = 0; i < N; i++) {
      //  printf("Solution S_B[%d] : %d\n", i, S_b.s[i]);
    //}
     //printf("Solution S_B cost : %f\n", S_b.cost);
    
    // Initial population generation
    int j, k;
    for (i = 0; i < beta_max; i++) {
        InitialSolDispersion(&CS);
        for (j = 0; j < N; j++) S[i].s[j] = CS.s[j];
        for (k = 0; k < K; k++) S[i].SizeG[k] = CS.SizeG[k];
        S[i].cost = CS.cost;
        if (S[i].cost > S_b.cost) {
            for (j = 0; j < N; j++) S_b.s[j] = S[i].s[j];
            for (k = 0; k < K; k++) S_b.SizeG[k] = S[i].SizeG[k];
            S_b.cost = S[i].cost;
        }
    }
    //printf("Initial Sol completed, best solution is now:\n");
  //  for(i = 0; i < N; i++) {
      //  printf("Solution S_B[%d] : %d\n", i, S_b.s[i]);
    //}
    // printf("Solution S_B cost : %f\n", S_b.cost);
 
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
   //         printf("Starting UndirectedPerturbation...\n");
           UndirectedPerturbationDispersion(eta, S[i].s, S[i].SizeG);
     //      printf("Starting LocalSearch...\n");
           DoubleNeighborhoodLocalSearchDispersion(S[i].s, S[i].SizeG, &S[i].cost);
            if (S[i].cost > S_b.cost) {
                for (j = 0; j < N; j++) S_b.s[j] = S[i].s[j];
                for (k = 0; k < K; k++) S_b.SizeG[k] = S[i].SizeG[k];
                S_b.cost = S[i].cost;
            }
        }
       // printf("UndirectedPerturbation completed, best solution is now:\n");
        //for(i = 0; i < N; i++) {
          //  printf("Solution S_B[%d] : %d\n", i, S_b.s[i]);
        //}
        //printf("Solution S_B cost : %f\n", S_b.cost);
        

        // Crossover and Local Search
        if (beta_max > 1) {
            for (i = 0; i < beta_max; i++) {
                pickedSolution = random_int(beta_max);
                do {
                    pickedSolution = (pickedSolution + 1) % beta_max;
                } while (pickedSolution == i);
            CrossoverDispersion(S[i].s, S[pickedSolution].s, O[i].s, O[i].SizeG);
            DoubleNeighborhoodLocalSearchDispersion(O[i].s, O[i].SizeG, &O[i].cost);
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
        
       // printf("Crossover completed, best solution is now:\n");
        //for(i = 0; i < N; i++) {
          //  printf("Solution S_B[%d] : %d\n", i, S_b.s[i]);
        //}
        //printf("Solution S_B cost : %f\n", S_b.cost);

        // Direct Perturbation and Local Search
        for (i = 0; i < beta_max; i++) {
            DirectPerturbationDispersion(eta_max, S[i].s, S[i].SizeG);
            DoubleNeighborhoodLocalSearchDispersion(S[i].s, S[i].SizeG, &S[i].cost);
            if (S[i].cost > S_b.cost) {
                for (j = 0; j < N; j++) S_b.s[j] = S[i].s[j];
                for (k = 0; k < K; k++) S_b.SizeG[k] = S[i].SizeG[k];
                S_b.cost = S[i].cost;
            }
        }

      //  printf("DirectPerturbation completed, best solution is now:\n");
        //for(i = 0; i < beta_max; i++) {
          //  printf("Solution i[%d] : %f\n", i, S[i].cost);
        //}
        //printf("Best Solution : %f\n", S_b.cost);

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
void InitialSolDispersion(Solution *S) {
    /* Algorithm 2: initializes a solution S */
    RandomInitialSol(S->s, S->SizeG);
    DoubleNeighborhoodLocalSearchDispersion(S->s, S->SizeG, &(S->cost));
}

void DoubleNeighborhoodLocalSearchDispersion(int partition[], int SizeGroup[], double* cost) {

   // printf("Start local search");
    int i, v, g, u;
    int imp;

    // Initialize the partition array
    for (i = 0; i < N; i++) s[i] = partition[i];

    // Initialize the delta_f value and tuple arrays
    double delta_f = -99999.0;
    fill_arrays(s, min_distance_tuple, min_distance_per_cluster);

    int g1, g2;
    double old_f1, old_f2;
    do { 
        imp = 0;  // Reset improvement flag

        // First loop: Move individual elements to improve partition
        for (v = 0; v < N; v++) {
            for (g = 0; g < K; g++) {
                // Check if moving `v` to `group` is valid and beneficial
                if ((s[v] != g) && (SizeGroup[s[v]] > LB[s[v]]) && (SizeGroup[g] < UB[g])) {
                    g1 = s[v];
                    old_f1 = min_distance_per_cluster[s[v]];
                    old_f2 = min_distance_per_cluster[g];
                    removing(v, s, min_distance_tuple, min_distance_per_cluster);  //remove 
                    	// do not forget to carry out the adding to another cluster afterwards using the procedure above
	                adding(v, g, s, min_distance_tuple, min_distance_per_cluster);

                    delta_f = min_distance_per_cluster[g] - old_f2 + min_distance_per_cluster[g1] - old_f1;
                    if (delta_f <= 0) {
                        // revert changes
                        s[v] = g1;
                        min_distance_per_cluster[g1] = old_f1;
                        min_distance_per_cluster[g] = old_f2;
                    } else {
                        // printf("LocalSearch-Push: an improvement was detected!\n");
                        // objective = evaluate_objective(min_distance_per_cluster);
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
                    g1 = s[v];
                    g2 = s[u];
                    old_f1 = min_distance_per_cluster[s[v]];
                    old_f2 = min_distance_per_cluster[s[u]];
                    swapping(u, v, s, min_distance_tuple, min_distance_per_cluster); 

                    delta_f = min_distance_per_cluster[g2] - old_f2 + min_distance_per_cluster[g1] - old_f1; // IMPORTANT: typo: g2 instead of g
                    if (delta_f <= 0) {
                        // revert changes
                        s[v] = g1;
                        s[u] = g2;
                        min_distance_per_cluster[g1] = old_f1;
                        min_distance_per_cluster[g2] = old_f2;
                    } else {
                        // printf("LocalSearch-Swap: an improvement was detected!\n");
                        // objective = evaluate_objective(min_distance_per_cluster);
                        imp = 1;
                    }
                }
            }
        }
    } while (imp == 1);  // Continue until no improvement is made
    
    // IMPORTANT: add the next line. Before writing objective to *cost, first update this value!
    objective = evaluate_objective(min_distance_per_cluster);

    // Update the partition array with the final assignments
    for (i = 0; i < N; i++) partition[i] = s[i];
 //   printf("Solution : %f\n", objective);
    *cost = objective;
}

void UndirectedPerturbationDispersion(int L, int partition[], int SizeGroup[]) {
    /* Algorithm 4: Undirected Perturbation. Applies a strong perturbation to the partition */

   // printf("Start Search Algortihm");
    int current_index;
    int v, g, x, y;
    int oldGroup, swap;

 //   printf("Starting UndirectedPerturbation, previous partition is:\n");
    for (int i = 0; i < N; i++) {
        s[i] = partition[i];
 //       printf("s[%d] = %d\n",i,s[i]);
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
       //         printf("\t\t\tThis message should never appear!\n");
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
 //   printf("Finished UndirectedPerturbation, new partition is:\n");
    for (int i = 0; i < N; i++) {
        partition[i] = s[i];
//        printf("s[%d] = %d\n",i,s[i]);
    }
}

void DirectPerturbationDispersion(int eta_max, int partition[], int SizeGroup[]) {
    /* Algorithm 6: Directed Perturbation. 
	Iteratively refines partitions to balance group sizes and minimize costs */

//    printf("Starting Directed Pertubation\n");

    int i, j, k, L, number, minElement;
    int new_ind, ind1, ind2, selectedElement;
    double objective_k, minDeltaValue, maxDeltaValue;

    // Initialize the partition and size groups
    for (i = 0; i < N; i++) s[i] = partition[i];
    for (j = 0; j < K; j++) SizeG[j] = SizeGroup[j];
    fill_arrays(s, min_distance_tuple, min_distance_per_cluster);
    
    for (k = 0; k < K; k++) {
        printf("min_distance_per_cluster[%d] = %f\n",k,min_distance_per_cluster[k]);
    }

    // Main loop for perturbation iterations
    for (L = 0; L < eta_max; L++) {

        // Reset tracking variables
        number = 0;
        for (i = 0; i < K; i++) {
            UnderLB[i] = 0;
            Rd[i] = -1;
        }

        for (k = 0; k < K; k++) {
            // Pseudo group 8 to 13
            /* Remove the element x, y from the tuple of cluster k which have the lowest dipserion.
             Calculate the minimal dispersion value after removing x and y in cluster k. 
             The element which leads to the lower minimum dispersion of the cluster k after its removal
             will be removed from the cluster k. */
          
            // store values before removing() in temporary values. This speeds-up the adding-process. In short:
            ind1 = min_distance_tuple[k][0];
            ind2 = min_distance_tuple[k][1];
            objective_k = min_distance_per_cluster[k];
            removing(ind1, s, min_distance_tuple, min_distance_per_cluster);
            minDeltaValue = min_distance_per_cluster[k];
            // restore changes simply from ind1 and objective_k:
            min_distance_tuple[k][0] = ind1;
            min_distance_tuple[k][1] = ind2;
            min_distance_per_cluster[k] = objective_k;
            s[ind1] = k; // restore class membership
            removing(ind2, s, min_distance_tuple, min_distance_per_cluster);
             if (min_distance_per_cluster[k] < minDeltaValue) {
                minDeltaValue = min_distance_per_cluster[k];
                minElement = ind2;
            } else {
                min_distance_tuple[k][0] = ind1;
                min_distance_tuple[k][1] = ind2;
                min_distance_per_cluster[k] = objective_k;
                s[ind2] = k; // restore class membership
                removing(ind1, s, min_distance_tuple, min_distance_per_cluster);
                minElement = ind1;
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
        
    //    for (k = 0; k < K; k++) {
      //      printf("min_distance_per_cluster[%d] = %f\n",k,min_distance_per_cluster[k]);
        //}

        // line 15 of pseudo group will be removed, since it is not necessary for dispersion     
		
        // Handle groups that are under the lower bound (LB)
        int selectedGroup;
        int nn = 0;
        int k;
        while (nn < number) {
            // pseudo code line 17 - 29
            k = random_int(K);

            // Find the element with the highest average connection to the group
            do {
                k = (k + 1) % K;
            } while (UnderLB[k] == 0);
            
            // IMPORTANT: add the next five lines here
            // IF cluster k has less than 2 elements its min_distance is set to INFINITY.
            // In this case, we can safely set it to 0
            objective_k = min_distance_per_cluster[k];
            bool is_infinite = objective_k == INFINITY;
            if (is_infinite) objective_k = 0.0;

            // store values before adding() in temporary values. This speeds-up the removing-process. In short:
            maxDeltaValue = -INFINITY; // min_delta cannot become positive by adding a new data point to a cluster for dispersion
            for (i = 0; i < K; i++) {
                new_ind = Rd[i];
                if (new_ind > -1) {
         //           printf("Found a new_ind!\n");
                    ind1 = min_distance_tuple[k][0];
                    ind2 = min_distance_tuple[k][1];
                    adding(new_ind, k, s, min_distance_tuple, min_distance_per_cluster);
        //            printf("min_distance[%d] = %f (before it was %f)\n",k,min_distance_per_cluster[k],objective_k);
                    if (min_distance_per_cluster[k] - objective_k > maxDeltaValue) {
                        maxDeltaValue = min_distance_per_cluster[k] - objective_k; // IMPORTANT: typo min --> max!
                        selectedElement = new_ind;
                        selectedGroup = i;
                    }
                    // restore changes simply from ind1 and objective_k:
                    min_distance_tuple[k][0] = ind1;
                    min_distance_tuple[k][1] = ind2;
                    // IMPORTANT: replace "min_distance_per_cluster[k] = objective_k;" by the if-else statement, see other IMPORTANT above
                    if (is_infinite) min_distance_per_cluster[k] = INFINITY;
                    else min_distance_per_cluster[k] = objective_k;
                    s[new_ind] = -1; // remove new_ind from class
                }
            }
         //   printf("Selected Element %d pushed to cluster %d\n",selectedElement,k);
            adding(selectedElement, k, s, min_distance_tuple, min_distance_per_cluster);  
            
            UnderLB[k] = 0;
            Rd[selectedGroup] = -1;
            nn++;
        }

        // Handle remaining elements for groups that are above LB
        nn = 0;
        while (nn < K - number) {
            // pseuo group line 32 - 43
            selectedGroup = rand() % K;
            do {
                selectedGroup = (selectedGroup + 1) % K;
            } while (Rd[selectedGroup] == -1);

            maxDeltaValue = -INFINITY;
            new_ind = Rd[selectedGroup];
            Rd[selectedGroup] = -1;
            for (k = 0; k < K; i++) {
                // check that upper-bound ub[k] is not yet reached
                    if (SizeG[k] < UB[k]) {
                    ind1 = min_distance_tuple[k][0];
                    ind2 = min_distance_tuple[k][1];
                    objective_k = min_distance_per_cluster[k];
                    adding(new_ind, k, s, min_distance_tuple, min_distance_per_cluster);
                    if (min_distance_per_cluster[k] - objective_k > maxDeltaValue) {
                        minDeltaValue = min_distance_per_cluster[k] - objective_k;
                        selectedGroup = k;
                    }
                    // restore changes simply from ind1 and objective_k:
                    min_distance_tuple[k][0] = ind1;
                    min_distance_tuple[k][1] = ind2;
                    min_distance_per_cluster[k] = objective_k;
                    s[new_ind] = -1; // remove new_ind from class
                }
            }
            adding(new_ind, selectedGroup, s, min_distance_tuple, min_distance_per_cluster); 

            nn += 1;
        }
    }

    for (i = 0; i < N; i++) partition[i] = s[i];
    for (j = 0; j < K; j++) SizeGroup[j] = SizeG[j];
}

void CrossoverDispersion(int partition1[], int partition2[], int solutionChild[], int scSizeGroup[]) {
    /* Algorithm 5: combines partitions in a way that maintains group constraints */

 //   printf("Start Crossover");
    int i, j, maxGroupDispersion, selectedGroup;
    int elementCount, groupCount;
    int targetGroup = -1;
    int processedCount;
    int selectedElement;
    int totalLowerBound, totalBelowLowerBound;

    // Initialize p1 and p2s with partition1
    for (i = 0; i < N; i++) {
        p1[i] = partition1[i];
        p2[i] = partition2[i];
    }

    fill_arrays(p1, min_distance_tuple_p1, min_distance_per_cluster_p1);
    fill_arrays(p2, min_distance_tuple_p2, min_distance_per_cluster_p2);
    for (int k = 0; k < K; k++) {
        min_distance_per_cluster_p1[k] = min_distance_per_cluster[k];  //groupDiversity
        min_distance_per_cluster_p2[k] = min_distance_per_cluster[k];
        min_distance_tuple_p1[k][0] = min_distance_tuple[k][0];
        min_distance_tuple_p2[k][1] = min_distance_tuple[k][0];
    }
    
    // Initialize arrays
    for (i = 0; i < N; i++) {
        vectorElement[i] = i;
        solutionChild[i] = -1; // ?
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
        if (((double)rand() / RAND_MAX) < 0.5) {
            // Process partition1
            //find group with highest dispersion
            maxGroupDispersion = -1;
            for (j = 0; j < K; j++) {
                if (min_distance_per_cluster_p1[j] > maxGroupDispersion) {
                    maxGroupDispersion = min_distance_per_cluster_p1[j];
                    selectedGroup = j;
                }
            }

            // select elements to move
            elementCount = 0;
            for (j = 0; j < N; j++ ) {
                if (p1[j] == selectedGroup) {
                    SelectEle[elementCount++] = j;
                }
            }

           // choose groups who have enough space to hole seleccted groups
            groupCount = 0;
            for (j = 0; j < K; j++) {
                if (ub[j] != -1 && ub[j] >= elementCount) {
                    SelectGroup[groupCount++] = j;
                }
            }

            if (groupCount == 0) { 
                // 16  from Pseudogroup Algotihm 5 
                int minDiff = 999999;
                for (j = 0; j < K; j++) {
                    if (ub[j] != -1 && elementCount - ub[j] < minDiff) {
                        minDiff = elementCount - ub[j];
                        targetGroup = j;
                    }
                }

                processedCount = 0;
                while (processedCount < elementCount - minDiff) {
                    // 17  from Pseudogroup Algotihm 5 
                    selectedElement = random_int(elementCount);
                    do {
                        selectedElement = (selectedElement + 1) % elementCount;
                    } while (SelectEle[selectedElement] == -1);

                    solutionChild[SelectEle[selectedElement]] = targetGroup;
                    SelectEleTemp[processedCount++] = SelectEle[selectedElement];
                    vectorElement[SelectEle[selectedElement]] = -1;
                    SelectEle[selectedElement] = -1;
                }
                elementCount = processedCount;
            } else {
                targetGroup = SelectGroup[random_int(groupCount)];  // 13  from Pseudogroup Algotihm 5 
                for (j = 0; j < elementCount; j++) {
                     // 14  from Pseudogroup Algotihm 5 
                    solutionChild[SelectEle[j]] = targetGroup;
                    vectorElement[SelectEle[j]] = -1;
                    SelectEleTemp[j] = SelectEle[j];
                }
            }
        } else {
            // Process partition2 
            maxGroupDispersion = -1;
            for (j = 0; j < K; j++) {
                if (min_distance_per_cluster_p2[j] > maxGroupDispersion) {
                    maxGroupDispersion = min_distance_per_cluster_p2[j];
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

                    solutionChild[SelectEle[selectedElement]] = targetGroup;
                    SelectEleTemp[processedCount++] = SelectEle[selectedElement];
                    vectorElement[SelectEle[selectedElement]] = -1;
                    SelectEle[selectedElement] = -1;
                }
                elementCount = processedCount;
            } else {
                targetGroup = SelectGroup[random_int(groupCount)];
                for (j = 0; j < elementCount; j++) {
                    solutionChild[SelectEle[j]] = targetGroup;
                    vectorElement[SelectEle[j]] = -1;
                    SelectEleTemp[j] = SelectEle[j];
                }
            }
        }

        // Update group dispersion
        for (j = 0; j < elementCount; j++) {
            removing(SelectEleTemp[j], p1, min_distance_tuple_p1, min_distance_per_cluster_p1);
            removing(SelectEleTemp[j], p2, min_distance_tuple_p2, min_distance_per_cluster_p2);
            p1[SelectEleTemp[j]] = -1;
            p2[SelectEleTemp[j]] = -1;
        }

        // group should not be available anymore
        min_distance_per_cluster_p1[targetGroup] = -1;
        min_distance_per_cluster_p2[targetGroup] = -1;
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
            if (solutionChild[j] == targetGroup) {
                SelectEle[elementCount++] = j;
            }
        }

        selectedElement = random_int(elementCount);
        solutionChild[SelectEle[selectedElement]] = -1;
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
        solutionChild[SelectEle[selectedElement]] = targetGroup;
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
        solutionChild[SelectEle[selectedElement]] = targetGroup;
        vectorElement[SelectEle[selectedElement]] = -1;
        scSizeGroup[targetGroup]++;
        if (scSizeGroup[targetGroup] == UB[targetGroup]) {
            UBGroup[targetGroup] = 0;
        }
        totalSize++;
    }
}

void AssignMemoryDispersion() {
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
    
    CS.s = (int*)malloc(N * sizeof(int));
    S_b.s = (int*)malloc(N * sizeof(int));
    
    CS.SizeG = (int*)malloc(K * sizeof(int));
    S_b.SizeG = (int*)malloc(K * sizeof(int));
    
    Neighbors = (Neighborhood*)malloc((N * (N - 1) / 2 + N * K) * sizeof(Neighborhood));
    
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

void ReleaseMemoryDispersion() {
    /* responsible for reading the input file, 
    initializing matrices, and setting constraints on group sizes. */ 
    
    free(s); s = NULL;
    free(SizeG); SizeG = NULL;

    
    Rprintf("relesee dispersion Until now runs trhough.");
    // IMPORTANT: releasing S and O like for TPSDP leads to error!
    //int i;
   // for (i = 0; i < beta_max; i++) {
        //free(S[i].s); S[i].s = NULL;
        //free(S[i].SizeG); S[i].SizeG = NULL;
       // free(O[i].s); O[i].s = NULL;
        //free(O[i].SizeG); O[i].SizeG = NULL;
   // }

    free(CS.s); CS.s = NULL;
    free(CS.SizeG); CS.SizeG = NULL;
    free(S_b.s); S_b.s = NULL;
    free(S_b.SizeG); S_b.SizeG = NULL;
    free(O); O = NULL;
    free(S); S = NULL; 
    Rprintf("relesee dispersion Until now runs trhough 2.");

       
    free(LB); LB = NULL;
    free(UB); UB = NULL;
    free(Neighbors); Neighbors = NULL;

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
