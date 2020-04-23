
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include "declarations.h"

/* Exchange Method for Anticlustering
 * param *data: vector of data points (in R, this is a data frame,
 *         the matrix structure must be restored in C)
 * param *N: The number of elements (i.e., number of "rows" in *data)
 * param *M: The number of features (i.e., number of "cols" in *data)
 * param *K: The number of clusters
 * param *frequencies: The number of elements per cluster, i.e., an array
 *         of length *K.
 * param *clusters: An initial assignment of elements to clusters,
 *         array of length *N (has to consists of integers between 0 and (K-1) 
 *         - this has to be guaranteed by the caller)
 * int *USE_CATS A boolean value (i.e., 1/0) indicating whether categorical 
 *         constraints
 * param *C: The number of categories
 * param *CAT_frequencies: The number of elements per category, i.e., an array
 *         of length *C.
 * param *categories: An assignment of elements to categories,
 *         array of length *N (has to consists of integers between 0 and (C-1) 
 *         - this has to be guaranteed by the caller)
 * 
 * The return value is assigned to the argument `clusters`, via pointer
 * 
 * 
 * ============================ Some explanations ============================
 * 
 * Throughout this method, some data structures are defined that are used
 * throughout. The function works by first initializing the required data 
 * structures on which an exchange optimization algorithm is concuted.
 * 
 * These are the data structures that are "global" throughout the function:
 * 
 * 1. `struct element POINTS[n]` - An array copy of the `n` data points in the 
 *    same order as the input. Stores the `m` values per data point and the 
 *    cluster assignment of the data point (initially passed via `*clusters`)
 * 2. `struct node *HEADS[k]` - An array of length `k` where each
 *    entry points to the head of a list that represents a cluster. Each cluster
 *    is implemented as a linked list where each node points to an `element` in 
 *    `struct element POINTS[n]`. Thus, there are `k` cluster lists that are 
 *    used during the algorithm to compute the variance by cluster. During
 *    the exchange method, elements are swapped between cluster lists 
 *    (implemented through the function `swap()`). 
 *    The function `initialize_cluster_heads()` sets up the pointer array; the 
 *    function `fill_cluster_lists()` fills the data points as nodes into the 
 *    lists.
 * 3. `struct node *PTR_NODES[n]` - Array of pointers to nodes, used for 
 *    iterating during the exchange method. Points to elements in the cluster 
 *    lists but can be used to iterate through the data in the original order
 *    of the data (across cluster lists, this order is lost).
 * 4. `double CENTERS[k][m]` - A matrix of cluster centers.
 * 5. `double OBJECTIVE_BY_CLUSTER[k]` - The variance objective by cluster.
 * 6. `double SUM_VAR_OBJECTIVE` - The total variance objective (i.e., sum 
 *    of all entries in `OBJECTIVE_BY_CLUSTER`).
 *    
 * Regarding the inclusion of categorical constraints: If the argument 
 * `USE_CATS` is `TRUE` (i.e, `1`), the exchange method restricts the 
 * exchange partners to elements of the same category. What elements are part 
 * of the same category is specified via the input argument `categories`. 
 * Generally, the arguments `C`, `CAT_frequencies`, `categories` have the same 
 * semantic as `K`, `frequencies` and `clusters`, respectively. However, they 
 * represent a fixed category and not a variable cluster affiliation (that is 
 * changed in this function to maximize the k-means criterion). If `USE_CATS` is
 * `FALSE` (i.e, `0`), the arguments `C`, `CAT_frequencies`, `categories` are 
 * not used.
 * 
 * To balance a categorical variable that is represented by `categories` across
 * clusters, it is necessary that these are already balanced when calling this
 * function. `c_anticlustering` will not obtain an initial balanced partitioning
 * itself - the caller is responsible.
 * 
 * ===========================================================================
*/

void c_anticlustering(double *data, int *N, int *M, int *K, int *frequencies,
        int *clusters, int *USE_CATS, int *C, int *CAT_frequencies,
        int *categories) {
        
        const size_t n = (size_t) *N; // number of data points
        const size_t m = (size_t) *M; // number of variables per data point
        const size_t k = (size_t) *K; // number of clusters
        
        // Set up array of data points, fill it, return if memory runs out
        struct element POINTS[n];
        if (fill_data_points(data, n, m, POINTS, clusters, USE_CATS, categories) == 1) {
                return;
        }
        
        // Set up array of exchange partners
        if (*USE_CATS) {
                const size_t c = (size_t) *C; // number of categories 
                size_t *C_HEADS[c]; // describe above in general explanation 
                // writes the indices into C_HEADS:
                if (category_indices(n, c, POINTS, C_HEADS, 
                                     categories, CAT_frequencies) == 1) {
                        free_points(n, POINTS, n);
                        return;
                }
        }
        
        // Set up array of pointer-to-cluster-heads, return if memory runs out
        struct node *HEADS[k];
        if (initialize_cluster_heads(k, HEADS) == 1) {
                free_points(n, POINTS, n);
                return; 
        }

        // Set up array of pointers-to-nodes, return if memory runs out
        struct node *PTR_NODES[n];
        if (fill_cluster_lists(n, k, clusters, POINTS, PTR_NODES, HEADS) == 1) {
                return;
        }
        
        // Set up matrix of cluster centers
        double CENTERS[k][m];
        for (size_t i = 0; i < k; i++) {
                compute_center(m, CENTERS[i], HEADS[i], frequencies[i]);
        }
        
        // Get variance objective of the initial cluster assignment
        double OBJ_BY_CLUSTER[k]; 
        objective_by_cluster(m, k, OBJ_BY_CLUSTER, CENTERS, HEADS);
        double SUM_VAR_OBJECTIVE = array_sum(k, OBJ_BY_CLUSTER);
        
        /* Some variables for bookkeeping during the optimization */
        
        size_t best_partner;
        double tmp_centers[k][m];
        double best_centers[k][m];
        double tmp_objs[k];
        double best_objs[k];
        double tmp_obj;
        double best_obj;
        
        /* Start main iteration loop for exchange procedure */
        
        /* 1. Level: Iterate through `n` data points */
        for (size_t i = 0; i < n; i++) {
                size_t cl1 = PTR_NODES[i]->data->cluster;
                
                // Initialize `best` variable for the i'th item
                best_obj = 0;
                copy_matrix(k, m, CENTERS, best_centers);
                copy_array(k, OBJ_BY_CLUSTER, best_objs);
                
                /* 2. Level: Iterate through `n` exchange partners */
                for (size_t j = 0; j < n; j++) {
                        
                        size_t cl2 = PTR_NODES[j]->data->cluster;
                        // no swapping attempt if in the same cluster:
                        if (cl1 == cl2) { 
                                continue;
                        }

                        // Initialize `tmp` variables for the exchange partner:
                        copy_matrix(k, m, CENTERS, tmp_centers);
                        copy_array(k, OBJ_BY_CLUSTER, tmp_objs);
                        
                        update_centers(
                                k, m, 
                                tmp_centers, 
                                PTR_NODES[i],
                                PTR_NODES[j],
                                frequencies
                        );
                        swap(n, i, j, PTR_NODES);
                        // Update objective
                        tmp_objs[cl1] = cluster_var(m, HEADS[cl1], tmp_centers[cl1]);
                        tmp_objs[cl2] = cluster_var(m, HEADS[cl2], tmp_centers[cl2]);
                        tmp_obj = array_sum(k, tmp_objs);
                        
                        // Update `best` variables if objective was improved
                        if (tmp_obj > best_obj) {
                                best_obj = tmp_obj;
                                copy_matrix(k, m, tmp_centers, best_centers);
                                copy_array(k, tmp_objs, best_objs);
                                best_partner = j;
                        }
                        
                        // Swap back to test next exchange partner
                        swap(n, i, j, PTR_NODES);
                }
                
                // Only if objective is improved: Do the swap
                if (best_obj > SUM_VAR_OBJECTIVE) {
                        swap(n, i, best_partner, PTR_NODES);
                        // Update the "global" variables
                        SUM_VAR_OBJECTIVE = best_obj;
                        copy_matrix(k, m, best_centers, CENTERS);
                        copy_array(k, best_objs, OBJ_BY_CLUSTER);
                }
        }

        // Write output
        for (size_t i = 0; i < n; i++) {
                clusters[i] = PTR_NODES[i]->data->cluster;
        }
        
        // in the end, free allocated memory:
        free_points(n, POINTS, n);
        free_nodes(k, HEADS);
        // TODO: Free memory allocated for category indexes!
}

/* 
 * Perform a swap between two elements.
 * 
 * param `size_t n`: Number of data points
 * param `size_t i`: Index of first element to be swapped
 * param `size_t j`: Index of second element to be swapped
 * param `struct node *PTR_NODES[n]`: pointer to nodes
 * 
 */

void swap(size_t n, size_t i, size_t j, struct node *PTR_NODES[n]) {
        
        struct node *one = PTR_NODES[i];
        struct node *two = PTR_NODES[j];
        
        // Get cluster indices
        size_t cl1 = one->data->cluster;
        size_t cl2 = two->data->cluster;
        
        // Update pointer in `PTR_NODES`
        size_t ID1 = one->data->ID;
        size_t ID2 = two->data->ID;
        PTR_NODES[ID1] = two;
        PTR_NODES[ID2] = one;
        
        // Update ID of the elements
        one->data->ID = ID2;
        one->data->ID = ID1;
        
        // Update the cluster affiliation
        one->data->cluster = cl2;
        two->data->cluster = cl1;
        
        // Update nodes in cluster lists 
        struct element *tmp = one->data;
        one->data = two->data;
        two->data = tmp;
        
}

/* Update cluster centers after a swap between two nodes in two cluster lists */
void update_centers(size_t k, size_t m, double CENTERS[k][m],
                    struct node *one, struct node *two, int *frequencies) {
        size_t cl1 = one->data->cluster;
        size_t cl2 = two->data->cluster;
        for (int i = 0; i < m; i++) {
                double change_cl1 = one->data->values[i] / frequencies[cl1];
                double change_cl2 = two->data->values[i] / frequencies[cl2];
                // Update first cluster center
                CENTERS[cl1][i] = CENTERS[cl1][i] + change_cl2;
                CENTERS[cl1][i] = CENTERS[cl1][i] - change_cl1;
                // Update second cluster center
                CENTERS[cl2][i] = CENTERS[cl2][i] - change_cl2;
                CENTERS[cl2][i] = CENTERS[cl2][i] + change_cl1;
        }
}

/* Compute variance for a cluster
 * param `size_t m`: Number of variables per data point
 * param `struct node *HEAD`: Pointer to a cluster list HEAD
 * param `double center[m]`: Array of mean feature values in the cluster
 */
double cluster_var(size_t m, struct node *HEAD, double center[m]) {
        double sum = 0;
        struct node *tmp = HEAD->next;
        while (tmp != NULL) {
                sum += euclidean_squared(center, tmp->data->values, m);
                tmp = tmp->next;
        }
        return sum;
}

/* Compute cluster center for one cluster
 * 
 * param `size_t m`: Number of variables per data point
 * param `double center[m]`: Empty array of cluster centers
 * param `struct node *HEAD`: Pointer to a cluster list HEAD
 * param `int freq`: Number of elements in the cluster
 * 
 * The input array `center` is changed through this function to represent
 * one cluster center.
 * 
 */
void compute_center(size_t m, double center[m], struct node *HEAD, int freq) {
        // Initialize center matrix as 0:
        for (size_t i = 0; i < m; i++) {
                center[i] = 0; 
        }
        
        struct node *tmp = HEAD->next; 
        while (tmp != NULL) {
                for (size_t i = 0; i < m; i++) {
                        center[i] = center[i] + tmp->data->values[i];
                }
                tmp = tmp->next;
        } 
        // To get cluster centers: Divide by number of elements 
        for (size_t i = 0; i < m; i++) {
                center[i] = center[i] / freq;
        }
}

/* Extracted method that fill data points into array of struct `element`
 * param `double *data` pointer to original data array of length n
 * param `size_t m`: Number of data points
 * param `size_t m`: Number of variables per data point
 * param `struct element points[n]`: Array to be filled with data points
 * param `int *clusters`: Cluster affiliation of the n data points
 * 
 * return: `0` if all data points were successfully stored; `1` if not.
 * 
 */
int fill_data_points(double *data, size_t n, size_t m, struct element POINTS[n], 
                     int *clusters, int *USE_CATS, int *categories) {
        // Create offset variable to correctly read out data points
        int m_ptr[m];
        for (size_t i = 0; i < m; i++) {
                m_ptr[i] = i * n;
        }
        
        // Size of a data vector per element:
        size_t data_size = m * sizeof(POINTS[0].values[0]);
        
        for (size_t i = 0; i < n; i++) {
                POINTS[i].cluster = clusters[i];
                if (*USE_CATS) {
                        POINTS[i].category = categories[i];
                } else {
                        POINTS[i].category = 0;
                }
                POINTS[i].ID = i;
                POINTS[i].values = (double*) malloc(data_size);
                if (POINTS[i].values == NULL) {
                        free_points(n, POINTS, i);
                        print_memory_error();
                        return 1;
                } 
                // Fill data into `element`:
                for (size_t j = 0; j < m; j++) {
                        POINTS[i].values[j] = data[m_ptr[j]++];
                }
        }
        return 0;
}

/* Squared Euclidean Distance between two arrays of same length
 * 
 * param *x: Array / pointer to first element
 * param *y: Array / pointer to second element
 * param m: length of the two arrays
 * 
 * return: The squared euclidean distance
 * 
 */

double euclidean_squared(double *x, double *y, size_t m) {
        double sum = 0;
        for (size_t i = 0; i < m; i++) {
                sum = sum + pow(x[i] - y[i], 2);
        }
        return sum;
}

/* Copy one array into another */
void copy_array(size_t n, double origin[n], double target[n]) {
        for (int i = 0; i < n; i++) {
                target[i] = origin[i];
        }
}

/* Copy one matrix into another */
void copy_matrix(size_t n, size_t m, double origin[n][m], double target[n][m]) {
        for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                        target[i][j] = origin[i][j];
                }
        }
}

/* After creation, initialize the HEAD of each cluster list
 * 
 * We access the cluster list using an array of length `k` pointing
 * to the `k` HEADs of the clusters. The first element pointed to 
 * in the list (i.e., the HEAD) is "empty". It just points to the 
 * first real element that is in that cluster.
 * 
 * param `size_t k`: The number of clusters
 * param `struct node *PTR_CLUSTER_HEADS[k]`: The array of pointers to 
 *     cluster HEADS.
 * 
 *  * return: `0` if the cluster list could be initialized successfully, `1` 
 *      if not (in that case, there was no memory that could be allocated).
 * 
 */

int initialize_cluster_heads(size_t k, struct node *HEADS[k]) {
        for (size_t i = 0; i < k; i++) {
                HEADS[i] = (struct node*) malloc(sizeof(struct node*));
                if (HEADS[i] == NULL) {
                        free_nodes(k, HEADS);
                        print_memory_error();
                        return 1;
                }
                HEADS[i]->next = NULL;
                HEADS[i]->data = NULL;
        }
        return 0;
}

/* After initialization, fill the cluster lists with data
 * This function does two things at the same time:
 * (a) add each data point as a node to a cluster list,
 * (b) store the pointer to each node in the array `PTR_NODES`
 */
int fill_cluster_lists(size_t n, size_t k, int *clusters,
                       struct element POINTS[n], struct node *PTR_NODES[n],
                       struct node *PTR_CLUSTER_HEADS[k]) {
        for (size_t i = 0; i < n; i++) {
                struct node *current_cluster = PTR_CLUSTER_HEADS[clusters[i]];
                PTR_NODES[i] = append_to_cluster(current_cluster, &POINTS[i]);
                if (PTR_NODES[i] == NULL) { // failed to allocate memory
                        free_points(n, POINTS, n);
                        free_nodes(k, PTR_CLUSTER_HEADS);
                        print_memory_error();
                        return 1; 
                }
        }
        return 0;
}

/* Append data point to linked list
 * 
 * param `struct node *HEAD` Pointer to HEAD of cluster list
 * param `struct element *data`: Pointer to the data point 
 *     that is appended to the list
 * 
 * return: Pointer to the `node` of the element that was appended to the
 *     cluster list
 * 
 */

struct node* append_to_cluster(struct node *HEAD, struct element *data) {
        struct node *tmp = HEAD->next; // may be NULL if list is empty
        HEAD->next = (struct node*) malloc(sizeof(struct node*));
        if (HEAD->next == NULL) {
                return NULL; // failed to allocate memory
        }
        // New element is right next to HEAD element:
        HEAD->next->data = data;
        HEAD->next->next = tmp; // tmp was next to HEAD before
        return HEAD->next; 
}

/* Compute sum of squared errors between cluster center and data points (i.e.
 * variance) for each of the k clusters. The input array `VAR_OBJECTIVE` is 
 * changed through this function. 
 */
void objective_by_cluster(size_t m, size_t k, double OBJ_BY_CLUSTER[k], 
                          double CENTERS[k][m], struct node *HEADS[k]) {
        for (size_t i = 0; i < k; i++) {
                OBJ_BY_CLUSTER[i] = cluster_var(m, HEADS[i], CENTERS[i]);
        }
}

/* Compute the sum of an array */
double array_sum(size_t k, double ARRAY[k]) {
        double sum = 0;
        for (size_t i = 0; i < k; i++) {
                sum += ARRAY[i];
        }
        return sum;
}


/* Get a structure that has multiple arrays, pointed to by entries in another array; 
 * each of the arrays represents a category
 */
int category_indices(size_t n, size_t c, struct element POINTS[n], 
                     size_t *C_HEADS[c], int *categories, 
                     int *CAT_frequencies) {
        
        struct node *HEADS[c]; // used for filling `C_HEADS` - convert list into array
        if (initialize_cluster_heads(c, HEADS) == 1) {
                return 1; 
        }
        
        // Set up array of pointers-to-nodes, return if memory runs out
        struct node *PTR_NODES[n];
        if (fill_cluster_lists(n, c, categories, POINTS, PTR_NODES, HEADS) == 1) {
                return 1;
        }
        
        // Initialize the index arrays
        size_t n_cats;
        struct node *tmp; 
        size_t j;
        for (size_t i = 0; i < c; i++) {
                n_cats = (size_t) CAT_frequencies[i];
                C_HEADS[i] = (size_t*) malloc(n_cats * sizeof(size_t));
                if (C_HEADS[i] == NULL) {
                        print_memory_error();
                        return 1;
                }
                // Now write `C_HEADS`! Fills all `c` arrays with indices, 
                // based on the category lists `HEAD[0], HEAD[1], ..., HEAD[c]`
                tmp = HEADS[i]->next;
                j = 0;
                while (tmp != NULL) {
                        C_HEADS[i][j] = tmp->data->ID;
                        printf("%zu ", C_HEADS[i][j]);
                        tmp = tmp->next;
                } 
                printf("\nmembers in category %zu: %zu\n", i, n_cats);
        }

        // free temporary category lists
        free_nodes(c, HEADS);
        return 0;
}

/* Free memory in the cluster lists
* param `size_t k`: The number of clusters
* param `struct node *PTR_CLUSTER_HEADS[k]`: The array of pointers to 
*     cluster HEADS
*/
void free_nodes(size_t k, struct node *PTR_CLUSTER_HEADS[k]) {
        struct node *ptr;
        struct node *prev; // using temp pointer for freeing
        for (size_t i = 0; i < k; i++) {
                ptr = PTR_CLUSTER_HEADS[i];
                while (ptr->next != NULL)
                {  
                        prev = ptr;
                        ptr = ptr->next;
                        free(prev);
                }
                free(ptr);
        }
        
}

/* Free memory in the data points
 * param `size_t n`: length of array `POINTS`
 * param `struct element POINTS[n]`: Array containing data points
 * param `size_t i`: The index up to which data points are freed
 */
void free_points(size_t n, struct element POINTS[n], size_t i) {
        for (size_t j = 0; j < i; j++) {
                free(POINTS[j].values);
        }
}

void print_memory_error() {
        fprintf(stderr, "Failed to allocate enough memory.");
}

