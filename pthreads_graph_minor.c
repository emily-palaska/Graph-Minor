#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include "mmio.c"
#include "mmio.h"
#define p 16        //number of threads we create at each parallel function
#define MODE 1      //mode 2 -> matrix already in sparse form (known I, J, V, M, N, c, v)
                    //mode 1 -> matrix market file
                    //mode 0 -> matrix to be converted to sparse format (known M, N, c, v)


// I,J,V: vectors for sparse format of A
// im, jm, vm: here we restore the result ins SPARSE FORMAT
// im_c, jm_c, vm_c; temporary pointers for calculations
// n/nd: size of square matrix A       m/md: number of non-zero entries in A
// c: number of clusters               vec: vector for configuration ids (size M/Md)


// define for MODE 0 or 2
#define nd 6
#define md 18

//define for MODE 1
#define filename "1n100.mtx"

//variables for all modes
int *I, *J, *im, *jm, *im_c, *jm_c;
double *V, *vm, *vm_c;
int *vec;
int m, n, c;
int count;
pthread_mutex_t mutex;


// just for MODE 0
double matrix[nd][nd] = {{1, 9, 0, 0, 0},
                         {0, 0, 5, 0, 0},
                         {0, 2, 0, 0, 0}};

// just for MODE 2
int II[md] = {1, 2, 6, 1, 2, 3, 5, 2, 3, 4, 6, 3, 4, 2, 5, 1, 3, 6};
int JJ[md] = {1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6};
double VV[md] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};


//////////////////////
//PARALLEL FUNCTIONS//
//////////////////////

void* par_matvec(void* arg) {
/* parallel call from matrix_to_vec() -> see there
*each thread has a unique index that correspons to a set number of elements in the matrix
*eg if matrix has 12 X 12 = 144 elements, each of the 16 threads handles 9 of them
*for each element, the sparse form is extracted if it's a non-zero value
*this means that with p processors and n X n elements the complexity is n ^ 2 / p
*/
    int index = *(int*)arg;
    while(index < n*n) {
        if (matrix[index / n][index % n] != 0) {
            pthread_mutex_lock(&mutex);
            I[count] = index / n + 1;
            J[count] = index % n + 1;
            V[count] = matrix[index / n][index % n];
            count++;
            pthread_mutex_unlock(&mutex);
        } 
        index += p;
    }
}

void* par_vecex(void* arg) {
/* parallel call from minor_vec() -> see there
*each thread has a unique index that correspons to a set number of edges in the original
*graph, eg if matrix has m = 80 edges, each of the 16 threads handles 5 of them
*for each edge, we caalculate the corresponding edge in the graph minor using the positions
*of the configuration ids vector vec.
*this means that with p processors and m edges the complexity is m / p
*/
    int index = *(int*)arg;
    while (index < m){
        im[index] = vec[I[index] - 1];
        jm[index] = vec[J[index] - 1];
        vm[index] = V[index];

        index += p;
    }
}

void* par_comp(void* arg) {
/* parallel call from vec_comp() -> see there
*each thread has a unique index that correspons to a set number of edges in the graph minor,
*eg if graph minor has c * c = 32 edges potentially, each of the 16 threads handles 2 of them.
*for each edge, we iterate the extracted graph minor vectors to find matching values.
*if none is found, the edge is not added to the result vectors but if one or more are located,
*the sum of their values are stored to the result vectors.
*the counter variable keeps track of the result vectors' length.
*this means that with p processors and m edges the complexity is c * c * m / p
*/
    int k = *(int*)arg;
    int pivi, pivj;
    int foundi = 0;
    double sum = 0;

    while(k < c*c) {
        foundi = 0;
        sum = 0;
        pivi = k / c + 1;
        pivj = k % c + 1;
        for(int i = 0; i < m; i++) {
            if(im[i] == pivi && jm[i] == pivj) {
                sum += vm[i];
                im[i] = 0;
                jm[i] = 0;
                vm[i] = 0;
                //keep first occurence and number of non-zero values
                if (!foundi) foundi = i + 1;
            } 
        }
        if(sum) {
            im[foundi - 1] = pivi;
            jm[foundi - 1] = pivj;
            vm[foundi - 1] = sum;
        }
        k += p;
    }
}


////////////////////////
//SEQUENTIAL FUNCTIONS//
////////////////////////

// function to initialize vec and find number of clusters
int init_vec() {
    c = n / 10;
    vec = malloc(sizeof(int) * n);
    if(!vec) {
        printf("Failed memory allocation.\n");
        return 1;
    }
    for(int b = 0; b < n; b++) {
        vec[b] = b % c + 1;
    }
    return 0;
}

// function from matrix market website for importing data
int import_matrix() {
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int rows, cols;   
    int i;

    if ((f = fopen(filename, "r")) == NULL) 
            exit(1);

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix and do necessary checks */

    if ((ret_code = mm_read_mtx_crd_size(f, &rows, &cols, &m)) !=0) exit(1);
    
    if(rows != cols) {
        printf("Matrix is not square, cannot run algorithm.\n");
        return 1;
    }
    else n = rows;

    if(m > n * n / 5) printf("Warning! Matrix is not sparse!\n");
    
    /* reseve memory for matrices */
    I = (int *) malloc(m * sizeof(int));
    J = (int *) malloc(m * sizeof(int));
    V = (double *) malloc(m * sizeof(double));
    if(!I || !J || !V) {
        printf("Failed memory allocattion.\n");
        exit(1);
    }
 
    // I, J are 1-based, needed for later in the program
    for (i = 0; i < m; i++)
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &V[i]);
    if (f !=stdin) fclose(f);

    return 0;
}

// function to convert a matrix from standard to sparse vector format
// used only in MODE 0
int matrix_to_vec() {  
    int i, nt;      //iteration index and number of threads
    int *thread_args;
    pthread_t *th;
    count = 0;
    n = nd;
    m = md;

    //error handling
    if(!th || !thread_args) {
        printf("Failed memory allocation.\n");
        exit(1);
    }

    //memory allocating for vectors
    I = malloc(sizeof(int) * m);
    J = malloc(sizeof(int) * m);
    V = malloc(sizeof(double) * m);

    //error handling
    if(!I || !J || !V) {
        printf("Failed memory allcoation.]\n");
        exit(1);
    }

    //deciding number of threads
    if(n * n <= p) nt = n * n;
    else nt = p;
    thread_args = malloc(sizeof(int) * nt);
    th = malloc(8 * nt);

    //creting threads
    for (i = 0; i < nt; i++) {
            thread_args[i] = i;
            if (pthread_create(th + i, NULL, &par_matvec, (void*)&thread_args[i]) !=0) {
                perror("Failed to create thread\n");
                return 1;
            }
    }
    //joining threaads
    for(i = 0; i < nt; i++) {
            if (pthread_join(th[i], NULL) != 0) return 2;
    }

    //memory de-allocation
    free(th);
    free(thread_args);
    return 0;
}

// function to find the UNCOMPRESSED vectors of the graph minor
int minor_vec() {
    int i, nt;      //iteration index and number of threads
    int *thread_args;
    pthread_t *th;

    //deciding number of threads
    if(m < p) nt = m;
    else nt = p;

    //memory allocation
    thread_args = malloc(sizeof(int) * nt);
    th = malloc(8 * nt);
    im = malloc(sizeof(int) * m);
    jm = malloc(sizeof(int) * m);
    vm = malloc(sizeof(double) * m);

    //error handliing
    if(!th || ! thread_args || !im || !jm || !vm) {
        printf("Failed memory allocation.\n");
        exit(1);
    }

    //create threads
    for(i = 0; i < nt; i++) {
        thread_args[i] = i;
        if (pthread_create(th + i, NULL, &par_vecex, (void*)&thread_args[i]) != 0) {
            perror("Failed to create thread");
            return 1;
        }
    }
    //join threads
    for(i = 0; i < nt; i++) {
        if(pthread_join(th[i], NULL) != 0) return 2;
    }

    //memory de-allocation
    free(th);
    free(thread_args);
    return 0;
}

// function to compress vectors of the graph minor
int vec_comp() {
    count = 0;
    int i, nt;      //iteration index and number of threads
    int *thread_args;
    pthread_t *th;

    //decide number of threads
    if(m < p) nt = m;
    else nt = p;

    //memory allocation
    thread_args = malloc(sizeof(int) * nt);
    th = malloc(8 * nt);

    //error handling
    if (!th || !thread_args) {
        printf("Could not allocate memory.\n");
        exit(1);
    }
    //thread creation
    for(i = 0; i < nt; i++) {
        thread_args[i] = i;
        if (pthread_create(th + i, NULL, &par_comp, (void*)&thread_args[i]) != 0) {
            perror("Failed to create thread");
            return 1;
        }
    }
    //join threads
    for(i = 0; i < nt; i++) {
        if(pthread_join(th[i], NULL) != 0) return 2;
    }
    
    //find length of result
    count = 0;
    for(i = 0; i < n; i++) {
        if (im[i]) count++;
    }
 
    // ifs for memory management
    if (count * 3 > c * c) printf("Warning! Not memory-friendly to store in vector format.\n");

    //create temporary vectors
    im_c = malloc(sizeof(int) * count);
    jm_c = malloc(sizeof(int) * count);
    vm_c = malloc(sizeof(double) * count);

    //error handling
    if(!im_c || !jm_c || !vm_c) {
        printf("Failed memory allocation.\n");
        exit(1);
    }

    //copy compressed form to temporary vectors
    int c = 0;
    for(int i = 0; i < m; i++) {
        if(im[i]) {
            im_c[c] = im[i];
            jm_c[c] = jm[i];
            vm_c[c] = vm[i];
            c++;
        }
    }

    //clear uncompressed vectors and store new values
    free(im);
    free(jm);
    free(vm);
    im = im_c;
    jm = jm_c;
    vm = vm_c;      

    //memory de-allocation
    free(th);
    free(thread_args);
    pthread_mutex_destroy(&mutex);
    return 0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char *argv[]) {
    int i;                              //iteration index
    pthread_mutex_init(&mutex, NULL);   //initialize mutex
    
    //INPUT HANDLING
    if(MODE == 0) {
        if(matrix_to_vec()) {
            printf("Conversion failed.\n");
            return 1;
        }
    }
    else if(MODE == 1) {
        if(import_matrix()) {
            printf("Importing failed.\n");
            return 1;
        }
    }
    else if(MODE == 2) {
        n = nd;
        m = md;
        I = II;
        J = JJ;
        V = VV;        
    }
    else{
        printf("Specify correct input MODE at the start of code.\n");
        return 1;
    }

    // INITIALIZATION OF CONFIGURATION VECTOR 
    if(init_vec()) printf("Couldn't find number of clusters.\n");

    // EXTRACTION OF RESULT VECTORS 
    if(minor_vec()) {
        printf("Extraction failed.\n");
        return 1;        
    }
    
    // COMPRESSION RESULT VECTORS
    if(vec_comp()) {
        printf("Vector compression failed.\n");
        return 1;
    }

    // RESULT PRINTING
    printf("\n~~~~~RESULTS~~~~~\n");
    printf("Number of non-zero elements in graph minor: %d\n", count);
    printf("i\tj\tv\n");
    for(i = 0; i < count; i++) printf("%d\t%d\t%lf\n", im[i], jm[i], vm[i]);

    // MEMORY DE-ALLOCATION AND TERMINATION
    pthread_mutex_destroy(&mutex);
    free(vec);
    free(im);
    free(jm);
    free(vm);
    if(MODE != 2) {
        free(I);
        free(J);
        free(V);
    }
    return 0;
}
