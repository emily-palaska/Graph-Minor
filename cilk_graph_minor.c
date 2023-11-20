#include <stdio.h>
#include <stdlib.h>
#include <cilk/cilk.h>
#include "mmio.c"
#include "mmio.h"
#include "ctimer.h"

int *I, *J, *im, *jm, *im_c, *jm_c;
double *V, *vm, *vm_c;
int *vec;

int m, n;
int c, p;
int count;

//parallel function for vector calculation
void par_vecex(int index) {
	while (index < m) {
		im[index] = vec[I[index] - 1];
		jm[index] = vec[J[index] - 1];
		vm[index] = V[index];
		
		index += p;
	}
}

//parallel function for vector compression
void par_comp(int k) {
	int pivi, pivj;		//here we store the coordenate set based on the index
	int foundi = 0;		//first occurence index of coordinate set
	double sum = 0;		//sum of all found values
	
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
				//keep first occurence
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

// sequential function to calculate (repeating) vectors of Graph Minor
int minor_vec() {
	int i, nt;
	
	//memory allocation	
	im = malloc(sizeof(int) * m);
	jm = malloc(sizeof(int) * m);
	vm = malloc(sizeof(double) * m);
	
	//error handling	
	if (!im || !jm || !vm) {
		printf("Memory allocation failed.\n");
		return 1;
	}
	
	//decide number of processors	
	if (m < p) nt = m;
	else nt = p;
	
	//extract in parallel
	cilk_scope {
		for(i = 0; i < nt; i++) cilk_spawn par_vecex(i);
	}
	return 0;
}

// sequential function to compress the calculated vectors of Graph Minor
int vec_comp() {
	int i, nt;	
	
	//decide number of spawns
	if (m < p) nt = m;
	else nt = p;
	
	//compress in parallel
	cilk_scope {
		for(i = 0; i < nt; i++) cilk_spawn par_comp(i);	
	}
	
	//measure length of result
	count = 0;
	for(i = 0; i < m; i++) if(im[i]) count++;
	if(count * 3 > c * c) printf("Warning! Not memory-friendly to store in vector format.\n");
	
	//create temporary vectors
	im_c = malloc(sizeof(int) * count);
	jm_c = malloc(sizeof(int) * count);
	vm_c = malloc(sizeof(double) * count);
	
	//error handling
	if(!im_c || !jm_c || !vm_c) {
		printf("Failed memory allocation.\n");
		return 1;	
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
	
	return 0;
}


//--------------------------------------------
//--------------INPUT VECTOR------------------
//--------------------------------------------


//initiliazes the configuration vector and find number of clusters
int init_vec() {
    	c = n / 10;
    	printf("Merging into %d clusters.\n", c);
    	vec = malloc(sizeof(int) * n);
    	if(!vec) {
    		printf("Failed memory allocation.\n");
    		return 1;
    	}
    	for(int b = 0; b < n; b ++){
        	vec[b] = b % c + 1;
    	}
    	return 0;
}

//imports matrix from matrix market format
int import_matrix(char *filename) {
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

//----------MAIN-----------
int main(int argc, char *agrv[]) {
	//file name and number of threads is given as input in terminal
	//eg ./gminor "file.mtx" "16"
	char* filename = agrv[1];
	p = atoi(agrv[2]);
	if(argc != 3) {
		printf("False arguments passed to main.\n.");
		return 1;
	}
	int i;
	
	
	//Importing Matrix
	if (import_matrix(filename)) {
		printf("Importing failed.\n");
		return 1;
	}
	
	//Initializing configuration vector
	if(init_vec()) {
		printf("Could't initialize vector.\n");
		return 1;
	}
	
	//Start timer
	ctimer_t t;
	ctimer_start(&t);
	
	//Extracting result vectors
	if(minor_vec()) {
		printf("Extraction failed.\n");
		return 1;
	}
	
	//Compressing result vectors
	if(vec_comp()) {
		printf("Compression failed.\n");
		return 1;
	}
	
	//Stop timer
	ctimer_stop(&t);
	ctimer_measure(&t);
	
	
	//print results (un-comment for vectors)
	printf("Number of edges in graph minor: %d\n", count);
	//printf("~~~RESULTS~~~\n");
	//printf("i\tj\tv\n");	
	//for(i = 0; i < count; i++)
	//	printf("%d\t%d\t%lf\n", im[i], jm[i], vm[i]);
	
	ctimer_print(t, "gminor");
	
	free(vec);
	free(im);
	free(jm);
	free(vm);	
	free(I);
	free(J);
	free(V);
	return 0;
	
}

