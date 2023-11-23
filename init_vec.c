//alternative vector initialization taking a txt file as input
int init_vec(char *filecluster) {
	FILE *f;
	int max = 0;	
    	
    	//memory allocation and error handling
    	vec = malloc(sizeof(int) * n);
    	if(!vec) {
    		printf("Failed memory allocation.\n");
    		return 1;
    	}
    	
    	if ((f = fopen(filecluster, "r")) == NULL)
    		exit(1);
    	
    	for (int i = 0; i < n; i++) {
    		fscanf(f, "%d\n", &vec[i]);
    		if (vec[i] > max) max = vec[i];
    	}
    	
  	if(max) c = max;
  	else return 1;
    	printf("Merging into %d clusters.\n", c);
    	return 0;
}
