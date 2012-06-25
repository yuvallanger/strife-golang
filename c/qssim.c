/////////////////////////////////////////////////
// quorum sensing simulation
// written by Elad Shtilerman and Yael Baran 
// (shtiler@gmail.com,yael.baran@mail.huji.ac.il)
/////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/////////////
// uncomment this flag if parameters are provided in a file
#define _CONF
/////////////

// constants
#define GRID_SIZE 100
#define COMP_PER_IT GRID_SIZE*GRID_SIZE
#define IT_NUM 10000
#define BASIC_COST 100 	// cost carried by all

#ifndef _CONF
#define D 1				// diffusion
#define SIG_RANGE 1 	// how far is singl absorbed
#define PG_RANGE 1		// how far is pg absorbed
#define PROD_THR 6		// how many signaling neighbors it takes for pg production
#define PG_THR 6 		// how many producing neighbors it takes to enjoy pg (n_q)
#define PROD_COST 30 	// additive to basic cost
#define PG_RELIEF 0.1 	// 0<=r<=1 and r<BASIC_COST/(BASIC_COST+PROD_COST)
						// so that cooperator with pg is better than defector without
#define P_MUT 5e-4		// mutation probability
#else
int D, SIG_RANGE, PG_RANGE, PROD_THR, PG_THR, PROD_COST;
double PG_RELIEF, P_MUT;
#endif




// general constants
#define STRL 200

// initialization modes
#define ALL_R0S0 1 // entire grid is WT strain

// strain definitions
int r4strain[4] = {0,0,1,1};
int s4strain[4] = {0,1,0,1};
int strain_spec[2][2] = {{0,1},{2,3}}; // index using [r][s]


// update signal, prod and pg according to specified change in grid
void update_arrays(char** grid, int*** signal, char** prod, int** pg,
					int oldstrain, int x, int y) {
		
	int sillyi,sillyj,sillyk,sillyl,i,j,k,l;
	int oldprod;
	
	for (sillyi=x-SIG_RANGE;sillyi<=x+SIG_RANGE;sillyi++) {
		for (sillyj=y-SIG_RANGE;sillyj<=y+SIG_RANGE;sillyj++) {
			i = (sillyi+GRID_SIZE)%GRID_SIZE;
			j = (sillyj+GRID_SIZE)%GRID_SIZE;
			// update signal level in sig range
			signal[s4strain[oldstrain]][i][j]--;
			signal[s4strain[(int)grid[x][y]]][i][j]++;
			// update producer status in sig range
			oldprod = prod[i][j];
			prod[i][j] = (PROD_THR<=signal[r4strain[(int)grid[i][j]]][i][j]) ? 1 : 0;
			// update pg level in sig+pg range
			if (oldprod != prod[i][j]) {
				for (sillyk=i-PG_RANGE;sillyk<=i+PG_RANGE;sillyk++) {
					for (sillyl=j-PG_RANGE;sillyl<=j+PG_RANGE;sillyl++) {
						k = (sillyk+GRID_SIZE)%GRID_SIZE;
						l = (sillyl+GRID_SIZE)%GRID_SIZE;
						pg[k][l] += (prod[i][j]-oldprod);
					}
				}	
			}
		}
	}
		
}		
							
// given a strain, mutates each of the loci with prob P_MUT
char mutate (int strain) {

	int r,s;
	
	r = r4strain[strain];
	if (drand48()<=P_MUT) {
		r = 1-r;
	}
	
	s = s4strain[strain];
	if (drand48()<=P_MUT) {
		s = 1-s;
	}
	
	return strain_spec[r][s];
	
}

// computes the fitness of a given cell
double comp_fit(char** grid, char** prod, int** pg, int x, int y) {
	
	double cost = BASIC_COST;
	if (prod[x][y]) {
		cost += PROD_COST;
	}
	if (PG_THR<=pg[x][y]) {
		cost *= PG_RELIEF;
	}
	return BASIC_COST/cost;

}

// updates the grid by iterating COMP_PER_IT times the following :
// - choose a random cell
// - try to mutate it, and let the mutant and the original cell compete
// - choose a random neighbor
// - let the neighbors compete
void update_grid(char** grid, int*** signal, char** prod, int** pg, FILE* fp) {

	int x1, y1, dx, dy, x2, y2, i, j;
	char oldstrain, newstrain;
	double f1, f2, curfit, newfit;
	
	static int counts[4] = {0,0,0,0};
	static int first = 1;
	if (first) {
		first = 0;
		for (i=0;i<GRID_SIZE;i++) {
			for (j=0;j<GRID_SIZE;j++) {
				counts[(int)grid[i][j]]++;
			}
		}
	}

	for (i=0;i<COMP_PER_IT;i++) {
	
		x1 = drand48()*GRID_SIZE;
		y1 = drand48()*GRID_SIZE;

		// try mutation current spot
		newstrain = mutate(grid[x1][y1]);
		if (grid[x1][y1] != newstrain) {
			curfit = comp_fit(grid,prod,pg,x1,y1);
			oldstrain = grid[x1][y1];
			grid[x1][y1] = newstrain;
			update_arrays(grid,signal,prod,pg,oldstrain,x1,y1);
			newfit = comp_fit(grid,prod,pg,x1,y1);
			if (drand48() < curfit/(curfit+newfit)) {
				grid[x1][y1] = oldstrain;
				update_arrays(grid,signal,prod,pg,newstrain,x1,y1);
			} else {
				counts[(int)oldstrain]--;
				counts[(int)newstrain]++;
			}
		}

		// look invading a direct neighbor
		do {
			dx = drand48()*3-1;;
			dy = drand48()*3-1;;
		} while (dx==0 && dy==0);
		x2 = (x1+dx+GRID_SIZE)%GRID_SIZE;
		y2 = (y1+dy+GRID_SIZE)%GRID_SIZE;
		if (grid[x1][y1]==grid[x2][y2]) {
			continue;
		}
		f1 = comp_fit(grid,prod,pg,x1,y1);
		f2 = comp_fit(grid,prod,pg,x2,y2);
		if (drand48() < f1/(f1+f2)) {
			oldstrain = grid[x2][y2];
			grid[x2][y2] = grid[x1][y1];
			update_arrays(grid,signal,prod,pg,oldstrain,x2,y2);
		} else {
			oldstrain = grid[x1][y1];
			grid[x1][y1] = grid[x2][y2];
			update_arrays(grid,signal,prod,pg,oldstrain,x1,y1);
		}
		counts[(int)oldstrain]--;
		counts[(int)grid[x1][y1]]++;
				
	}
	
	for (i=0;i<4;i++) {
		fprintf(fp,"%10d",counts[i]);
	}
	fprintf(fp,"\n");

}

// initialize grid according to specified mode
void init_grid(char** grid, int mode) {

	int i,j;
	
	if (mode==ALL_R0S0) { // initialize all to strain 0
		for (i=0;i<GRID_SIZE;i++) {
			for (j=0;j<GRID_SIZE;j++) {
				grid[i][j] = 0;
			}
		}
	}

}


// initialize signal, prod, pg according to grid
void init_arrays(char** grid, int*** signal, char** prod, int** pg) {


	int i,j,k,l,sillyk,sillyl;
	
	// compute signal level
	for (i=0;i<GRID_SIZE;i++) {
		for (j=0;j<GRID_SIZE;j++) {
			for (sillyk=i-SIG_RANGE;sillyk<=i+SIG_RANGE;sillyk++) {
				for (sillyl=j-SIG_RANGE;sillyl<=j+SIG_RANGE;sillyl++) {
					k = (sillyk+GRID_SIZE)%GRID_SIZE;
					l = (sillyl+GRID_SIZE)%GRID_SIZE;
					signal[s4strain[(int)grid[i][j]]][k][l]++;			
				}
			}
		}
	}

	// compute prod status
	for (i=0;i<GRID_SIZE;i++) {
		for (j=0;j<GRID_SIZE;j++) {
			prod[i][j] = (PROD_THR<=signal[r4strain[(int)grid[i][j]]][i][j]) ? 1 : 0;
		}
	}
	
	// compute pg level
	for (i=0;i<GRID_SIZE;i++) {
		for (j=0;j<GRID_SIZE;j++) {
			if (!prod[i][j]) {
				continue;
			}	
			for (sillyk=i-PG_RANGE;sillyk<=i+PG_RANGE;sillyk++) {
				for (sillyl=j-PG_RANGE;sillyl<=j+PG_RANGE;sillyl++) {
					k = (sillyk+GRID_SIZE)%GRID_SIZE;
					l = (sillyl+GRID_SIZE)%GRID_SIZE;
					pg[k][l]++;			
				}
			}		
		}
	}


}

// perform diffusion steps in grid
void diffuse (char** grid, int*** signal, char** prod, int** pg) {

	char keep[2][2] = {{0,0},{0,0}};
	char new[2][2] = {{0,0},{0,0}};
	int howmany = D*COMP_PER_IT/4; // one cell is changed per udpate iteration
	int it,x,y,clockwise;
	
	for (it=0;it<howmany;it++) {
	
		x = drand48()*GRID_SIZE;
		y = drand48()*GRID_SIZE;
		keep[0][0] = grid[x][y];
		keep[0][1] = grid[x][(y+1)%GRID_SIZE];
		keep[1][0] = grid[(x+1)%GRID_SIZE][y];
		keep[1][1] = grid[(x+1)%GRID_SIZE][(y+1)%GRID_SIZE];
		
		clockwise = drand48()+0.5;
		if (clockwise) {
			new[0][0] = keep[1][0];
			new[1][0] = keep[1][1];		
			new[1][1] = keep[0][1];
			new[0][1] = keep[0][0];		
		} else {
			new[0][0] = keep[0][1];
			new[0][1] = keep[1][1];
			new[1][1] = keep[1][0];
			new[1][0] = keep[0][0];
		}
		
		grid[x][y] = new[0][0];
		grid[x][(y+1)%GRID_SIZE] = new[0][1];
		grid[(x+1)%GRID_SIZE][y] = new[1][0];
		grid[(x+1)%GRID_SIZE][(y+1)%GRID_SIZE] = new[1][1];
		
		update_arrays(grid,signal,prod,pg,keep[0][0],x,y);
		update_arrays(grid,signal,prod,pg,keep[0][1],x,(y+1)%GRID_SIZE);
		update_arrays(grid,signal,prod,pg,keep[1][0],(x+1)%GRID_SIZE,y);
		update_arrays(grid,signal,prod,pg,keep[1][1],(x+1)%GRID_SIZE,(y+1)%GRID_SIZE);
	
	}
		
}


// output grid to file while counting strains
void print_grid(FILE* fp, char** grid) {

	int i,j;
	int counts[4] = {0,0,0,0};
	for (i=0;i<GRID_SIZE;i++) {
		for (j=0;j<GRID_SIZE;j++) {
			counts[(int)grid[i][j]]++;
			fprintf(fp,"%2d",grid[i][j]);
		}
		fprintf(fp,"\n");
	}

}


int main (int argc, char **argv) {


	#ifdef _CONF
	FILE* fpconf = fopen (argv[1], "r");
	if (fpconf == NULL) {
    	printf ("cannot open file %s\n", argv[1]);
    	exit(-1);
  	}	
	char s[STRL];
	char* tok;
	char outfile[STRL];
	int printGrid;
	while (fgets(s,STRL,fpconf) != NULL) {
    	tok = strtok (s, " \t");
    	if (!strcmp(tok,"d")) {
      		tok = strtok (NULL, " \t");
      		tok = strtok (NULL, " \t\n");
      		D = atoi(tok);
    	} else if (!strcmp(tok,"sig_range")) {
      		tok = strtok (NULL, " \t");
      		tok = strtok (NULL, " \t\n");
      		SIG_RANGE = atoi(tok);
    	} else if (!strcmp(tok,"pg_range")) {
      		tok = strtok (NULL, " \t");
      		tok = strtok (NULL, " \t\n");
      		PG_RANGE = atoi(tok);
     	} else if (!strcmp(tok,"prod_thr")) {
      		tok = strtok (NULL, " \t");
      		tok = strtok (NULL, " \t\n");
      		PROD_THR = atoi(tok);
     	} else if (!strcmp(tok,"pg_thr")) {
      		tok = strtok (NULL, " \t");
      		tok = strtok (NULL, " \t\n");
      		PG_THR = atoi(tok);
    	} else if (!strcmp(tok,"prod_cost")) {
      		tok = strtok (NULL, " \t");
      		tok = strtok (NULL, " \t\n");
      		PROD_COST = atoi(tok);
    	} else if (!strcmp(tok,"pg_relief")) {
      		tok = strtok (NULL, " \t");
      		tok = strtok (NULL, " \t\n");
      		PG_RELIEF = atof(tok);
    	} else if (!strcmp(tok,"p_mut")) {
      		tok = strtok (NULL, " \t");
      		tok = strtok (NULL, " \t\n");
      		P_MUT = atof(tok);
    	} else if (!strcmp(tok,"outfile")) {
      		tok = strtok (NULL, " \t");
      		tok = strtok (NULL, " \t\n");
      		strcpy(outfile,tok);
    	} else if (!strcmp(tok,"printgrid")) {
      		tok = strtok (NULL, " \t");
      		tok = strtok (NULL, " \t\n");
      		printGrid = atoi(tok);
    	} else {
    		printf("unknown parameter in confile - %s; quitting\n",tok);
    		exit(1);
    	}
    }
    fclose(fpconf);
 	#endif
	
	// init rand
	unsigned short tmps[3];
	unsigned long t = time(NULL);
	tmps[0] = 0; tmps[1] = (t >> 16); tmps[2] = t & 0xFFFF;
	seed48(tmps);
	
	// allocate grid, prod, pg
	// grid - 0/1/2/3 - contains strain number
	// prod - 0/1 - is the cell producing pg ?
	// pg - how much pg is measured in this position
	int i,j;
	char** grid = (char**) malloc (sizeof(char*)*GRID_SIZE);
	char** prod = (char**) malloc (sizeof(char*)*GRID_SIZE);
	int** pg = (int**) malloc (sizeof(int*)*GRID_SIZE);
	for (i=0;i<GRID_SIZE;i++) {
		grid[i] = (char*) malloc (GRID_SIZE);
		prod[i] = (char*) malloc (GRID_SIZE);
		pg[i] = (int*) malloc (sizeof(int)*GRID_SIZE);
		for (j=0;j<GRID_SIZE;j++) {
			grid[i][j] = -1;
			prod[i][j] = 0;
			pg[i][j] = 0;		
		}
	}

	// allocate signal
	// how much signal of this type is measured in this position
	// index using [signal(0/1)][row][col]
	int k;
	int*** signal = (int***) malloc (sizeof(int**)*2);
	for (i=0;i<2;i++) {
		signal[i] = (int**) malloc (sizeof(int*)*GRID_SIZE);
		for (j=0;j<GRID_SIZE;j++) {
			signal[i][j] = (int*) malloc (sizeof(int)*GRID_SIZE);
			for (k=0;k<GRID_SIZE;k++) {
				signal[i][j][k] = 0;	
			}
		}
	}
	
	printf("initializing arrays\n");
	fflush(stdout);
	init_grid(grid,ALL_R0S0);
	init_arrays(grid,signal,prod,pg);

	char filename[STRL];
	strcpy(filename,outfile);
	strcat(filename,".count");
	FILE* fpc = fopen(filename,"w");
	FILE* fpg;
	
	printf("iterating ...\n");
	fflush(stdout);
	for (i=0;i<IT_NUM;i++) {
		diffuse(grid,signal,prod,pg);
		update_grid(grid,signal,prod,pg,fpc);
		if (printGrid) {
			if (!(i%printGrid)) {
				sprintf(filename,"%s.grid.%d",outfile,i);
				fpg = fopen(filename,"w");
				print_grid(fpg,grid);
				fclose(fpg);
			}
		}
	}
	fclose(fpc);
	printf("completed %d iterations\n",i);
	
	printf("freeing memory and quitting\n");
	
	// freeing memory

	for (i=0;i<GRID_SIZE;i++) {
		free(grid[i]);
		free(prod[i]);
		free(pg[i]);
	}	
	free(grid);
	free(prod);
	free(pg);
	
	for (i=0;i<2;i++) {
		for (j=0;j<GRID_SIZE;j++) {
			free(signal[i][j]);
		}
		free(signal[i]);
	}
	free(signal);
	
	return (0);
	
}
