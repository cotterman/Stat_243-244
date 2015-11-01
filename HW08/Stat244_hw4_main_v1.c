#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tree1.h"
#include "btree1.h"
#define LOC(a,b)  ((a) * p + (b))
#include "Stat244_hw4_header1.h" //this basically means to copy and paste code from the .h file into the top of this code document
#include "timer.h"
FILE *infile;


/******************************************************************************/
/******************************************************************************/

struct OBS {
	double value;
	int count;
	int level;
};

int comp(struct OBS *mydat1, struct OBS *mydat2){
	if ((*mydat1).value<(*mydat2).value) return 1;
	else if ((*mydat1).value>(*mydat2).value) return -1;
	else return 0;
}
void printtraverse(struct OBS *myobs){
	printf("Level %d: value = %lf, count = %d \n", (*myobs).level, (*myobs).value, (*myobs).count);
}

int lvlcount = 0; //beware -- this is a global variable
void leveltraverse(struct OBS *myobs){
	(*myobs).level = lvlcount;
	lvlcount++;
}

int main(int argc, char *argv[]){

	//exercise control over which parts of program are run;
	  int run_Q1 = 1;
	  int run_Q2 = 1;

	int i, n, k, t, j, r, c, tnobs, anobs, nvars;
	double myd;
	tnobs = 1000; //number of elements in vector for tree exercise
	double matT[tnobs]; //vector for tree exercise, unsorted
	double matTS[tnobs]; //vector for tree exercise, sorted
	anobs = 1236; //number of observations in ANOVA data
	nvars = 9; //number of variables in ANOVA data
	double matA[anobs*nvars]; //data for ANOVA
	double matAS[anobs*nvars]; //sorted data

	//Read in the anova data generated in R
		 //number of observations
		infile=fopen("babies_for_C.txt", "r");
			if (!infile) {
			printf("Cannot Open File");
			return 1;
		}
		for (i=0; i<anobs; i++) {
			for (n=0; n<nvars; n++) {
				fscanf(infile, "%lf", &myd); //fscanf reads row-by-row
				matA[i+anobs*n] = myd;
			}
		}
		fclose(infile);
		//printf("Babies Data Matrix (from file)\n");
		//printmat(matA,anobs,nvars);

	//Extract just 1000 obs from one of the columns for the tree exercise
		for (i=0; i<tnobs; i++) {
			for (n=0; n<nvars; n++) {
				if (n==0) matT[i] = matA[i+anobs*n];
			}
		}
		//printf("Birthweights for use with tree functions, unsorted \n");
		//printmat(matT,1,tnobs);

	//Read in the sorted list of 1000 obs for tree exercise
		/*to sort data from unix, type the following:
			cat original.txt | sort -n > newdata.txt */
		infile=fopen("babies_for_C_sorted.txt", "r");
			if (!infile) {
			printf("Cannot Open File");
			return 1;
		}
		for (i=0; i<anobs; i++) {
			for (n=0; n<nvars; n++) {
				fscanf(infile, "%lf", &myd); //fscanf reads row-by-row
				matAS[i+anobs*n] = myd;
			}
		}
		fclose(infile);
		for (i=0; i<tnobs; i++) {
			for (n=0; n<nvars; n++) {
				if (n==0) matTS[i] = matAS[i+anobs*n];
			}
		}
		//printf("Birthweights for use with tree functions, sorted \n");
		//printmat(matTS,1,tnobs);


if (run_Q1) {
printf("\n\n*************** Q1: BINARY TREES ****************\n\n");

	struct OBS *mypoint = NULL;
	struct OBS myobs;
	int mytimes=10000; //number of times to run tree algorithms (for timing purposes)
	int timeme;

	// Use binary trees to order and count factor levels -- use unsorted data
	TIME0
	for (timeme=0; timeme<mytimes; timeme++){
		char *myroot = NULL;
		for(i=0;i<tnobs;i++){
			myobs.value=matT[i];
			myobs.count=0;
			myobs.level=0;
			mypoint = myinsert(&myroot,&myobs,sizeof(myobs),comp);
			(*mypoint).count +=1;
		}
		mytraverse(myroot,leveltraverse);
		//mytraverse(myroot,printtraverse);
	}
	printf("\n");
	printf("Time using regular binary trees %d times, unsorted", timeme);
	TIME1(" ")
 	printf("\n \n");

	// Use balanced binary trees to order and count factor levels -- use unsorted data
	TIME0
	for (timeme=0; timeme<mytimes; timeme++){
		char *root = NULL;
		for(i=0;i<tnobs;i++){
			myobs.value=matT[i];
			myobs.count=0;
			myobs.level=0;
			mypoint = insert(&root,&myobs,sizeof(myobs),comp);
			(*mypoint).count +=1;
		}
		traverse(root,leveltraverse);
		//traverse(root,printtraverse);
	}
	printf("\n");
	printf("Time using balanced binary trees %d times, unsorted", timeme);
	TIME1(" ")
 	printf("\n \n");

	// Use binary trees to order and count factor levels -- use unsorted data
	TIME0
	for (timeme=0; timeme<mytimes; timeme++){
		char *myroot = NULL;
		for(i=0;i<tnobs;i++){
			myobs.value=matTS[i];
			myobs.count=0;
			myobs.level=0;
			mypoint = myinsert(&myroot,&myobs,sizeof(myobs),comp);
			(*mypoint).count +=1;
		}
		mytraverse(myroot,leveltraverse);
		//mytraverse(myroot,printtraverse);
	}
	printf("\n");
	printf("Time using regular binary trees %d times, sorted", timeme);
	TIME1(" ")
 	printf("\n \n");

	// Use balanced binary trees to order and count factor levels -- use unsorted data
	TIME0
	for (timeme=0; timeme<mytimes; timeme++){
		char *root = NULL;
		for(i=0;i<tnobs;i++){
			myobs.value=matTS[i];
			myobs.count=0;
			myobs.level=0;
			mypoint = insert(&root,&myobs,sizeof(myobs),comp);
			(*mypoint).count +=1;
		}
		traverse(root,leveltraverse);
		//traverse(root,printtraverse);
	}
	printf("\n");
	printf("Time using balanced binary trees %d times, sorted", timeme);
	TIME1(" ")
 	printf("\n \n");



} /*end run_Q1*/

if (run_Q2) {
printf("\n\n*************** Q2: ANOVA ****************\n \n");



} /*end run_Q2*/


}
