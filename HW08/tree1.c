#include <stdio.h>
#include "tree1.h"

/****************************************************************
                REGULAR BINARY TREE ROUTINE
				(code inspired by/stolen from btree1.c)
*********************************************************************/

static struct MYNODE *ret;

char *myinsert(char **pt,char *udata,int usize,int (*comp)())
{
  struct MYNODE *mydo_insert();
  static short zero = 0;

 return((char*)((mydo_insert((struct MYNODE **)pt,udata,usize,comp))->udata));
}


struct MYNODE *mydo_insert(struct MYNODE **pt,char *udata,int usize, int (*comp)())
{
   struct MYNODE *t1,*t2;
   int cc,a,i;


   if(*pt == NULL){
	if((ret = *pt = (struct MYNODE*)calloc(1,sizeof(struct MYNODE))) == NULL){
		fprintf(stderr,"insert: No memory available.  Exiting ...\n");
		exit(1);
		}

	if(((*pt)->udata = malloc((unsigned)usize)) == NULL){
		fprintf(stderr,"insert: No memory available.  Exiting ...\n");
		exit(1);
		}

	for(i=0;i<usize;i++)(*pt)->udata[i] = udata[i]; //this copies, byte by byte, member udata of whatever pt points to into udata
	return(ret);
	}

   cc = (*comp)((*pt)->udata,udata);


   if(cc > 0){
	(void)mydo_insert(&((*pt)->right),udata,usize,comp);
	return(ret);
	}
   else if(cc < 0){
	(void)mydo_insert(&((*pt)->left),udata,usize,comp);
	return(ret);
   	}
   ret = *pt;
   return(ret);
}

void mytraverse(char *t,int (*func)()) //this is a "depth-first" tranversal
{
  struct MYNODE *s = (struct MYNODE*)t; //s is a struct that was passes in to traverse

  if(s->left != NULL)mytraverse((char*)s->left,func); //if there is no left-child, then call traverse with member left of struct s
  (*func)((char*)(s->udata)); //will need a function that takes udata as the parameter
  if(s->right != NULL)mytraverse((char*)s->right,func);
}

void mydelete(char **pt,int (*delfn)())
{

  struct MYNODE **pn = (struct MYNODE **)pt;

  if(*pt){
	mydelete((char**)&((*pn)->left),delfn);
	mydelete((char**)&((*pn)->right),delfn);
	if(delfn)(*delfn)((*pn)->udata);
	free((char*)(*pn)->udata);
	free(*pt);
	*pt = NULL;
	}
}
