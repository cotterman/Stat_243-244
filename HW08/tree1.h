 struct MYNODE {
	struct MYNODE *left,*right;
   	char *udata;
        } ;

/* function prototypes for functions in btree1.c */

char *myinsert(char **pt,char *udata,int usize,int (*comp)());
void mytraverse(char *t,int (*func)());
void mydelete(char **pt,int (*delfn)());
