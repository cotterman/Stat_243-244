CFLAGS=-g

hw3: armalk.o as154.o arma_amoeba.o arma_powell.o arma_linmin.o arma_brent.o arma_mnbrak.o nrutil.o Stat244_hw3_main_v2.o drfftb.o drffti.o drfftf.o drftb1.o drfti1.o drftf1.o dradb2.o dradb3.o dradb4.o dradb5.o dradf2.o dradf3.o dradf4.o dradf5.o dradfg.o dradbg.o 
#cfftf.o cfftb.o cffti.o cfftf1.o cfftb1.o cffti1.o rfftf.o rfftb.o rffti.o rfftf1.o rfftb1.o rffti1.o
	cc -o hw3 armalk.o as154.o arma_amoeba.o arma_powell.o arma_linmin.o arma_brent.o arma_mnbrak.o nrutil.o Stat244_hw3_main_v2.o drfftb.o drffti.o drfftf.o drftb1.o drfti1.o drftf1.o dradb2.o dradb3.o dradb4.o dradb5.o dradf2.o dradf3.o dradf4.o dradf5.o dradfg.o dradbg.o 
	#cfftf.o cfftb.o cffti.o cfftf1.o cfftb1.o cffti1.o rfftf.o rfftb.o rffti.o rfftf1.o rfftb1.o rffti1.o
	
# 1st line says that armalk.o as154.o Stat244_hw3_main_v1.o etc. must exist in order to execute hw3
# 2nd line says to create the executable hw3 that consists of armalk.o as154.o Stat244_hw3_main_v1.o etc.
# this file will also compile each of the corresponding .c or .f files to create the .o files
# to have these lines executed, type make in the command prompt (PC will by default look for and run the file called makefile)
# if you want to both create the hw3 executable and also run it, you can type make && ./hw3

clean:
	rm -rf *.o hw3