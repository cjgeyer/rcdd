
%.Rout : %.R
	R CMD BATCH --vanilla --debugger=valgrind $<

all : scdd.Rout redund.Rout linear.Rout
