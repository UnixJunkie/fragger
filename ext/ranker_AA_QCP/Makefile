all:    tags
	./build.sh

tags:
	find . -name "*.cc" >  src_files
	find . -name "*.h"  >> src_files
	cat src_files | xargs etags -a -o src/TAGS
	rm -f src_files

clean:
	rm -rf current.out TAGS build very_few_pdbs_rmsd_current out *.o \
	ranker_aa
