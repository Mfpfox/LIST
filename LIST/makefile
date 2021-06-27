all:
	g++ -O3 source/LIST.c source/*.cpp -o list
	g++ -O3 source/Process_DB.c -o process_db
unpack:
	-gunzip UniProt/uniprot_*.fasta.gz
	-tar -zxvf Taxa/taxdump.tar.gz -C Taxa/
	-cat Taxa/nodes.dmp| awk '{print $$1, $$3}' > Taxa/tID_ptID.txt
	./process_db UniProt/uniprot_sprot.fasta UniProt/uniprot_trembl.fasta 
