DIFF = grep -v "^[ ]*$$" $@.out | diff $@.expect -

CHECK = if [ ` $(DIFF) | wc -l ` -eq 0 ]; then echo "*** $@ OK ***"; else echo "*** $@ FAILED ***"; echo "Here's a diff: "; $(DIFF); exit 1; fi

SQLITE_SELECT = sqlite3 $@.sqlite3 "select * from accessions order by accno; select * from bestscoring order by accno, profile; select * from hmm_profiles order by profile; select * from taxa order by ncbi_taxon_id;" > $@.out

all: hmmsearch2classification.00 hmmsearch2classification.01 hmmsearch2classification.02

Nrd.test.tar.gz:
	scp rnrdb:projects/profiles/profiles/Nrd.test.tar.gz .

hmmsearch2classification.00:
	../R/hmmsearch2classification.r --verbose --profilehierarchies=$@.phier.tsv --singletable=$@.out --taxflat=taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(CHECK)

hmmsearch2classification.01:
	../R/hmmsearch2classification.r --verbose --profilehierarchies=$@.phier.tsv --singletable=$@.out --taxflat=taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(CHECK)

hmmsearch2classification.02:
	../R/hmmsearch2classification.r --verbose --profilehierarchies=$@.phier.tsv --sqlitedb=$@.sqlite3 --taxflat=taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(SQLITE_SELECT)
	@$(CHECK)
