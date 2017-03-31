DIFF = grep -v "^[ ]*$$" $@.out | diff $@.expect -

CHECK = if [ ` $(DIFF) | wc -l ` -eq 0 ]; then echo "*** $@ OK ***"; else echo "*** $@ FAILED ***"; echo "Here's a diff: "; $(DIFF); exit 1; fi

all: hmmsearch2classification.00

Nrd.test.tar.gz:
	scp rnrdb:projects/profiles/profiles/Nrd.test.tar.gz .

hmmsearch2classification.00:
	../R/hmmsearch2classification --verbose $@.d/*.tblout > $@.out
	@$(CHECK)
