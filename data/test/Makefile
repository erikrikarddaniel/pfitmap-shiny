all_sql3dumps: bioentry.data.sql3dump taxon.data.sql3dump taxon_name.data.sql3dump

all_pgdumps: taxon.data.pgdump.gz taxon_name.data.pgdump.gz taxon.schema.pgdump.gz taxon_name.schema.pgdump.gz

%.data.pgdump.gz:
	pg_dump --inserts --data-only -t $(basename $(basename $(basename $@))) $$DB | gzip -c > $@

%.schema.pgdump.gz:
	pg_dump --schema-only -t $(basename $(basename $(basename $@))) $$DB | gzip -c > $@

%.sql3dump: %.pgdump.gz
	echo "BEGIN;" > $@
	gunzip -c $< | sed '/^SET/d' | sed '/SELECT pg_catalog.setval/d' | sed "s/true/'t'/g" | sed "s/false/'f'/g" >> $@
	echo "END;" >> $@
