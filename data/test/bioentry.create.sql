CREATE TABLE bioentry (
    bioentry_id integer NOT NULL,
    biodatabase_id integer NOT NULL,
    taxon_id integer,
    name character varying(40) NOT NULL,
    accession character varying(128) NOT NULL,
    identifier character varying(40),
    division character varying(6),
    description text,
    version integer NOT NULL
);
