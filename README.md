# fragger
protein fragments picker

Fragger is a protein fragment picker allowing to create and query protein fragment databases.

All fragment lengths are supported and any set of PDB files can be used to create a database.

Fragger can efficiently search a fragment database with a query fragment and an RMSD threshold.
The query fragment can have structural gaps and the allowed amino acid sequences matching a query can be constrained via a regular expression of one-letter amino acid codes.

Fragger also incorporates a tool to compute the RMSD of one versus many fragments in high throughput.

[![DOI](https://zenodo.org/badge/30861587.svg)](https://zenodo.org/badge/latestdoi/30861587)

