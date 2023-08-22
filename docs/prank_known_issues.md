# PRANK: known issues
### Known bugs and problems with the latest version of the program

[Back to PRANK home.](../README.md)  


– PRANK infers ancestral character states with program BppAncestor. That program seems to fail with very large datasets (more than 500 sequences) and causes a premature termination of the run. The use of BppAncestor can be disabled with option -nobppa. If one does not care about ancestral states and alignment score, using this option makes the analysis slightly faster.