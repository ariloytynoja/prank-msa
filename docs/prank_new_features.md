# PRANK: new features
### New features in the latest version of PRANK

[Back to PRANK home.](../README.md)  
 
#### v.170427

Changes by Nikola Hecker to make the temp files thread safe, a fix of a compile issue and reduce of STDOUT traffic. Probably lots more.

#### v.140603

A bug causing BppAncestor to crash due to weird sequence names fixed. Added codon ancestral reconstruction also for translated alignments. New option ‘-njgaps’ to consider gaps as mismatches for NJ distances.

#### v.140110

Search order of external tools changed and bundled programs now preferred. Ancestor inference for translated alignments of DNA sequences. A fix for rare crashes due to overly long branch lengths.

#### v.130820

Disabled computation of alignment score for the guide tree alignment.  
More information about iteration with a user-provided tree.  
Path to other binaries no more affected by renaming of the program file.

#### v.130708

Significant bug fixes — update \*strongly\* recommended:  
\* option -F was mistakenly turned off by the new iterative approach  
\* without option -F, ancestor reconstruction (and scoring) were incorrect

#### v.130410

New optimization score and iteration strategy. Ancestral character inference using maximum likelihood (!BppAncestor) and output of events per tree branch. Changes in program options.

#### v.121218

More detailed information about unmatching names. New option “-prunedata”.

#### v.121212

Support for some NHX tags. Mainly for internal usage now.

#### v.121210

Fixed underflow errors affecting ancestral reconstruction of large alignments.

#### v.121018

Ancestral sequences can differently indicate insertions and deletions.  
Can update an alignment, recomputing nodes with tag “\[&&NHX:XN=realign\]”

#### v. 121002

Alignment merge now accepts trees such as “(t1:#.#,t2:#.#);”. The two-node tree can be provided with “-t=filename” or “-tree=tree-string”.

#### v. 120814

Can now also merge “alignments” of one sequence. New option {{{-mergedist=#}}} to define the distance for two alignments.

#### v. 120717

All input data now converted to upper-case.

#### v. 120716

Fixed the translated alignment option that had been broken in recent clean up, and corrected the output order for ancestral sequences.

#### v. 120712

Bug fixes: the new speed-up with MAFFT often failed with codon sequence data. That has now been fixed. Also, for codon data, the guide tree is now computed from a MAFFT alignment of translated protein sequences.

#### v. 120626

The new version of PRANK automatically uses programs MAFFT and Exonerate to speed up the alignment. If these programs are available, the analysis starts with a quick MAFFT alignment and estimation of a NJ guide tree from that. Exonerate is used to anchor the pairwise alignments and thus massively reduce the computation time. The anchoring works for DNA, protein and codon alignments, the last one using translation to proteins in the anchoring step. The guide tree estimation and anchoring work automatically and can be disabled with options {{{-nomafft}}} and {{{-noanchors}}}, respectively.

In addition to the alignment anchoring, the new version also provides an option to merge two pre-existing alignments. When merging two alignments, PRANK first infers the ancestors and the sequence at the root for each of them and then aligns the two root sequences. The new option is based on a previous (undocumented) function to finish a partially aligned alignment. See the documentation page for details.

Also the previous versions of PRANK have been able to infer the ancestral sequences for a new or an existing alignment. However, the output format for this was somewhat unclear and a new option {{{-showanc}}} has now been added with a more straightforward output.

Finally, the default options of the program have been slightly altered. The old behaviour and output (e.g., guide tree and xml file) can be achieved with appropriate options. See \`prank -h\` for details.