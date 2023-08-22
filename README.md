# PRANK

PRANK is a probabilistic multiple alignment program for DNA, codon and amino-acid sequences. It’s based on a novel algorithm that treats insertions correctly and avoids over-estimation of the number of deletion events. In addition, PRANK borrows ideas from maximum likelihood methods used in phylogenetics and correctly takes into account the evolutionary distances between sequences. Lastly, PRANK allows for defining a potential structure for sequences to be aligned and then, simultaneously with the alignment, predicts the locations of structural units in the sequences.

* * *

![](docs/data/prank_logo.png)

*   [PRANK binaries](binaries)

Packed source code and pre-compiled binaries for Linux, Windows and OSX.

*   [Using PRANK](#using-prank)

Instructions for the use of PRANK with example commands and test data.

*   [Installing PRANK](docs/prank_installation.md)

Instructions to install and start using PRANK on different operating systems. Pre-compiled executables are provided for Windows, Mac OSX and Linux; the source code should compile on any system supporting the GNU C compiler package.

*   [New features](docs/prank_new_features.md)

Brief description of the new features in the latest versions of the program.

*   [Known issues](docs/prank_known_issues.md)

Known bugs and problems with the latest versions of the program.

*   [Structure models](docs/prank_structure_models.md)

Generation of PRANK structure models.

*   [HSAML format](docs/hsaml_format.md)

Definition of the XML-format used by PRANK.


* * *

### Using PRANK

At the simplest, PRANK can be run with command:  
 
```
prank input_file
```

where ```input_file``` contains sequences in FASTA format.

See below for the description of the most central program options.

*   [Using PRANK](#using-prank)
*   [About PRANK alignments](#about-prank-alignments)
*   [Main program options](#main-program-options)
*   [Other program options](#other-program-options)
    *   [Translated alignment and codon alignment](#translated-alignment-and-codon-alignment)
    *   [Inference of ancestral sequences and ancestral events](#inference-of-ancestral-sequences-and-ancestral-events)
    *   [Finishing an alignment and merging two alignments](#finishing-an-alignment-and-merging-two-alignments)
    *   [Alignment reformatting and back-translation](#alignment-reformatting-and-back-translation)
*   [Bug reporting](#reporting)
*   [Methods](#methods)

### About PRANK alignments

PRANK aims at an evolutionarily correct alignment and the alignments inferred with PRANK can be expected to look different from ones generated with other alignment methods. There are, however, cases where the different look is caused by violations of the method’s assumptions. To understand why things may go wrong and how to avoid that, read this [explanation of differences](docs/prank_differences.md) between PRANK and traditional progressive alignment methods.

The reconstruction of evolutionary homology — including the correct placement of insertion and deletion events — is only feasible for rather closely-related sequences. PRANK is not meant for the alignment of very diverged sequences. If sequences are very different, the correct homology cannot be reconstructed with confidence and PRANK may simply refuse to match them. There are several methods developed for structural matching of very distant protein sequences. One should not consider the resulting alignment a proper inference of evolutionary homology, though.

Often there are ties in the alignment, i.e. positions with many equally good solutions. The common practice among sequence aligners is always to pick the same solution among many different ones and produce consistent alignments. We believe that this gives the user a wrong impression and too high confidence on the resulting alignment (for example, 10 iterations of the alignment always produce the same result -> the alignment must be correct). Our decision is to break the ties randomly. This means that different runs of the program with the same data may give different alignments and PRANK alignment may not be reproducible. However, the practice we have chosen also tells the user that the regions not consistently aligned in a similar manner simply cannot be resolved reliably. Be aware that the reproducability in many other methods does not mean higher confidence — they also reproduce exactly the same errors!

### Main program options

A typical command for PRANK could be:

```
prank -d=input_file -t=tree_file -o=output_file -F -showxml
```

Here, ```input_file``` is the name of the file with input sequences in FASTA format, ```tree_file``` is the name of the file with a phylogeny with branch lengths relating these sequences, and ```output_file``` is the name of the file (with extension .best.fas) where the resulting alignment will be written in FASTA format; ```-F``` specifies that the inference of insertions should be trusted and sites appearing as insertions should not be aligned at the later stages of the process; and ```-showxml``` specifies that the output is also written in HSAML format, compatible for Wasabi import (including the alignment phylogeny).

Of these parameters, only ```-d=input_file``` is strictly required (and, when given alone, can be given in the form ```prank input_file```). If the option ```-t=tree_file``` is omitted, an alignment guide tree is inferred from the input data (see below). If the option ```-o=output_file``` is omitted, the results are written to file named output.best.fas. Finally, if option ```-F``` is omitted, the algorithm does not incorrectly penalise for gaps for insertions but it will also not guarantee to leave the insertions unmatched.  
 

By default, PRANK iterates the alignment five times, inferring a new guide tree after each alignment. For each solution, it computes a phylogeny-aware alignments score and, after the full analysis, outputs the best-scoring solution. The number of iterations performed can be changed with option ```-iterate=#``` and the iteration completely disabled with option ```-once```. See [Methods](#methods) for more details.  
 

Without additional options, PRANK outputs only the alignment. More information can be outputted with the following options:  
 
```
-showtree: output the guide tree used  
-showxml: output the results in PRANK’s own XML format  
-showanc: output ancestral sequences and guide tree  
-showevents: output events per branch  
-showall: output all listed above  
-showiter: output results for each iteration  
```

If a guide tree is provided, PRANK assumes that the sequence names in the tree and data file match exactly. Special characters and white spaces in the names may cause problems and may have to be removed/replaced with underscores. With option ```-shortnames``` only the first part (until the first white space) of the names in the data file are matched against the names in the guide tree. Extra names in the guide tree/data file can be removed with options ```-prunetree``` and ```-prunedata```.  
 

A list of possible program options is shown with command:  
 
```
prank -help
```

By default, PRANK writes the output alignments in FASTA format. With the option \-f=format\_name some other popular formats can be selected. See prank -help for details.

[back to top](#using-prank)

### Other program options

#### Translated alignment and codon alignment

PRANK can do translated alignments of protein-coding DNA sequences or align them using the codon model. Translation is selected with the options ```-translate``` (standard code) or ```-mttranslate``` (mitochondrial code), and the codon alignment with the option ```-codon```. Using the example data [input_dna.fas](docs/data/input_dna.fas), the following commands make a translated alignment:  
 
```
prank -d=input_dna.fas -o=output_translated -translate -F
```

and a codon alignment:  
 
```
prank -d=input_dna.fas -o=output_codon -codon -F
```
   
See [Methods](#methods) for more details.

#### Inference of ancestral sequences and ancestral events

The key aim of the phylogeny-aware alignment algorithm is to correctly model the insertion and deletion events and accurately infer the ancestral sequences, especially the presence and absence of characters. The ancestral sequences, that PRANK internally will anyway infer, can be of great interest and can also be outputted. The earlier output format for the ancestral sequences was not clear and new option ```-showanc``` has been added. This option now inserts the ancestral sequences in the position they appear in the alignment tree structure and, when viewed together with a phylogenetic tree with internal nodes labelled accordingly, make it easy to assign evolutionary events to specific tree branches. The two output files are named ```*.best.anc.fas``` and ```*.best.anc.dnd```; with the option ```-showevents```, the inferred events are listed per branch and outputted to file ```*.best.events```.  
 

Ancestral sequences can be outputted for alignments generated with PRANK but they can also be inferred for existing alignments. While PRANK by default removes all the gap characters in the input data, this can be disabled with option ```-keep``` and it then reads the alignment as it is. Thus with command:  
 
```
prank -d=alignment_pep.fas -showanc -showevents -keep -o=output_ancestors -njtree [OR -t=treefile]
```

PRANK will infer the ancestral sequences for an exiting alignment (this one available at [alignment\_pep.fas](docs/data/alignment_pep.fas)).  
 

If the aim is infer ancestral codon sequences but the codon alignment itself is found too slow, one can make the alignment using the protein model and then infer the ancestors for the back-translated DNA alignment:  
 
```
prank -d=input.fas -o=output_translate -showall -translate
```

The ancestral sequences will be outputted to ```output_translate_codon.[nuc|pep\].anc.fas```.  
 

See [Methods](#methods) for more details.

#### Finishing an alignment and merging two alignments

PRANK has new features that allow both aligning two alignments and finishing an alignment of a partially resolved dataset, consisting of one of more aligned sub-trees.  
 

In the first case, the two alignments and the corresponding trees (for example, [alignment_1.fas](docs/data/alignment_1.fas), [alignment_2.fas](docs/data/alignment_2.fas), [tree_1.tre](docs/data/tree_1.tre) and [tree_2.tre](docs/data/tree_2.tre)) are defined as:  
 
```
prank -d1=alignment_1.fas -d2=alignment_2.fas -t1=tree_1.tre -t2=tree_2.tre -tree="(t1:0.1,t2:0.2);"
```

If no trees (```-t1=``` and ```-t2=```) are provided, PRANK will compute NJ trees for the two alignments and infer ancestors based on those. If no merge-tree (```-t=``` or ```-tree=```) is provided, the default branch length is used; that can be modified with option ```-mergedist=```. The merging of the two alignments is based on the alignment of the two ancestors.  
 

In the second case, the sequences belonging to the same group are indicated by ” group_XX” in the end of the name (where XX is unique for each group) and PRANK automatically finds the maximally large sub-trees consisting of sequences from the same group.  
 

This feature is best explained by an example. If you download files [example_partaligned.fas](docs/data/example_partaligned.fas) and [example.tre](docs/data/example.tre), the alignment can be finished with the command:  
 
```
prank -d=example_partaligned.fas -t=example.tre -partaligned
```

When finishing an alignment, PRANK looks for maximally large sub-trees of pre-aligned sequences belonging to the same group. If a sub-tree includes a sequence not from the same group, that sequence is normally aligned and also sequences outside that sequence (even if they belong to that group) will be re-aligned. This is again better explained by an example. If you download files [example_partaligned_2.fas](docs/data/example_partaligned_2.fas) and run the command:  
 
```
prank -d=example_partaligned_2.fas -t=example.tre -partaligned
```

you notice that many more nodes need to be aligned in order to finish the complete alignment. Note that PRANK will fail if the group annotation does not match actual data and sequences assigned to a specific group are not truly aligned.

#### Alignment reformatting and back-translation

PRANK supports several different alignment formats and can translate and back-translate sequence data between DNA and protein. These features can be exploited also without performing alignment of sequences. Example data used below can be found at [alignment\_codon.fas](docs/data/alignment_codon.fas), [alignment\_pep.fas](docs/data/alignment_pep.fas) and [input\_dna.fas](docs/data/input_dna.fas).  
 

Option ```-convert``` tells to convert between formats. This is typically accompanied with option ```-f=format``` where format is either ```fasta```, ```phylipi```, ```phylips```, ```paml```, ```nexus``` or ```raxml```. The input and output files are specified as normally with options ```-d=filename``` and ```-o=filename```; with options ```-translate``` and ```-mttranslate``` the DNA sequences are additionally translated to proteins. The command:  
 
```
prank -convert -d=alignment_codon.fas -translate -f=phylipi -o=alignment_pep -keep
```

converts a codon alignment in FASTA format to a protein alignment in interleaved PHYLIP format.  
 

If option ```-dna=filename``` is included, PRANK attempts to back-translate the input protein alignment to the corresponding DNA alignment. This assumes that the input alignment (```-d=filename```) is proteins and the the other file (```-dna=filename```) contains matching unaligned DNA sequences with identical names. The command:  
 
```
prank -convert -d=alignment_pep.fas -dna=input_dna.fas -o=alignment_dna -keep
```

converts a protein alignment in FASTA format to a DNA alignment in FASTA format.

[back to top](#using-prank)


Methods
-------

#### Models

For DNA data, PRANK by default uses HKY model with empirical base frequencies and kappa=2. With the optional command parameters, it supports TN [(Tamura and Nei, 1993)](http://mbe.oxfordjournals.org/content/10/3/512.short) and models below it (JC, K2P, FEL, HKY). For example, JC model is defined as ```-kappa=1 -dnafreqs=0.25,0.25,0.25,0.25```. WAG [(Whelan and Goldman, 2001)](http://mbe.oxfordjournals.org/content/18/5/691.short) is used for protein alignments.  
 

PRANK can do translated alignments of protein-coding DNA. It translates protein-coding DNA to protein sequences, aligns these sequences as proteins, and back-translates the resulting alignment to DNA such that the gaps are maintained. (PRANK can also back-translate protein alignments produced with external alignment software; see below for examples.)  
 

In addition to translated alignment, PRANK can also align codon sequences using a codon substitution matrix [(Kosiol, Holmes and Goldman, 2007)](http://mbe.oxfordjournals.org/content/24/7/1464.short). Translation into amino acids and codons is done in the first forward frame without any error-checking. Based on independent benchmarks, codon alignment produces more accurate alignments than alignment of translated protein sequences.  
 

Simulation studies with nucleotide sequences containing high numbers of insertions and deletions showed that the option ```-F``` gives the most accurate results and should be used when the guide phylogeny can be trusted [(Löytynoja and Goldman, 2008)](http://www.sciencemag.org/content/320/5883/1632.short).  
 

#### Guide tree

Progressive alignment requires a guide tree, and the algorithm of PRANK can be especially sensitive to errors in that. If a trusted pre-existing phylogeny with branch lengths is available for the input sequences (as often is in comparative genomic studies), it is recommended to use that.  
 

If no tree is provided, PRANK constructs one: it first call MAFFT [(Katoh et al, 2005)](http://nar.oxfordjournals.org/content/33/2/511.short) to make a quick alignment and infers an NJ tree from the evolutionary distances based on that; an alternative (and much slower) approach is to estimate the evolutionary distances from pairwise alignments generated by PRANK itself.  
 

By default, PRANK iterates the alignment five times and writes the solution with the best score to the file named filename.best.fas. To prevent alignment iteration (e.g. when providing a guide tree that you trust), use the flag ```-once``` or ```-iterate=1```.  
 

#### Anchoring

The standard PRANK algorithm is based on an exhaustive search of the best pairwise solution, and for long sequences this soon becomes too time consuming. By default PRANK uses Exonerate [(Slater and Birney, 2005)](http://www.biomedcentral.com/1471-2105/6/31) to anchor the pairwise alignments and thus speed up the process.  
 

#### Ancestral reconstruction

PRANK’s own ancestral reconstruction happens during progressive alignment and is based on descendant sequences only. To improve the ancestral reconstruction, PRANK uses BppAncestor from the BppSuite package [(Dutheil and Boussau, 2008)](http://www.biomedcentral.com/1471-2148/8/255) to infer the character states and then combines these with the inferred character presence/absence information. PRANK’s optimization score is also based on these maximum likelihood inferences.

[back to top](#using-prank)
