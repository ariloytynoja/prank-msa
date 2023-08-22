# PRANK: structure models

[Back to PRANK home.](../README.md)  

PRANK can align sequences using a “structure model” that describes multiple evolutionary processes and infer regions evolving under each process along with the alignment. PRANK imports these alignment models from structured flat files with parameter -m=model_name. 

## Format

This is a simple Two-state model describing regions of fast and slowly evolving sites. Comments are ignored by the program.

### Two-states-model
```
# Number of states   (obvious, isn't it? This is a two-state model!)
2
# Alphabet   (specifies DNA, RNA or protein -- or any other alphabet)
ACGT
# equilibrium distributions   (character frequencies for each state [1 and 2 here])
0.25 0.25 0.25 0.25 
0.25 0.25 0.25 0.25 
# Scaled Q matrices   (Q is instantaneous rate matrix, separately for each state)
# Model 1 ; scale = 0.855871   ('scale' can be seen as divider for branch length;
-1.1684 0.389467 0.389467 0.389467    as it is less than 1, branches are longer
0.389467 -1.1684 0.389467 0.389467    and sites are expected to evolve faster)
0.389467 0.389467 -1.1684 0.389467 
0.389467 0.389467 0.389467 -1.1684 
# Model 2 ; scale = 1.33333   ('scale' is greater than 1, so branches are shorter
-0.75 0.25 0.25 0.25            and sites are expected to evolve more slowly)
0.25 -0.75 0.25 0.25 
0.25 0.25 -0.75 0.25 
0.25 0.25 0.25 -0.75 
# Structure background probabilities   (starting probability -- equal so no info)
0.5 0.5
# Structure transition probabilities   (probability to move from one structure to
0.9667 0.0333                         another or to stay in a structure [diagonal])
0.0500 0.9500
# Indel rates   (indel rates in the two structures; 2nd seems to have fewer indels)
0.0500 0.0250
# Indel extension probabilities   (indel lengths; 2nd seems to have shorter indels)
0.8000 0.5000
# Match extension probabilities   (for completeness; better kept 0?)
0 0
# Codon position   (two structures are non-coding; otherwise 1, 2 or 3)
0 0
```


Below are examples of a Two-states+Codon model describing the fast and slow regions and the periodicity of the three codon sites, and a a Four-state model.

### Two-states+Codon-model
```
# Number of states
5
# Alphabet
ACGT
# equilibrium distributions
0.25 0.25 0.25 0.25 
0.25 0.25 0.25 0.25 
0.262295 0.262295 0.262295 0.213115 
0.229508 0.262295 0.245902 0.262295 
0.229508 0.262295 0.245902 0.262295 
# Scaled Q matrices
# Model 1 ; scale = 0.855871
-1.1684 0.389467 0.389467 0.389467 
0.389467 -1.1684 0.389467 0.389467 
0.389467 0.389467 -1.1684 0.389467 
0.389467 0.389467 0.389467 -1.1684 
# Model 2 ; scale = 1.33333
-0.75 0.25 0.25 0.25 
0.25 -0.75 0.25 0.25 
0.25 0.25 -0.75 0.25 
0.25 0.25 0.25 -0.75 
# Model 3; scale = 1 after codon adjustment
# Codon model for position[1]: omega = 0.5
-0.468687 0.179497 0.159553 0.129637 
0.179497 -0.488631 0.159553 0.149581 
0.159553 0.159553 -0.448743 0.129637 
0.159553 0.184099 0.159553 -0.503205 
# Model 4; scale = 1 after codon adjustment
# Codon model for position[2]: omega = 0.5
-0.478659 0.159553 0.159553 0.159553 
0.139609 -0.448743 0.149581 0.159553 
0.148916 0.159553 -0.468022 0.159553 
0.139609 0.159553 0.149581 -0.448743 
# Model 5; scale = 1 after codon adjustment
# Codon model for position[3]: omega = 0.5
-0.831954 0.262123 0.307709 0.262123 
0.229357 -0.77782 0.229357 0.319106 
0.287195 0.244648 -0.776491 0.244648 
0.229357 0.319106 0.229357 -0.77782 
# Structure background probabilities
0.5 0.5 0 0 0
# Structure transition probabilities
0.9667 0.0333 0 0 0
0.0487 0.9513 0.0254 0 0
0 0 0 1 0
0 0 0 0 1
0 0.0600 0.9400 0 0
# Indel rates
0.0500 0.0250 0.0250 0 0
# Indel extension probabilities
0.8000 0.5000 0 0 0.0323
# Match extension probabilities
0 0 0 0 0
# Codon position
0 0 1 2 3
# Draw pattern
2 2 2 2 2
# Draw color
2 1 3 3 4
# Draw offset
1 1 1 1 1
# State name
process1 process2 codon1 codon2 codon3
```

### Four-states-model

```
# Number of states
4
# Alphabet
ACGT
# equilibrium distributions
0.25 0.25 0.25 0.25 
0.25 0.25 0.25 0.25 
0.25 0.25 0.25 0.25 
0.25 0.25 0.25 0.25 
# Scaled Q matrices
# Model 1 ; scale = 0.568958
-1.7576 0.585867 0.585867 0.585867 
0.585867 -1.7576 0.585867 0.585867 
0.585867 0.585867 -1.7576 0.585867 
0.585867 0.585867 0.585867 -1.7576 
# Model 2; scale = 1
-1 0.333333 0.333333 0.333333 
0.333333 -1 0.333333 0.333333 
0.333333 0.333333 -1 0.333333 
0.333333 0.333333 0.333333 -1 
# Model 3 ; scale = 1.33333
-0.75 0.25 0.25 0.25 
0.25 -0.75 0.25 0.25 
0.25 0.25 -0.75 0.25 
0.25 0.25 0.25 -0.75 
# Model 4 ; scale = 2
-0.5 0.166667 0.166667 0.166667 
0.166667 -0.5 0.166667 0.166667 
0.166667 0.166667 -0.5 0.166667 
0.166667 0.166667 0.166667 -0.5 
# Structure background probabilities
0.25 0.25 0.25 0.25
# Structure transition probabilities
0.8000 0.0667 0.0667 0.0667
0.0667 0.8000 0.0667 0.0667
0.0667 0.0667 0.8000 0.0667
0.0667 0.0667 0.0667 0.8000
# Indel rates
0.0500 0.0250 0.0250 0.0250
# Indel extension probabilities
0.8000 0.8000 0.5000 0.5000
# Match extension probabilities
0 0 0 0
# Codon position
0 0 0 0
# Draw pattern
2 2 2 2
# Draw color
2 1 3 4
# Draw offset
1 1 1 1
# State name
process1 process2 process3 process4
```