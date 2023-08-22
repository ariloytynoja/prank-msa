#!/usr/bin/perl

# This file gives an example how to access and manipulate prank XML alignment
# files. The code given here is not professional and may contain bugs. In this
# example case it seems to work fine. The code can be freely copied and used
# elsewhere.

# See the documentation for libXML to learn more of the XML module.
# A tutorial is available at http://www.xml.com/pub/a/2001/11/14/xml-libxml.html

use XML::LibXML;

# Node.pm is a small module for manipulating the tree.
# Some functions are used and explained here, the rest should be self-explanatory

use Node;

# Here we call the XML module and read the file '20prim.xml'.

my $parser = XML::LibXML->new();
my $tree = $parser->parse_file("20prim.xml");
my $root = $tree->getDocumentElement;

# We get the newick tree.. 

$newick = ${$root->getElementsByTagName('newick')}[0]->getFirstChild->getData;

# .. the root node using Node.pm

my $rootNode = Node->new();
$rootNode = $rootNode->parseTree($newick);

# The tree is using id's, so we need a mapping between id's and real names.
# This mapping is stored as '%nameToLeaf' hash. Getting and storing the sequence 
# data is also trivial.

my %nameToLeaf = ();
my %seqsByName = ();
foreach my $lid (@{$root->getElementsByTagName('leaf')}) {
	$nameToLeaf{$lid->getAttribute('name')} = $lid->getAttribute('id');
	my $seq = $lid->findvalue('sequence');
	$seq =~ s/\s//g;
	$seqsByName{$lid->getAttribute('name')} = $seq;
}

# Posterior probabilities are named according to the nodes id and the process id.
# We store them in a hash ($postprob{$node}{$prob}) as raw data (string ot numbers).

%postprob = ();
foreach my $nid (@{$root->getElementsByTagName('node')}) {
	my $node = $nid->getAttribute('id');
	foreach my $pid (@{$nid->getElementsByTagName('probability')}) {
		my $prob = $pid->getAttribute('id');
		my $data = $pid->getFirstChild->getData;
		$data =~ s/\s//g;
		$postprob{$node}{$prob} = $data;
	}
}

# In order to use the human-readable process names, we need the mapping between
# these names and the model state id's.

%nameToState = ();
foreach my $sid (${$root->getElementsByTagName('model')}[0]->getElementsByTagName('probability')) {
	$nameToState{$sid->getAttribute('name')} = $sid->getAttribute('id');
}

# Let's define a group of nodes that end either with 'H.sapiens', 'P.paniscus' or 
# 'P.troglodytes' and then convert this to an array of sequence id's.

my @ingroup = grep /H.sapiens$/,keys %nameToLeaf;
push @ingroup, grep /P.paniscus$/,keys %nameToLeaf;
push @ingroup, grep /P.troglodytes$/,keys %nameToLeaf;

my @ingroupids = ();
foreach my $name(@ingroup) {
	push @ingroupids,$nameToLeaf{$name};
}

# Other functions for the manipulation of the tree.

#print join " ","rootNode->nodeNames",$rootNode->nodeNames,"\n";
#print join " ","rootNode->leafNames",$rootNode->leafNames,"\n";

#print "rootNode->treeDepth ",$rootNode->treeDepth,"\n"; 
#print "rootNode->leafNumber ",$rootNode->leafNumber,"\n";

# 'maxIngroup' gives the internal node that maximises the number of child
# nodes without including any nodes not mentioned on the list (it may not contain
# all the names in the list!). Here, we find the ancestor node for human and chimp
# sequences. The return values are the maximum number of leaf nodes included and 
# the name of the ancestor node maximising that.

(my $leafs,my $node) = $rootNode->maxIngroup(0,0,@ingroupids);
# print "rootNode->maxIngroup $leafs $node \n";

# One can do the same with names NOT included in the selection.

#($leafs,$node) = $rootNode->maxNotOutgroup(0,0,@ingroupids);


# As an example, we can find the sites for which the support at the nodes in the 
# human-chimp subtree is below 90%. To do that, we select 'postprob' for each 
# ancestral node and split the string into sitewise support values. First we call
# another function that return a list of nodes below the given node.

my @nodes = $rootNode->nodeNamesBelow($node,0);

my @sites = ();
my $threshold = 50;
foreach my $node(@nodes) {
	my @score = split /,/,$postprob{$node}{$nameToState{"postprob"}};
	for(my $i=0;$i<=$#score;$i++) {
		if($score[$i]<$threshold && $score[$i]>=0) {
			unless(grep /^$i$/,@sites){
				push @sites,$i;
			}
		}
	}
}

# Print the selected sites nicely (add 1 to correct for real numbering)...

print join " ","Sites in human-chimp subtree with support < $threshold%:\n";

@sites = sort numerically @sites;
for(my $i=0;$i<=$#sites;$i++) {
	if($i==0 || $sites[$i-1]+1<$sites[$i]){
		print " ",$sites[$i]+1;
		if(($i==0 || $sites[$i-1]+1<$sites[$i]) && $sites[$i]+1==$sites[$i+1]){
			print "-";
		}
		next;
	}
	if($sites[$i]+1<$sites[$i+1]){
		print $sites[$i]+1;
	}
}
print $sites[$#sites]+1,"\n";

# ..and print the sequence sites corresponding to the coordinates.

foreach my $name(@ingroup) {
	my $seq = $seqsByName{$name};
	my $p = -1;
	for(my $i=0;$i<=$#sites;$i++) {
		print " " if($p+1<$sites[$i]);
		$p = $sites[$i];
		print substr($seq,$sites[$i],1);
	}
	print "\n";
}



sub numerically { $a <=> $b; }