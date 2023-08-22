package Node;

# This file gives help functions to access and manipulate prank XML alignment
# files. The code given here is not professional and may contain bugs. In the
# examples tested,the code seems to work fine. The code can be freely copied and
# used elsewhere.

sub new {
	my $self  = {};
	$self->{NAME}   = undef;
	$self->{LEFT}    = undef;
	$self->{RIGHT}  = undef;
	bless($self);
	return $self;
}

sub name {
	my $self = shift;
	if (@_) { $self->{NAME} = shift }
	return $self->{NAME};
}

sub left {
	my $self = shift;
	if (@_) { $self->{LEFT} = shift }
	return $self->{LEFT};
}

sub right {
	my $self = shift;
	if (@_) { $self->{RIGHT} = shift }
	return $self->{RIGHT};
}

sub nodeNames {
	my $self = shift;
	my @nodes = ();
	if($self->left) {
		push @nodes,$self->left->nodeNames;
	}
	if($self->name =~ /#\d+#/) {
		push @nodes,$self->name;
	}
	if($self->right) {
		push @nodes,$self->right->nodeNames;
	}
	return @nodes;
}

sub nodeNamesBelow {
	my $self = shift;
	my $anc = shift;
	my $found = shift;

	$found = 1 if($self->name =~ /^$anc$/);

	my @nodes = ();
	if($self->left) {
		push @nodes,$self->left->nodeNamesBelow($anc,$found);
	}
	if($found && $self->name =~ /#\d+#/) {
		push @nodes,$self->name;
	}
	if($self->right) {
		push @nodes,$self->right->nodeNamesBelow($anc,$found);
	}
	return @nodes;
}

sub leafNames {
	my $self = shift;
	if($self->name =~ /#\d+#/) {
		return ($self->left->leafNames,$self->right->leafNames);
	} else {
		return $self->name;
	}
}

sub treeDepth {
	my $self = shift;
	my $td = 0;
	if($self->left) {
		$td = $self->left->treeDepth + 1;
	}
	if($self->right) {
		my $rd = $self->right->treeDepth + 1;
		$td = $rd if($rd>$td);
	}
	return $td;
}

sub leafNumber {
	my $self = shift;
	if($self->name =~ /#\d+#/) {
		return ($self->left->leafNumber + $self->right->leafNumber);
	} else {
		return 1;
	}
}

sub maxIngroup {
	my $self = shift;
	my $mln = shift;
	my $mnn = shift;
	my @yes = @_;
	if($self->left) {
		($lln,$lnn) = $self->left->maxIngroup($mln,$mnn,@yes);
		if($lln > $mln) {
			$mln = $lln;
			$mnn = $lnn;
		}
	}
	if($self->right) {
		($rln,$rnn) = $self->right->maxIngroup($mln,$mnn,@yes);
		if($rln > $mln) {
			$mln = $rln;
			$mnn = $rnn;
		}
	}
	if($self->leafNumber > $mln){
		@below = $self->leafNames;
		$fine = 1;
		foreach $n(@below) {
			unless(grep /^$n$/,@yes){
				$fine = 0;
				last;
			}
		}
		if($fine) {
			$mln = $self->leafNumber;
			$mnn = $self->name;
		}
	}
	return ($mln,$mnn);
}

sub maxNotOutgroup {
	my $self = shift;
	my $mln = shift;
	my $mnn = shift;
	my @not = @_;
	if($self->left) {
		($lln,$lnn) = $self->left->maxNotOutgroup($mln,$mnn,@not);
		if($lln > $mln) {
			$mln = $lln;
			$mnn = $lnn;
		}
	}
	if($self->right) {
		($rln,$rnn) = $self->right->maxNotOutgroup($mln,$mnn,@not);
		if($rln > $mln) {
			$mln = $rln;
			$mnn = $rnn;
		}
	}
	if($self->leafNumber > $mln){
		@below = $self->leafNames;
		$fine = 1;
		foreach $n(@not) {
			$fine = 0 if(grep /^$n$/,@below);
		}
		if($fine) {
			$mln = $self->leafNumber;
			$mnn = $self->name;
		}
	}
	return ($mln,$mnn);
}

sub parseTree {
	my $self = shift;
	my $newick = shift;

	my $rootNode = 0;
	my %nodes = {};
	
	while((my $pair) = ($newick =~ /(\([\w#]+:\d+.\d+,[\w#]+:\d+.\d+\)#\d+#)/)){
		(my $n1,my $d1,my $n2,my $d2,my $p) = ($pair =~ /\(([\w#]+):(\d+.\d+),([\w#]+):(\d+.\d+)\)(#\d+#)/);
		$newick =~ s/\($n1:$d1,$n2:$d2\)//;
	
		$pn = Node->new();
		$pn->name($p);
		$nodes{$p} = $pn;
	
		if($nodes{$n1}) {
			$pn->left($nodes{$n1});
		} else {
			$ln = Node->new();
			$ln->name($n1);
			$pn->left($ln);
		}
		if($nodes{$n2}) {
			$pn->right($nodes{$n2});
		} else {
			$rn = Node->new();
			$rn->name($n2);
			$pn->right($rn);
		}
	
		$rootNode = $pn;
	}
	return $rootNode;
}

1;
