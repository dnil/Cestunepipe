#!/usr/bin/perl -w

use Bio::Tree::TreeI;
use Bio::Tree::NodeI;
use Bio::TreeIO;

my $DEBUG = 0;
my $VERBOSE = 1;

if(@ARGV < 3) {
    print "get_subtrees.pl treefile selectlist outfile\n";
    exit;
}

my $treefile = $ARGV[0];
my $subtree_selectlist = $ARGV[1];
my $outfile = $ARGV[2];

$DEBUG && print "arg0 $treefile\narg1 $subtree_selectlist\narg2 $outfile\n";

open(SELECT, $subtree_selectlist);
while ($select = <SELECT>) {
    chomp $select;
    $select =~ s/\|/\\|/g;
    push @subtreeid, $select;
}

$VERBOSE && print "Select subtree containing leaves with ids matching ", join " ", @subtreeid, "\n\n";

my $treeio = new Bio::TreeIO('-format' => 'newick',
			     '-file'   => $treefile);

my $tree = $treeio->next_tree;

$VERBOSE && print "Input tree is ", $tree->number_nodes, " nodes large (".scalar($tree->get_leaf_nodes)." leaves).\n\n";

# choose the subnodes for finding the lca somehow

my @nodes = $tree->get_nodes();
my @selectnodes;

foreach my $node (@nodes) {
    #check
    foreach my $subtreeid (@subtreeid) {
	if($node->id =~ /$subtreeid/) {
	    push @selectnodes, $node;
	    $VERBOSE && print "Select for $subtreeid: ", $node->id, "\n";
	}
    }
}

my @minimum_selected_node_distance;

for (my $a = 0; $a<@selectnodes-1; $a++) {
    
    my $select_node_a = $selectnodes[$a];    
    $$minimum_selected_node_distance[$a][$b] = -1;
    
    for my $b= $a+1; $b<@selectnodes; $b++) {
    $select_node_b (@selectnodes) {
	my $dist= $tree->distance([$select_node_a, $select_node_b]);
	if($$minimum_selected_node_distance[$a][$b] == -1) {
	    
	    $$minimum_selected_node_distance[$a][$b] = $dist;

	} elsif( $dist < $$minimum_selected_node_distance[$a][$b] ) {
	    $$minimum_selected_node_distance[$a][$b] = $dist;
	} elsif( $dist == $$minimum_selected_node_distance[$a][$b]) {
	    # surely useful..
	}

    }
    
    
}

# get the least common ancestor
my $lcanode = $tree->get_lca( -nodes => \@selectnodes );

# build a new tree with its root at $node
my $subtree = Bio::Tree::Tree->new(-root => $lcanode, -nodelete => 1);

$VERBOSE && print "\nSelected subtree is ", $subtree->number_nodes, " nodes large (".scalar($subtree->get_leaf_nodes)." leaves):\n\n";

if($DEBUG)
{ 
    my $out = Bio::TreeIO->new('-format'=>'newick');
    $out->write_tree($subtree);
}

open OUTFILE, ">$outfile";
for my $node ( $subtree->get_leaf_nodes ) {
    print OUTFILE $node->id, "\n";
    $VERBOSE && print $node->id, "\n";
}
close OUTFILE;
