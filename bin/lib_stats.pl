#!/bin/bash
for file in `cut -f1 db/db_initials` ; do echo $file; ./bin/count_fasta_seq_len.pl <db/$file |cut -f2|perl -e 'use POSIX; $sum=0; $n=0;while($r=<STDIN>) {chomp $r; $sum+=$r; $n++; push @val,$r;} if ($n%2==0){ $median = ($val[floor($n/2)-1]+$val[floor($n/2)])/2 } else { $median = $val[$n/2-1] } print "$sum\t$n\t".$sum/$n."\t$median\n";' ; done

