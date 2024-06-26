#!/usr/local/bin/perl -w
#
########################################
#
# File Name:
#   oneLineFa.pl
# 
# Description:
#   
# 
# Usage:
#   
# 
# Author:
#   Xi Wang, Xi.Wang@dkfz.de
#
########################################

use strict;
my $usage = "$0 <infile> <outfile>\n";
my $infile = shift || die $usage;
my $outfile = shift || die $usage;
open(IN, $infile) || die "Can't open $infile for reading!\n";
open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";

my $first = 1;
while(<IN>){
  if(/^>/) {
    my @a = split;
    if($first) {
      $first = 0; 
      print OUT $a[0]."\n";
    } else {
      print OUT "\n".$a[0]."\n";
    }
  }
  else {
    chomp;
    print OUT $_;
  }
}

close IN;
close OUT;
