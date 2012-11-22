#!/usr/bin/perl -w

# Check the mass on the AMR grid from gridding a disc of SPH particles

use strict;

my $massRequired=0.011;
# Based on figure 2 of Acreman, Harries and Rundle
my $tolerance=0.03;

open RUNLOG, "<$ARGV[0]";

my @line;
my $mass;
my $fracDiff;
while(<RUNLOG>){
  if (/Mass of envelope:/){
    @line=split;
    $mass=$line[4];
    print "Mass on grid is $mass solar masses \n";
    print "Expected value is $massRequired solar masses \n";
    $fracDiff= ($mass - $massRequired) / $massRequired;
    print "Fractional difference = $fracDiff \n";
    if ( $fracDiff < $tolerance ){
      print "TORUS: Test successful\n"}
    else
      {
      print "TORUS: Test failed\n"}
  }
}

close RUNLOG;

