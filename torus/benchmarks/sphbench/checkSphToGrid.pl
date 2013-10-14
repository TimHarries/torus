#!/usr/bin/perl -w

# Check the mass on the AMR grid from gridding a disc of SPH particles

use strict;

my $massRequired=0.011;
# Based on figure 2 of Acreman, Harries and Rundle
my $tolerance=0.03;

# Central source of one solar mass
my $msol=1.9891e33;
my $massSource=$msol;
my $toleranceSource=0.001;

open RUNLOG, "<$ARGV[0]";

my @line;
my $mass;
my $fracDiff;
while(<RUNLOG>){

# Check gas mass
  if (/Mass of envelope:/){
    @line=split;
    $mass=$line[4];
    print "\nGas mass on grid is $mass solar masses \n";
    print "Expected value is $massRequired solar masses \n";
    $fracDiff= ($mass - $massRequired) / $massRequired;
    print "Fractional difference = $fracDiff \n";
    if ( abs($fracDiff) < $tolerance ){
      print "Gas mass agrees within tolerance of $tolerance\n\n"
    } else
      {
	print "TORUS: test failed (gas mass)\n";
	exit 1
      }
  }

# Check point source mass
    if (/Sink Particle/){
      @line=split;
      $mass=$line[7];
      print "Point mass is $mass solar masses \n";
      $fracDiff= ($mass - $massSource) / $massSource;
      print "Fractional difference = $fracDiff \n";
      if ( abs($fracDiff) < $toleranceSource ){
      print "Point mass agrees within tolerance of $toleranceSource\n\n"
    } else
      {
	print "TORUS: test failed (point source)\n";
	exit 1
      }
    }
}

  print "TORUS: Test successful\n";

close RUNLOG;

