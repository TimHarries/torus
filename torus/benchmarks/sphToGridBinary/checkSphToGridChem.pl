#!/usr/bin/perl -w

# Check the mass on the AMR grid from gridding HI and CO

use strict;

# These are taken from info_sph.dat
my $massRequired=1725900235.8178394;
my $tolerance=0.12;
my $COmassRequired=12660.161777656678;
my $COtolerance=0.10;

open RUNLOG, "<$ARGV[0]" or die "Cannot open log file";

my @line;
my $mass;
my $fracDiff;

while(<RUNLOG>){

# Check HI gas mass. With convertRhoHI true the reported envelope mass will be HI.
  if (/Mass of envelope:/){
    @line=split;
    $mass=$line[4];
    print "\nHI mass on grid is $mass solar masses \n";
    print "Expected value is $massRequired solar masses \n";
    $fracDiff= ($mass - $massRequired) / $massRequired;
    print "Fractional difference = $fracDiff \n";
    if ( abs($fracDiff) < $tolerance ){
      print "HI mass agrees within tolerance of $tolerance\n"
    } else
      {
	print "TORUS: test failed (HI mass)\n";
	exit 1
      }
  }

# Check CO gas mass.
  if (/Molecular mass of envelope:/){
    @line=split;
    $mass=$line[5];
    print "\nCO mass on grid is $mass solar masses \n";
    print "Expected value is $COmassRequired solar masses \n";
    $fracDiff= ($mass - $COmassRequired) / $COmassRequired;
    print "Fractional difference = $fracDiff \n";
    if ( abs($fracDiff) < $COtolerance ){
	print "CO mass agrees within tolerance of $tolerance\n";
	print "TORUS: Test successful\n";
    } else
      {
	print "TORUS: test failed (CO mass)\n";
	exit 1
      }
  }
}

close RUNLOG;

