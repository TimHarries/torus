#!/usr/bin/perl -w

# Check the mass on the AMR grid from gridding HI and CO

use strict;

print "Checking mass values \n\n";

my $tolerance=0.12;
my $TRAPtolerance=0.30; 
my $COtolerance=0.10;

# Read in required masses from info_sph.dat
open SPHINFO, "info_sph.dat" or die "Cannot open info_sph.dat";

my @line;
my $massRequired;
my $COmassRequired;
while(<SPHINFO>){

    if (/Total HI mass/){
	@line=split;
	$massRequired=$line[4];
    }

    if (/Total molecular mass/){
	@line=split;
	$COmassRequired=$line[4];
    }

}
close SPHINFO;

print "Mass required is $massRequired \n";
print "CO mass required is $COmassRequired \n\n";;

# Parse log file to determine mass values from the AMR grid

open RUNLOG, "<$ARGV[0]" or die "Cannot open log file";

my $mass;
my $fracDiff;

while(<RUNLOG>){

# Check HI gas mass. With convertRhoHI true the reported envelope mass will be HI.
  if (/Mass of envelope:/){
    @line=split;
    $mass=$line[4];
    print "HI mass on grid is $mass solar masses \n";
    print "Expected value is $massRequired solar masses \n";
    $fracDiff= ($mass - $massRequired) / $massRequired;
    print "Fractional difference = $fracDiff \n";
    if ( abs($fracDiff) < $tolerance ){
      print "HI mass agrees within tolerance of $tolerance\n\n"
    } else
      {
	print "TORUS: test failed (HI mass)\n";
	exit 1
      }
  }

  if (/TRAP/){
    @line=split;
    $mass=$line[5];
    print "HI mass on grid (TRAP) is $mass solar masses \n";
    print "Expected value is $massRequired solar masses \n";
    $fracDiff= ($mass - $massRequired) / $massRequired;
    print "Fractional difference = $fracDiff \n";
    if ( abs($fracDiff) < $TRAPtolerance ){
      print "HI mass (TRAP) agrees within tolerance of $TRAPtolerance\n\n"
    } else
      {
	print "TORUS: test failed (HI mass TRAP)\n";
	exit 1
      }
  }



# Check CO gas mass.
  if (/Molecular mass of envelope:/){
    @line=split;
    $mass=$line[5];
    print "CO mass on grid is $mass solar masses \n";
    print "Expected value is $COmassRequired solar masses \n";
    $fracDiff= ($mass - $COmassRequired) / $COmassRequired;
    print "Fractional difference = $fracDiff \n";
    if ( abs($fracDiff) < $COtolerance ){
	print "CO mass agrees within tolerance of $tolerance\n\n";
	print "TORUS: Test successful\n";
    } else
      {
	print "TORUS: test failed (CO mass)\n";
	exit 1
      }
  }
}

close RUNLOG;

