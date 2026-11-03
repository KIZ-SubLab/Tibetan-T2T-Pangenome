#!/usr/bin/perl 
use strict;
use warnings;

my ($all, $my, $out) = @ARGV;
my %hash;
open IN1, "<$all";
open IN2, "<$my";
open OUT, ">$out";

while(<IN1>){
		chomp;
    my @line = split/\t/, $_;
    $hash{$line[0].$line[1]} = join "\t", @line[3..$#line];
}

while(<IN2>){
    chomp;
    my @row = split/\t/, $_;
	  if (exists $hash{$row[0].$row[1]}){
	      print OUT "chr$row[0]\t$row[1]\t$row[0]_$row[1]\t$hash{$row[0].$row[1]}\n";
    }
		else {
		    next;}
}
close IN1;
close IN2;
close OUT;
