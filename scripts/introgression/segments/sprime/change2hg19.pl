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
    $hash{$line[2].$line[3]} = "$line[0]\t$line[1]";
}

print OUT "CHROM\tPOS\tID\tREF\tALT\tSEGMENT\tALLELE\tSCORE\n";
while(<IN2>){
    chomp;
		my @row = split/\t/, $_;
	  if (exists $hash{$row[0].$row[1]}){
	      print OUT "$hash{$row[0].$row[1]}\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\n";
    }
		else {
		    next;}
}
close IN1;
close IN2;
close OUT;
