#!/usr/bin/perl

open(PROC,"cat mouse.gff3 |grep -v primary_trans|");
while(<PROC>){
	print "$_\n";
	@array=split("\t",$_);

	if (/Name=([^;]+)/){
		$name=$1;
	}

	$length=$array[4]-$array[3];
	print "$length $name\n";
}
