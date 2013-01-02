#!/usr/bin/perl

while (<>) {
	chomp;
	s/^ +//;
	@arr = split /\s+/;
	@arr = split /->/,$arr[0];
	print "$arr[1] ";
}
