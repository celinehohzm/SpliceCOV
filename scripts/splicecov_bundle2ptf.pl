#!/usr/bin/perl

use strict;

my ($bundlefile) = @ARGV;

my $chr;

open(F, $bundlefile) or die "Cannot open $bundlefile: $!";
while (<F>) {
    chomp;
    my @a = split;
    if ($a[0] eq "bundle") {
        $chr = $a[1];
    }
    elsif ($a[0] eq "tstart") {
        my $new_value = 1 - $a[2];
        print $chr, "\t", $a[1], "\t.\tTSS\t", sprintf("%.4f", $new_value), "\t", join("\t", @a[3..$#a]), "\n";
    }
    elsif ($a[0] eq "tend") {
        my $new_value = 1 - $a[2];
        print $chr, "\t", $a[1], "\t.\tCPAS\t", sprintf("%.4f", $new_value), "\t", join("\t", @a[3..$#a]), "\n";
    }
    elsif ($a[0] eq "jstart") {
        my $max_value = -1;
        my $max_field = '';
        my $maxsign = '.';
        for (my $i = 5; $i < @a - 2; $i++) {
            my @b = split(/:/, $a[$i]);
            my $value = $b[-1];  # The last element is the coverage value
            if ($value > $max_value) {
                $max_value = $value;
                $max_field = $a[$i];
                $maxsign = $b[1];
            }
        }
        my $new_value = 1 - $a[2];
        print $chr, "\t", $a[1], "\t", $maxsign, "\tJSTART\t", sprintf("%.4f", $new_value), "\t", $a[3], "\t", $a[4], "\t", $max_field, "\t", $a[-2], "\t", $a[-1], "\n";
    }
    elsif ($a[0] eq "jend") {
        my $max_value = -1;
        my $max_field = '';
        my $maxsign = '.';
        for (my $i = 5; $i < @a - 2; $i++) {
            my @b = split(/:/, $a[$i]);
            my $value = $b[-1];  # The last element is the coverage value
            if ($value > $max_value) {
                $max_value = $value;
                $max_field = $a[$i];
                $maxsign = $b[1];
            }
        }
        my $new_value = 1 - $a[2];
        print $chr, "\t", $a[1], "\t", $maxsign, "\tJEND\t", sprintf("%.4f", $new_value), "\t", $a[3], "\t", $a[4], "\t", $max_field, "\t", $a[-2], "\t", $a[-1], "\n";
    }
    else {
        print STDERR "Line error: $_\n";
        exit;
    }
}
close(F);
