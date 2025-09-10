#!/usr/bin/perl
use strict;
use warnings;

# Fast version for sorted junction BED:
# Input columns: chr, start, end, name, count, strand(+|-|.)
# Output: track header, then BED with last field "percs-percd-percl-percr" (4 decimals)

# ------------- helpers -------------
sub strand_idx {
    my ($s) = @_;
    return $s eq '+' ? 1 : $s eq '-' ? 0 : 2;  # 0:neg, 1:pos, 2:neutral
}

sub build_logs {
    my ($n) = @_;
    my @lg = (0) x ($n+1);
    for (my $i=2; $i <= $n; $i++) { $lg[$i] = $lg[$i>>1] + 1 }
    return \@lg;
}

# Build sparse table for RMQ (max), returns [table, logs]
sub build_sparse {
    my ($arr) = @_;
    my $n = scalar(@$arr);
    my $lg = build_logs($n);
    my $K = $$lg[$n];
    my @st;
    $st[0] = [ @$arr ];
    for (my $k=1; $k <= $K; $k++) {
        my $len = 1 << $k;
        my $half = $len >> 1;
        my $prev = $st[$k-1];
        my @cur;
        for (my $i=0; $i + $len - 1 < $n; $i++) {
            my $a = $prev->[$i];
            my $b = $prev->[$i + $half];
            $cur[$i] = ($a > $b) ? $a : $b;
        }
        $st[$k] = \@cur;
    }
    return (\@st, $lg);
}

# Query max on inclusive range [l, r] using sparse table
sub rmq_max {
    my ($st, $lg, $l, $r) = @_;
    return 0 if $l > $r;  # empty
    my $len = $r - $l + 1;
    my $k = $lg->[$len];
    my $a = $st->[$k]->[$l];
    my $b = $st->[$k]->[$r - (1<<$k) + 1];
    return ($a > $b) ? $a : $b;
}

# recompute max key from a freq hash (when current max becomes empty)
sub recompute_max_from_freq {
    my ($freq) = @_;
    my $m = 0;
    for my $k (keys %$freq) {
        next unless $freq->{$k} > 0;
        $m = $k if $k > $m;
    }
    return $m;
}

# ------------- main per-chrom processor -------------
sub process_chrom {
    my ($chr, $recs) = @_;
    return if !@$recs;

    # Events & per-edge maxima (for percl/percr)
    my (%starts_at, %ends_at);
    # Junction arrays to keep output order
    my (@S, @E, @C, @Z, @NM); # start, end, count, strand_idx, name

    # Gather events
    for my $r (@$recs) {
        my ($s, $e, $name, $cnt, $strand) = @$r; # start, end, name, count, strand_char
        my $z = strand_idx($strand);
        push @S, $s; push @E, $e; push @C, $cnt; push @Z, $z; push @NM, $name;

        push @{ $starts_at{$s} }, [$z, $cnt];
        push @{ $ends_at{$e}   }, [$z, $cnt];
    }

    # All unique boundary coordinates (union of starts/ends), sorted
    my %keys = map { $_ => 1 } (keys %starts_at, keys %ends_at);
    my @coords = sort { $a <=> $b } keys %keys;
    return if @coords < 2; # nothing to do

    # Map coord -> index
    my %idx_of;
    for (my $i=0; $i<@coords; $i++) { $idx_of{$coords[$i]} = $i }

    # Precompute per-start / per-end maxima (for percl / percr)
    my %start_max; # start_max{x}[z] = max count on strand z starting at x
    my %end_max;   # end_max{x}[z]   = max count on strand z ending at x
    for my $x (keys %starts_at) {
        my @m = (0,0,0);
        for my $ev (@{ $starts_at{$x} }) {
            my ($z,$cnt)=@$ev; $m[$z] = $cnt if $cnt > $m[$z];
        }
        $start_max{$x} = \@m;
    }
    for my $x (keys %ends_at) {
        my @m = (0,0,0);
        for my $ev (@{ $ends_at{$x} }) {
            my ($z,$cnt)=@$ev; $m[$z] = $cnt if $cnt > $m[$z];
        }
        $end_max{$x} = \@m;
    }

    # Sweep to build segment-wise active maxima per strand
    # Active sets per strand as freq maps: count -> how many active with that count
    my @freq = ({}, {}, {});
    my @curmax = (0,0,0);
    my @dirty  = (0,0,0);

    # Process first coordinate events BEFORE first segment so intervals starting at coords[0] are active in [coords[0], coords[1])
    my $prev = $coords[0];
    if (exists $ends_at{$prev}) { # typically none
        for my $ev (@{ $ends_at{$prev} }) {
            my ($z,$cnt)=@$ev;
            my $h = $freq[$z];
            if (exists $h->{$cnt}) {
                $h->{$cnt}--;
                delete $h->{$cnt} if $h->{$cnt} <= 0;
            }
            $dirty[$z] = 1 if $curmax[$z] == $cnt && !exists $h->{$cnt};
        }
        for my $z (0..2) {
            if ($dirty[$z]) { $curmax[$z] = recompute_max_from_freq($freq[$z]); $dirty[$z]=0; }
        }
    }
    if (exists $starts_at{$prev}) {
        for my $ev (@{ $starts_at{$prev} }) {
            my ($z,$cnt)=@$ev;
            my $h = $freq[$z];
            $h->{$cnt} = 0 unless exists $h->{$cnt};
            $h->{$cnt}++;
            $curmax[$z] = $cnt if $cnt > $curmax[$z];
        }
    }

    my @seg_max_neg; my @seg_max_pos; my @seg_max_neu;
    for (my $i=1; $i<@coords; $i++) {
        my $x = $coords[$i];

        # Segment [prev, x): store current maxima
        push @seg_max_neg, $curmax[0];
        push @seg_max_pos, $curmax[1];
        push @seg_max_neu, $curmax[2];

        # Handle ends at x (drop first for [prev, x))
        if (exists $ends_at{$x}) {
            for my $ev (@{ $ends_at{$x} }) {
                my ($z,$cnt)=@$ev;
                my $h = $freq[$z];
                if (exists $h->{$cnt}) {
                    $h->{$cnt}--;
                    delete $h->{$cnt} if $h->{$cnt} <= 0;
                }
                $dirty[$z] = 1 if $curmax[$z] == $cnt && !exists $h->{$cnt};
            }
        }
        # Recompute max if dirty
        for my $z (0..2) {
            if ($dirty[$z]) { $curmax[$z] = recompute_max_from_freq($freq[$z]); $dirty[$z]=0; }
        }
        # Handle starts at x (be active for next segment)
        if (exists $starts_at{$x}) {
            for my $ev (@{ $starts_at{$x} }) {
                my ($z,$cnt)=@$ev;
                my $h = $freq[$z];
                $h->{$cnt} = 0 unless exists $h->{$cnt};
                $h->{$cnt}++;
                $curmax[$z] = $cnt if $cnt > $curmax[$z];
            }
        }
        $prev = $x;
    }

    # Build RMQ structures for each strand
    my ($st_neg, $lg_neg) = build_sparse(\@seg_max_neg);
    my ($st_pos, $lg_pos) = build_sparse(\@seg_max_pos);
    my ($st_neu, $lg_neu) = build_sparse(\@seg_max_neu);

    # Emit per junction in original order
    for (my $i=0; $i<@S; $i++) {
        my $s = $S[$i]; my $e = $E[$i]; my $cnt = $C[$i]; my $z = $Z[$i]; my $name = $NM[$i];
        my $ls = $idx_of{$s}; my $re = $idx_of{$e};

        # segments covered are [ls .. re-1]
        my ($same_max, $diff_max);
        if    ($z == 0) { # neg
            $same_max = rmq_max($st_neg, $lg_neg, $ls, $re-1);
            my $m1    = rmq_max($st_pos, $lg_pos, $ls, $re-1);
            my $m2    = rmq_max($st_neu, $lg_neu, $ls, $re-1);
            $diff_max = ($m1 > $m2) ? $m1 : $m2;
        } elsif ($z == 1) { # pos
            $same_max = rmq_max($st_pos, $lg_pos, $ls, $re-1);
            my $m1    = rmq_max($st_neg, $lg_neg, $ls, $re-1);
            my $m2    = rmq_max($st_neu, $lg_neu, $ls, $re-1);
            $diff_max = ($m1 > $m2) ? $m1 : $m2;
        } else {           # neutral
            $same_max = rmq_max($st_neu, $lg_neu, $ls, $re-1);
            my $m1    = rmq_max($st_neg, $lg_neg, $ls, $re-1);
            my $m2    = rmq_max($st_pos, $lg_pos, $ls, $re-1);
            $diff_max = ($m1 > $m2) ? $m1 : $m2;
        }

        $same_max = $cnt if $same_max < $cnt;          # include itself (safety)
        my $percs = $same_max > 0 ? $cnt / $same_max : 1.0;
        my $percd = $diff_max > 0 ? $cnt / $diff_max : 1.0;

        # percl/percr from per-start/per-end maxima on same strand
        my $cl = 1.0;
        if (exists $start_max{$s} && $start_max{$s}->[$z] > 0) {
            my $m = $start_max{$s}->[$z];
            $cl = $cnt / $m;
        }
        my $cr = 1.0;
        if (exists $end_max{$e} && $end_max{$e}->[$z] > 0) {
            my $m = $end_max{$e}->[$z];
            $cr = $cnt / $m;
        }

        printf "%s\t%d\t%d\t%s\t%d\t%s\t%.4f-%.4f-%.4f-%.4f\n",
            $chr, $s, $e, $name, $cnt, ($z==1?'+':$z==0?'-':'.'),
            $percs, $percd, $cl, $cr;
    }
}

# ------------- reading & driving -------------
my ($infile) = @ARGV;
die "Usage: $0 <sorted_junctions.bed>\n" unless defined $infile;

open my $F, '<', $infile or die "Cannot open $infile: $!";

my $printed_header = 0;
my $cur_chr = '';
my @buf = (); # per-chrom tuples: [start, end, name, count, strand_char]

while (my $line = <$F>) {
    chomp $line;
    if ($line =~ /^(track|\#)/) {
        unless ($printed_header) {
            print "# track name=junctions type=bedDetail description=\"percsame-percdifferent-percleft-percright\"\n";
            $printed_header = 1;
        }
        next;
    }
    my @a = split /\t/, $line;
    next unless @a >= 6;
    my ($chr, $s, $e, $name, $cnt, $strand) = @a[0..5];

    if ($cur_chr ne '' && $chr ne $cur_chr) {
        process_chrom($cur_chr, \@buf);
        @buf = ();
    }
    $cur_chr = $chr;
    push @buf, [$s+0, $e+0, $name, $cnt+0, $strand];
}
close $F;

# flush last chrom
if (@buf) {
    print "# track name=junctions type=bedDetail description=\"percsame-percdifferent-percleft-percright\"\n" unless $printed_header;
    process_chrom($cur_chr, \@buf);
}
