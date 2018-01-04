#!/usr/bin/perl

$in_file  = shift @ARGV;
$fwd_file = shift @ARGV;
$rev_file = shift @ARGV;

@fwd_buf = ();
@rev_buf = ();

$pair = 'f';
$line_cnt = -1;
open (IN, $in_file);
while (<IN>) {
    chomp;
    ++$line_cnt;
    if ($line_cnt > 0 && ($line_cnt % 4) == 0) {
	$pair = ($pair eq 'f') ? 'r' : 'f';
    }
    if ($pair eq 'f') {
	push (@fwd_buf, $_);
    } else {
	push (@rev_buf, $_);
    }
}
close (IN);

open (FWD, '>'.$fwd_file);
print FWD join("\n", @fwd_buf)."\n";
close (FWD);

open (REV, '>'.$rev_file);
print REV join("\n", @rev_buf)."\n";
close (REV);
