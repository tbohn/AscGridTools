#!/usr/bin/perl

$infile = shift;
$miny = shift;
$maxy = shift;
$minx = shift;
$maxx = shift;
$op = shift; # lt, le, eq, ge, gt, ne
$opval = shift; # value
$include_nodata = shift; # 0 = do not operate on nodatas; 1 = paint them too
$newval = shift;

open(FILE, $infile) or die "$0: ERROR: cannot open file $infile for reading\n";
$row = 0;
foreach (<FILE>) {
  chomp;
  @fields = split /\s+/;
  if ($fields[0] =~ /nrows/i) {
    $nrows = $fields[1];
  }
  elsif ($fields[0] =~ /ncols/i) {
    $ncols = $fields[1];
  }
  elsif ($fields[0] =~ /xllcorner/i) {
    $xllcorner = $fields[1];
  }
  elsif ($fields[0] =~ /yllcorner/i) {
    $yllcorner = $fields[1];
  }
  elsif ($fields[0] =~ /cellsize/i) {
    $cellsize = $fields[1];
  }
  elsif ($fields[0] =~ /nodata/i) {
    $nodata = $fields[1];
  }
  else {
    $y = $yllcorner + ($nrows-1-$row+0.5)*$cellsize;
    for ($col=0; $col<$ncols; $col++) {
      $x = $xllcorner + ($col+0.5)*$cellsize;
      if ($y >= $miny && $y <= $maxy && $x >= $minx && $x <= $maxx) {
        if ($include_nodata || $fields[$col] != $nodata) {
          if (   ($op eq "lt" && $fields[$col] <  $opval)
              || ($op eq "le" && $fields[$col] <= $opval)
              || ($op eq "eq" && $fields[$col] == $opval)
              || ($op eq "ge" && $fields[$col] >= $opval)
              || ($op eq "gt" && $fields[$col] >  $opval)
              || ($op eq "ne" && $fields[$col] != $opval) ) {
            $fields[$col] = $newval;
          }
        }
      }
    }
    $row++;
  }
  $line = join " ", @fields;
  print "$line\n";
}
close(FILE);
