#!/usr/bin/perl

$indir = shift;
$outdir = shift;

`mkdir -p $outdir`;

opendir (DIR,$indir) or die "$0: ERROR: cannot open dir $indir for reading\n";
@files = grep /\.c$/, readdir(DIR);
closedir(DIR);

foreach $file (sort(@files)) {
  $outfile = $file;
  $outfile =~ s/\.c$//;
  $cmd = "gcc $indir/$file -lm -o $outdir/$outfile";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
}
