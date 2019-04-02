# Instructions for Compiling C-language Grid Tools

The C-language tools can be compiled via the script `wrap_compile.pl`. This script uses the `gcc` compiler, which must be installed on your system. `wrap_compile.pl` is a Perl script that requires Perl 5 or later. Before running, make sure the script is executable by running:

   `chmod +x wrap_compile.pl`

Also, make sure to add `$PROJECT/tools/bin` to your `$PATH` variable.

Once these are complete, run the script. Usage:

   `wrap_compile.pl $INDIR $OUTDIR`

where:

 - `$INDIR` = directory containing the `.c` files, i.e. `$PROJECT/tools/src` where `$PROJECT` is the path to the directory containing your clone of this GitHub repo.
 - `$OUTDIR` = directory where the output executable files should be stored, i.e. `$PROJECT/tools/bin`.
