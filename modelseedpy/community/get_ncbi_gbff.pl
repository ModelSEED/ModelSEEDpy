use strict;

while (<>){
    chomp ($_);
    next if ($_=~/^\s*$/);
    my $val = `grep $_  assembly_summary_refseq.txt |cut -f 20`;
    chomp ($val);
    my @p = split ("/", $val);
    my $n = $p[-1];
    my $url = "${val}/${n}_genomic.gbff.gz";
    my $fpath = "${n}_genomic.gbff.gz ";
    print "curl $url -o $fpath" . "\n";
}
