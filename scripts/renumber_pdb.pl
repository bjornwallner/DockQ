#!/usr/bin/perl -w 
#require 'bjornlib.pl';

if(scalar(@ARGV)==0) {
    print STDERR "\nusage: ./renumber_pdb.pl <PDB>\n OUTFILE: PDB.renum\n\n";
    exit;

}

$diff=1;
if(defined($ARGV[1])) {
    $diff=$ARGV[1];
}
#$file=number_pdb($ARGV[0],205+20);
$file=number_pdb($ARGV[0],$diff);
$outfile=$ARGV[0].'.renum';
open (OUT,">$outfile");
print OUT $file;


sub number_pdb #Takes a pdbfile and a number to start with, and starts renumbering the file from there
{
    my($file,$number)=@_;
    my $in_number=$number;
   # $number--;
    my $oldresno="";
    my $old_chain=" ";
    my $temp="";
    my $new_file="";
    my $atomcount=1;
    my $c="A";
#    my $c=" ";
    open(PDB,"$file");
    while(<PDB>)
    {
	my $line=$_;
	if($line=~/^ATOM/)
	{
	    my $chain=substr($line,21,1);
	    my $resno=substr($line, 22, 4);
	    if($oldresno ne $resno)
	    {
		$number++;
		if($old_chain ne $chain) { 
		    $number=$in_number;
		    $new_file.="TER\n" if(length($new_file)>0);
		}
		$temp=sprintf("%4d",$number);
	    }
	    substr($line,6,5)=sprintf("%5d",$atomcount);
	    $atomcount++;
	    substr($line, 22, 4)=$temp;

	    substr($line, 26, 1)=" ";
	    #substr($line, 21, 1)="$c";
	    $oldresno=$resno;
	    $old_chain=$chain;
	    $new_file.=$line;
	}
	if(/ENDMDL|^TER/) {
	 #   $number=$in_number;
	    #$c="B";


	}

    }
    return $new_file;
}
