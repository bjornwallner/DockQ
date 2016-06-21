#!/usr/bin/perl 
use File::Temp qw/ tempfile /;


if(scalar(@ARGV)==0) {
    print STDERR "\nusage: ./fix_numbering.pl <model.pdb> <template.pdb> <read_seq_from_atom_in_residue (if def)>\n OUTFILE: model.fixed\n\n ";
    
    $chneedle=`which needle`;
    print $chneedle,"\n";
    if ($chneedle =~ m/needle/)
    {
	print "needle found to be installed in the path\n";
    }
    else
    {
	print "needle NOT found in the path\n";
	print "YOU NEED TO HAVE 'needle' INSTALLED AS PART OF THE EMBOSS PACKAGE: http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html\n";
    }
    exit;
}






$pdb_model=$ARGV[0];
$pdb_template=$ARGV[1];
$use_CA=1; #to get sequence
if(defined($ARGV[2])) {
    $use_CA=0;
}
$a=0;
if($use_CA) {
    ($seq_model,$resnum_model)=aa321_resnum($pdb_model);
} else {
    ($seq_model,$resnum_model)=aa321_resnumANY($pdb_model);
}
$a=1 if(defined($resnum_model));

if($use_CA) {
    ($seq_template,$resnum_template)=aa321_resnum($pdb_template);
} else {
    ($seq_template,$resnum_template)=aa321_resnumANY($pdb_template);

}
($ali_model,$ali_template)=align($seq_model,$seq_template);
if(length($ali_model)==0) {

    exit;
}
@ali_model=split(//,$ali_model);
@ali_template=split(//,$ali_template);
@ali_resnum_model=(); 
my $pos=0;
my $insert_num=9000;
for(my $i=0;$i<scalar(@ali_template);$i++) {
    
    if($ali_model[$i] ne '-') { 
	if($ali_template[$i] ne '-') {
	    push(@ali_resnum_model,${$resnum_template}[$pos]);
	    #print "$ali_model[$i] $ali_template[$i] ${$resnum_template}[$pos]\n";
	} else {
	    push(@ali_resnum_model,$insert_num."X");
	    $insert_num++;
	    #print "$ali_model[$i] $ali_template[$i] 999X\n";
	}
    }
    $pos++ if($ali_template[$i] ne '-');
}
#exit;
#open(OUT,>"$pdb_model.num");
print  ">MODEL:\n$ali_model\n\n>NATIVE:\n$ali_template\n";
#print STDERR ">MODEL:\n$ali_model\n\n>NATIVE:\n$ali_template\n";
#exit;
my $old_resnum="whatever";
$pos=-1;
#ATOM   1226  N   GLY A 188A     20.672  19.160  17.606  1.00 26.27  
open(MODEL,$pdb_model);
open(OUT,">$pdb_model.fixed");
while(<MODEL>) {
 #   print "A: ";
  #  print;
    if(/^ATOM/) {
#	my $atomno=substr($_, 7, 4);
#	my $atomtypeq=substr($_, 13, 3);
	my $resnum=substr($_,22,5);
	$resnum=~s/\s+//g;
	#print "$resnum $old_resnum $atomtype\n";
	if($old_resnum ne $resnum)
	{
	    $pos++;
#	    print "POS $ali_resnum_model[$pos] $_";
	}
	$old_resnum=$resnum;
	#print $pos."\n";
	if($ali_resnum_model[$pos]=~/X/) {
	    substr($_,22,5)=sprintf("%-5s",$ali_resnum_model[$pos]);
	} else {
	    substr($_,22,5)=sprintf("%5s",$ali_resnum_model[$pos]);
	}
	
    }
#    print "B: ";
    print OUT;
}



sub align   # Takes two strings removes all dashes and returns the alignment. 
{
    my ($seq1,$seq2)=@_;
    my $needle_linux="/usr/bin/needle";
    my $needle_mac="/opt/local/bin/needle";
    my $osname = $^O;
    my $input1=$seq1;
    my $input2=$seq2;
    $seq1=~s/-//g;
    $seq2=~s/-//g;
    #$seq1=remove_dashes($seq1);
    #$seq2=remove_dashes($seq2);
    $seq1=~s/\n//g;
    $seq2=~s/\n//g;
    $seq1=~s/\s+//g;
    $seq2=~s/\s+//g;
    
    my ($fh1,$file1)=tempfile("/tmp/seq.XXXXXXXX");
    my ($fh2,$file2)=tempfile("/tmp/seq.XXXXXXXX");
    my ($fh3,$file3)=tempfile("/tmp/ali.XXXXXXXX");
    close($fh3);
    print $fh1 ">seq1\n$seq1\n";
    close($fh1);
    print $fh2 ">seq2\n$seq2\n";
    close($fh2);

    if($osname eq "linux" && -e $needle_linux) {
	$needle=$needle_linux;
    }
    if($osname eq "darwin" && -e $needle_mac) {
	$needle=$needle_mac;
    }

    #print "needle -aseq $file1 -bseq $file2 -gapopen 10 -gapextend 0.5 -outfile $file3\n";
    `needle -aseq $file1 -bseq $file2 -gapopen 1 -gapextend 0.5 -outfile $file3 > /dev/null 2>&1`;
    #print $file3."\n";
    ($ali_return1,$ali_return2)=parse_needle_output($file3);
    `rm $file1 $file2 $file3`;
   
   
    return ($ali_return1,$ali_return2);
    
}

sub parse_needle_output
{
    my $file=shift;
    my $seq1="";
    my $seq2="";
    my $header1="";
    my $header2="";
    
    my $isFirst=1;
    open(FILE,$file);
    while(<FILE>){
	next if (/^#/);
        if(/^(.{13})(.{6}\d)\ ([\w\-]+)/){
	  #  print "header:$1, seqnro:$2, seq:$3|\n";
	    my $header = $1;
	    my $seq = $3;
	    if ($isFirst){
		$seq1.=$seq;
		$header1 = $header;
		$isFirst=0;
	    }
	    else {
		$seq2.=$seq;
		$header2 = $header;
		$isFirst=1;
	    }
        }
	
    }
    close(FILE);
    if(length($seq1) == 0) {

	print STDERR "needle from the EMBOSS package needs to be installed.\n";
    }
    return($seq1,$seq2);
}
sub aa321_resnum
{
    my $file=shift;
    my $seq="";
    my $old_resnum="whatever";
    my @resnum=();
    open(PDB,"$file");
    while(<PDB>)
    {
	if(/^ATOM\s/)
	{
	    my $atomno=substr($_, 7, 4);
	    my $atomtype=substr($_, 12, 4);
	    my $resnum=substr($_,22,5);
#	    $resnum=~s/\s+//g;
	    #print "$resnum $old_resnum $atomtype\n";
	    if($atomtype=~/CA/ && $old_resnum ne $resnum)
	    {
		$res=substr($_,17, 3);
		$seq.=aa321($res);
	#	print $table{$res};
		push(@resnum,$resnum);
		$old_resnum=$resnum;
	    }
	}
	
	last if(/^ENDMDL/);
    }
    close(PDB);
    if(scalar(@resnum)==0) { #read sequence
	$seq=`grep -v '>' $file`;
	$seq=~s/\n//g;
	my @seq=split(//,$seq);
	for(my $i=0;$i<scalar(@seq);$i++) {
	    push(@resnum,$i+1);
	}
    }

   # print "\n";
    return ($seq,\@resnum);
}

sub aa321_resnumANY
{
    my $file=shift;
    my $seq="";
    my $old_resnum="whatever";
    my @resnum=();
    open(PDB,"$file");
    while(<PDB>)
    {
	if(/^ATOM/)
	{
	    my $atomno=substr($_, 7, 4);
	    my $atomtype=substr($_, 12, 4);
	    my $resnum=substr($_,22,5);
#	    $resnum=~s/\s+//g;
	    #print "$resnum $old_resnum $atomtype\n";
	    if($old_resnum ne $resnum)
	    {
		$res=substr($_,17, 3);
		$seq.=aa321($res);
	#	print $table{$res};
		push(@resnum,$resnum);
		$old_resnum=$resnum;
	    }
	}
	
	last if(/^ENDMDL/);
    }
    close(PDB);
   # print "\n";
    return ($seq,\@resnum);
}
sub aa321
{
    my $aa=shift;
    my %aa321=('ALA', 'A',
	       'ARG', 'R',
	       'ASN', 'N',
	       'ASP', 'D',
	       'CYS', 'C',
	       'GLN', 'Q',
	       'GLU', 'E',
	       'GLY', 'G',
	       'HIS', 'H',
	       'ILE', 'I',
	       'LEU', 'L',
	       'LYS', 'K',
	       'MET', 'M',
	       'PHE', 'F',
	       'PRO', 'P',
	       'SER', 'S',
	       'THR', 'T',
	       'TRP', 'W',
	       'TYR', 'Y',
	       'VAL', 'V',
	       'ASX', 'B',
	       'GLX', 'Z',
	       'XXX', 'A',
	       'MSE', 'M',
	       'FME', 'M',
	       'PCA', 'E',
	       '5HP', 'E',
	       'SAC', 'S',
	       'CCS', 'C');
    my $aa1=0;
    $aa1=$aa321{$aa} if(defined($aa321{$aa}));
    return($aa1);

}
