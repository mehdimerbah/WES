#!/usr/bin/env perl




my ($input) = @ARGV; 

open(FASTQ, '<', $input) or die "Can't open $input, $!\n"; # Open a file for reading.

my @Ids = ();
my $count = 0;
my $flag = 0;
my $curr_ID;
while (<FASTQ>){
	$flag = 0;
    chomp;
    if(/^@[A][0][0][7]/) {
    	$curr_ID = $input;
    	$flag=0}

    if(/[+]/) {$flag=1}
    else{

    if ($flag==1 && $count < 10)
    {
        if(length($input) < 36){
        		@Ids[$count] = $curr_ID;  
        		$count++;
        }
    }
}
print $count;
}


close (FASTQ);