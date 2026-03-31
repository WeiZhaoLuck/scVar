#!/usr/bin/perl -w
open(IN,$ARGV[0]);
#use Getopt::Long;
#Getopt::Long::GetOptions( 'file=s' => \$path, 't=i' => \$t)
#open(DATA, $path)
#foreach $line (<IN>) {
#    $test=chomp($line);
#    print $test."\t".'VAF';
#}
#close(IN);
while(<IN>){
    chomp($_);
    $raw=$_;
    if(/^[\d|\w]/){
      if(/:(\d+\,\d*):/){
        my @result=split /\,/,$1;if(($result[0]+$result[1])>0){
            my $vaf=$result[1]*100/($result[0]+$result[1]);if($vaf<=$ARGV[1]){
                print $raw."\t".$vaf."\n";
        }
      }
    }
    }else{
        print $_."\t"."VAF"."\n";
    }
}
