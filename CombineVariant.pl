#!/usr/bin/perl -w
use strict;

my ($samplelist,$indir,$outdir,$maf,$stat,$suffix)=@ARGV;

# Combine variants obtained from WGS and WES into one MAF

$suffix ||="tsv.gz";

`mkdir -p $outdir` unless (-e $outdir);

if($maf=~/\.gz$/){open OT,"|gzip >$maf" or die $!;}else{open OT,">$maf" or die $!;}
my $head="Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\tAAChange_refGene\tn_depth\tt_depth\tt_ref_count\tt_alt_count\tTumorVAF";
print OT "$head\n";
open STAT,">$stat" or die $!;
print STAT "CaseID\tType\tNumber\tRatio\n";
open INL,"$samplelist" or die $!;
while (<INL>){
	# Format (tab-delimited text file)
	# GlimsID Comparisons     CaseID  Gender
	# Bca_10  10N-vs-10T      CCGA-UBC-010    Male
	chomp;
	next if $.==1;
	my ($Comparisons,$CaseID)=(split /\t/)[1,2];
	$Comparisons=~s/\s+//g;
	$CaseID=~s/\s+//g;
	my $sample=$Comparisons;

	my (%line,%v,%pos,%detail);
	my @f=glob("$indir/*/$sample.$suffix");
	foreach my $f(@f){
		print "Reading: $f\n";
		my $type=(split /\//,$f)[-2];
		if($f=~/\.gz$/){open IN,"gunzip -cd <$f|" or die $!;}else{open IN,"<$f" or die $!;}
		while (<IN>){
			chomp;
			next if $.==1;
			my @t=split /\t/;
			my $line=join("\t",@t[0..@t-6]);
			$line{$line}++;
			$v{$line}{'n_depth'}+=$t[-5];
			$v{$line}{'t_depth'}+=$t[-4];
			$v{$line}{'t_ref_count'}+=$t[-3];
			$v{$line}{'t_alt_count'}+=$t[-2];	

			$detail{$line}{$type}="";
			my $pos="$t[4]:$t[5]";
			$pos{$pos}++;	
		}
		close IN;
	}
	foreach my $k(sort keys %pos){
		print "$sample\t$k\t$pos{$k}\n" if ($pos{$k} >2);
	}

	my $maf_sample="$outdir/$sample.tsv.gz";
	if($maf_sample=~/\.gz$/){open Sub,"|gzip >$maf_sample" or die $!;}else{open Sub,">$maf_sample" or die $!;}
	print Sub "$head\n";
	my ($both,$wgs,$wes)=("0","0","0");
	foreach my $k(sort keys %line){
		
		my $n_depth=$v{$k}{'n_depth'};
		my $t_depth=$v{$k}{'t_depth'};
		my $t_ref_count=$v{$k}{'t_ref_count'};
		my $t_alt_count=$v{$k}{'t_alt_count'};
		my $TumorVAF=$t_alt_count/$t_depth;

		print Sub "$k\t$n_depth\t$t_depth\t$t_ref_count\t$t_alt_count\t$TumorVAF\n";
		print OT "$k\t$n_depth\t$t_depth\t$t_ref_count\t$t_alt_count\t$TumorVAF\n";

		if (exists $detail{$k}{'WGS'} and exists $detail{$k}{'WES'}){
			$both++;
		}
		elsif (exists $detail{$k}{'WGS'}){
			$wgs++;
		}
		elsif (exists $detail{$k}{'WES'}){
			$wes++;
		}
	}
	close Sub;
	my $sum=$both+$wgs+$wes;
	my $ratio_wgs=$wgs/$sum;
	my $ratio_wes=$wes/$sum;
	my $ratio_both=1-$ratio_wgs-$ratio_wes;
	print STAT "$CaseID\tboth\t$both\t$ratio_both\n";
	print STAT "$CaseID\tWGS\t$wgs\t$ratio_wgs\n";
	print STAT "$CaseID\tWES\t$wes\t$ratio_wes\n";
}
close INL;
close OT;
close STAT;
exit;
