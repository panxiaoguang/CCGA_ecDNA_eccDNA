#!/usr/bin/perl -w
use strict;

my ($samplelist,$indir,$outdir,$maf,$suffix)=@ARGV;

# To convert Annovar output into Mutation Annotation Format (MAF)
# Creator: Jinrong Huang <huangjinrong@genomics.cn>
# Date: Aug 2022

# NCBI_Build: Default GRCh38

$suffix ||="hg38_multianno.txt.gz";

`mkdir -p $outdir` unless (-e $outdir);

my %Annotation =(
	'frameshift insertion' => 'Frame_Shift_Ins',
	'frameshift deletion' => 'Frame_Shift_Del',
	'nonframeshift insertion' => 'In_Frame_Ins',
	'nonframeshift deletion' => 'In_Frame_Del',
	'nonsynonymous SNV' => 'Missense_Mutation',
	'synonymous SNV' => 'Silent',
	'stopgain' => 'Nonsense_Mutation',
	'stoploss' => 'Nonstop_Mutation',
	'startloss' => 'Translation_Start_Site',
	'splicing' => 'Splice_Site');

my %unknown;
if($maf=~/\.gz$/){open OT,"|gzip >$maf" or die $!;}else{open OT,">$maf" or die $!;}
my $head="Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\tAAChange_refGene\tn_depth\tt_depth\tt_ref_count\tt_alt_count\tTumorVAF";
print OT "$head\n";
open INL,"$samplelist" or die $!;
while (<INL>){
	# Format (tab-delimited text file)
	# GlimsID Comparisons     CaseID  Gender
	# Bca_10  10N-vs-10T      CCGA-UBC-010    Male
	chomp;
	next if $.==1;
	my ($Comparisons,$CaseID,$Gender)=(split /\t/)[1,2,3];
	$Comparisons=~s/\s+//g;
	$CaseID=~s/\s+//g;
	$Gender=~s/\s+//g;

	my ($column_normal,$column_tumor);
	if ($Comparisons =~/^\d+N-vs-\d+T$/){$column_normal="Otherinfo13";$column_tumor="Otherinfo14";}
	elsif ($Comparisons =~/^\d+T-vs-\d+N$/){$column_normal="Otherinfo14";$column_tumor="Otherinfo13";}
	else {print "Pls check: $Comparisons\n";}

	print "$CaseID\t$Comparisons\t$column_normal\t$column_tumor\n";

	my $sample=$Comparisons;
	my %column;
	my $maf_sample="$outdir/$sample.tsv.gz";
	if($maf_sample=~/\.gz$/){open Sub,"|gzip >$maf_sample" or die $!;}else{open Sub,">$maf_sample" or die $!;}
	print Sub "$head\n";
	my $f="$indir/$sample/$sample.$suffix";
	if($f=~/\.gz$/){open IN,"gunzip -cd <$f|" or die $!;}else{open IN,"<$f" or die $!;}
	while (<IN>){
		chomp;
		my @t=split /\t/;
		if ($.==1){
			foreach my $i(0..@t-1){
				$column{$t[$i]}=$i;
			}
			next;
		}
		if ($Gender eq "Female"){
			if ($t[$column{'Chr'}] eq "chrY"){
				next;			
			}
		}
		
		next if $t[$column{'Func_refGene'}]=~/intergenic|upstream|downstream|ncRNA_|intronic|UTR3|UTR5/;
		next unless $t[$column{'Otherinfo10'}]=~/PASS/;

		my $Chromosome =$t[$column{'Chr'}];
		my $Start_Position =$t[$column{'Start'}];
		my $End_Position =$t[$column{'End'}];

		my $Hugo_Symbol =$t[$column{'Gene_refGene'}];
		my $Entrez_Gene_Id ="-";

		my $Center="BGI";
		my $NCBI_Build="GRCh38";

		my ($Tumor_Seq_Allele1,$Tumor_Seq_Allele2);
		my $Reference_Allele =$t[$column{'Ref'}];
		my $Alt =$t[$column{'Alt'}];

		my $Refs =$t[$column{'Otherinfo7'}];
		my $Alts =$t[$column{'Otherinfo8'}];	
		if ($Alts=~/,/){
			print "$CaseID\t$Refs\t$Alts\t$_\n";
			next;
		}
		else {
			$Tumor_Seq_Allele1=$Reference_Allele;
			$Tumor_Seq_Allele2=$Alt;
		}
		my $Strand="+";

		my $Func_refGene =$t[$column{'Func_refGene'}];
		my $ExonicFunc_refGene = $t[$column{'ExonicFunc_refGene'}];

		my $AAChange_refGene = $t[$column{'AAChange_refGene'}];


		my $Variant_Classification;
		if ($ExonicFunc_refGene=~/\./){
			if (exists $Annotation{$Func_refGene}){
				$Variant_Classification=$Annotation{$Func_refGene};
			}
			else {
				$unknown{$Func_refGene}="";
				$Variant_Classification=$Func_refGene;
			}
		}
		else {
			if (exists $Annotation{$ExonicFunc_refGene}){
				$Variant_Classification=$Annotation{$ExonicFunc_refGene};
			}
			else {
				$unknown{$ExonicFunc_refGene}="";
				$Variant_Classification=$ExonicFunc_refGene;
			}
		}

		my $Variant_Type;
		my $length_Reference_Allele=length($Reference_Allele);
		my $length_Tumor_Seq_Allele2=length($Tumor_Seq_Allele2);
		if ($Tumor_Seq_Allele2 eq "-"){$Variant_Type="DEL";}
		elsif ($Reference_Allele eq "-"){$Variant_Type="INS";}
		elsif ($length_Reference_Allele eq $length_Tumor_Seq_Allele2){
			if ($length_Reference_Allele eq "1"){$Variant_Type="SNP";}
			elsif ($length_Reference_Allele eq "2"){$Variant_Type="DNP";}
			elsif ($length_Reference_Allele eq "3"){$Variant_Type="TNP";}
			elsif ($length_Reference_Allele > 3){$Variant_Type="ONP";}
			else {
				print "$Reference_Allele\t$Tumor_Seq_Allele2\t$_\n";
			}
		}
		else {
			print "$Reference_Allele\t$Tumor_Seq_Allele2\t$_\n";
		}
		my $DP_normal=(split /:/,$t[$column{$column_normal}])[3];
		my $AD_tumor=(split /:/,$t[$column{$column_tumor}])[1];

		my ($t_ref_count,$t_alt_count)=split /,/,$AD_tumor;
		my $t_depth=$t_ref_count+$t_alt_count;
		my $n_depth=$DP_normal;
		my $TumorVAF=$t_alt_count/$t_depth;

		my $Tumor_Sample_Barcode = $CaseID."T";

		my $line="$Hugo_Symbol\t$Entrez_Gene_Id\t$Center\t$NCBI_Build\t$Chromosome\t$Start_Position\t$End_Position\t$Strand\t$Variant_Classification\t$Variant_Type\t$Reference_Allele\t$Tumor_Seq_Allele1\t$Tumor_Seq_Allele2\t$Tumor_Sample_Barcode\t$AAChange_refGene\t$n_depth\t$t_depth\t$t_ref_count\t$t_alt_count\t$TumorVAF";
		print OT "$line\n";
		print Sub "$line\n"
	}
	close IN;
	close Sub;
}
close INL;
close OT;

print "The unknown classification:\t";
foreach my $k(sort keys %unknown){
	print "$k; ";
}
print "\n";
exit;
