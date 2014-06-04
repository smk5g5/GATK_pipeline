use lib qw(/home/skhan/local_lib/);
use Data::Dumper;
use File::Slurp;
use List::MoreUtils qw{uniq};
use Parallel::ForkManager;

my $num_args = $#ARGV + 1;
if ($num_args != 3) {
  print "\nUsage: gatk_new.pl genome_prefix old_or_new_genome directory_containing_fastq_files\n";
  # exp_file
  exit;
}
my $gen_prefix = $ARGV[0];
my $gentyp = $ARGV[1];
my $reference;
if($gentyp=~/old/){$reference = "/home/santosj/data/Gmax_v9.0/Gmax_v_0.9.fa";}
else{$reference = "/home/skhan/bio/Gmax_assembly/Gmax_275_v2.0.fa";}

print $reference,"\n";

my $dir = $ARGV[2];
chomp $dir;
chdir($dir);
my $outdir = "GATK_3.0_".$gen_prefix;
mkdir $outdir;
my @files = grep/\.(fq|fastq)/,glob("*.*");
# my @gt = map { (split/_/,$_)[0];} @files;
# my @HN = uniq(@gt);
print join("\n",@files),"\n";
my @file_types = ('Sam','Sorted_sam','MarkDup','AddOrRep','intervals','Indelrealign','UnifiedGenotyper','HaplotypeCaller');

# exit;
# my $pm = Parallel::ForkManager->new(5);

# foreach my $gen_prefix(@HN){
# $pm->start and next; # do the fork
# 'SNP_UG','Indel_UG','filtered_snp','filtered_indel'
my @fastq = sort grep/$gen_prefix/,@files;
print join("\t",@fastq),"\n";
my %hash_files = ();
%hash_files =  map{
my $output;
if($_=~/(Sam|MarkDup|Sorted_sam)/){
$output = "$outdir/".$_."_".$gen_prefix.".sam";
($_,$output);
}
elsif($_=~/(AddOrRep|Indelrealign|PrintReads|ReduceReads)/){
$output = "$outdir/".$_."_".$gen_prefix.".bam";
($_,$output);
}
elsif($_=~/(UnifiedGenotyper|HaplotypeCaller)/){
$output = "$outdir/".$_."_".$gen_prefix.".vcf";
($_,$output);
}
else{
$output = "$outdir/".$gen_prefix."_realigner.intervals";
('intervals',$output);
}
} @file_types;

print Dumper(\%hash_files),"\n";
# java -Xmx10g -XX:+UseSerialGC -jar
# -XX:+UseSerialGC
@commands_argv = (
"
bwa mem -M $reference  $fastq[0] $fastq[1] > $hash_files{'Sam'};
java -Xmx10g -XX:+UseSerialGC -jar /home/skhan/bio/picard-tools-1.107/SortSam.jar CREATE_INDEX=TRUE MAX_RECORDS_IN_RAM=5000000  TMP_DIR= /home/skhan/tmp/  I=$hash_files{'Sam'} O=$hash_files{'Sorted_sam'} SO=coordinate VALIDATION_STRINGENCY=LENIENT;
java -Xmx10g -XX:+UseSerialGC -jar /home/skhan/bio/picard-tools-1.107/MarkDuplicates.jar CREATE_INDEX=TRUE MAX_RECORDS_IN_RAM=5000000 TMP_DIR=/home/skhan/tmp/  I=$hash_files{'Sorted_sam'} O=$hash_files{'MarkDup'} METRICS_FILE=$outdir/Duplicates_$gen_prefix VALIDATION_STRINGENCY=LENIENT;
java -Xmx10g -XX:+UseSerialGC -jar /home/skhan/bio/picard-tools-1.107/AddOrReplaceReadGroups.jar  MAX_RECORDS_IN_RAM=5000000  TMP_DIR=/home/skhan/tmp/  I=$hash_files{'MarkDup'} O=$hash_files{'AddOrRep'} RGID=$gen_prefix LB=$gen_prefix PL=\"Illumina\" SM=$gen_prefix CN=BGI RGPU=$gen_prefix VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate CREATE_INDEX=TRUE;
"
);
system(@commands_argv) == 0
        or die "system @args failed: $?";

		print join("\n\n\n",@commands_argv),"\n";

delete_files($outdir,$gen_prefix,('Sam','MarkDup','Sorted_sam'));

# $gen_prefix

my $outfile = "$outdir/$gen_prefix.txt";

`java -Xmx10g -XX:+UseSerialGC -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -I $hash_files{'AddOrRep'} -o $hash_files{'intervals'} 2> $outfile`;

my $out = read_file("$outfile");
if(($out=~/extremely\s+high\s+quality\s+score/)||($out=~/wrong\s+encoding/)){
print "extremely high quality score detected\n";
my @gatk_sys =
(
"
java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -I $hash_files{'AddOrRep'} -o $hash_files{'intervals'} --fix_misencoded_quality_scores;
java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -I $hash_files{'AddOrRep'} -targetIntervals $hash_files{'intervals'} -o $hash_files{'Indelrealign'} --fix_misencoded_quality_scores;
"
);
system(@gatk_sys) == 0
        or die "system @args failed: $?";
}
else{
print "no fallacies detected\n";
my @gatk_sys =
(
"
java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -I $hash_files{'AddOrRep'} -targetIntervals $hash_files{'intervals'} -o $hash_files{'Indelrealign'};
"
);
system(@gatk_sys) == 0
        or die "system @args failed: $?";
}

##steps removed
# java -Xmx10g -XX:+UseSerialGC  -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reference -I $hash_files{'Indelrealign'} -glm SNP -o $hash_files{'SNP_UG'};
# java -Xmx10g -XX:+UseSerialGC  -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reference -I $hash_files{'Indelrealign'} -glm INDEL -o $hash_files{'Indel_UG'};


# system "java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V $hash_files{'SNP_UG'} --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0\" --filterName \"my_snp_filter\" -o $hash_files{'filtered_snp'}";
# system "java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V $hash_files{'Indel_UG'} --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filterName \"my_indel_filter\" -o $hash_files{'filtered_indel'}";
# my $temp_pass = "$hash_files{'filtered_snp'}"."temp";
# my $temp_header = "$hash_files{'filtered_snp'}"."temp_header";
# system "grep PASS $hash_files{'filtered_snp'}>$temp_pass";
# print "grep \# $hash_files{'filtered_snp'}> $temp_header\n";
# my $regex = '\#';
# `grep $regex $hash_files{'filtered_snp'}> $temp_header`;
# `cat $temp_header $temp_pass >$hash_files{'passed_snp'}`;
# my $dump = $hash_files{'passed_snp'}.".dump";
# system "vcf-stats $hash_files{'passed_snp'} > $dump";
# system "rm $temp_pass $temp_header";
# delete_files($outdir,$gen_prefix,('PrintReads'));

# unlink("$outdir/$gen_prefix.txt");
# $pm->finish; # do the exit in the child process
		# exit;
# }
# $pm->wait_all_children;



sub delete_files{
my $dir_name = shift;
my $gtype = shift;
my(@df) = @_;
foreach my $type(@df){
if($type=~/(PrintReads)/){
my $file = "$dir_name/".$type."_".$gtype.'.bam';
unlink "$file" or warn "Could not unlink $file: $!";
my $file = "$dir_name/".$type."_".$gtype.'.bai';
unlink "$file" or warn "Could not unlink $file: $!";
}
elsif($type=~/(Sam|Sorted_sam|MarkDup)/){
my $file = "$dir_name/".$type."_".$gtype.'.sam';
unlink "$file" or warn "Could not unlink $file: $!";
}
else{next;}
}
}
