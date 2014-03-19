use lib qw(/home/skhan/local_lib/);
use Data::Dumper;
use File::Slurp;
use List::MoreUtils qw{uniq};
use Parallel::ForkManager;

my $num_args = $#ARGV + 1;
if ($num_args != 2) {
  print "\nUsage: gatk_run.pl genome_prefix directory_containing_fastq_files\n";
  # exp_file
  exit;
}
my $gen_prefix = $ARGV[0];
my $dir = $ARGV[1];
chdir($dir);

my $reference = "/home/skhan/bio/Gmax_assembly/Gmax_275_v2.0.fa";
my @files = grep/\.(fq|fastq)/,glob("*.*");
# my @gt = map { (split/_/,$_)[0];} @files;
# my @HN = uniq(@gt);
print join("\n",@files),"\n";
my @file_types = ('Sam','Sorted_sam','MarkDup','AddOrRep','intervals','Indelrealign','bqsr','plot','recalib_bqsr','recalib_plot','csv_plot','csv_recalib','PrintReads','ReduceReads','UnifiedGenotyper','HaplotypeCaller');

# exit;
# my $pm = Parallel::ForkManager->new(5);

# foreach my $gen_prefix(@HN){
# $pm->start and next; # do the fork

my @fastq = sort grep/$gen_prefix/,@files;
print join("\t",@fastq),"\n";
my %hash_files = ();
%hash_files =  map{
my $output;
if($_=~/(Sam|MarkDup|Sorted_sam)/){
$output = $_."_".$gen_prefix.".sam";
($_,$output);
}
elsif($_=~/(AddOrRep|Indelrealign|PrintReads|ReduceReads)/){
$output = $_."_".$gen_prefix.".bam";
($_,$output);
}
elsif($_=~/(bqsr|recalib_bqsr)/){
$output = $_."_".$gen_prefix.".grp";
($_,$output);
}
elsif($_=~/(plot|recalib_plot)/){
$output = $_.$gen_prefix."_grp.pdf";
($_,$output);
}
elsif($_=~/(csv_plot|csv_recalib)/){
$output = $_.$gen_prefix."_grp.csv";
($_,$output);
}
elsif($_=~/(UnifiedGenotyper|HaplotypeCaller)/){
$output = $_."_".$gen_prefix.".vcf";
($_,$output);
}
else{
$output = $gen_prefix."_realigner.intervals";
($_,$output);
}
} @file_types;

print Dumper(\%hash_files),"\n";


@commands_argv = (
"
bwa mem -M $reference  $fastq[0] $fastq[1] > $hash_files{'Sam'};
java -Xmx10g -jar /home/skhan/bio/picard-tools-1.92/SortSam.jar CREATE_INDEX=TRUE MAX_RECORDS_IN_RAM=5000000  TMP_DIR= /home/skhan/tmp/  I=$hash_files{'Sam'} O=$hash_files{'Sorted_sam'} SO=coordinate VALIDATION_STRINGENCY=LENIENT;
java -Xmx10g -jar /home/skhan/bio/picard-tools-1.92/MarkDuplicates.jar CREATE_INDEX=TRUE MAX_RECORDS_IN_RAM=5000000 TMP_DIR=/home/skhan/tmp/  I=$hash_files{'Sorted_sam'} O=$hash_files{'MarkDup'} METRICS_FILE=Duplicates_$gen_prefix VALIDATION_STRINGENCY=LENIENT;
java -Xmx10g -jar /home/skhan/bio/picard-tools-1.92/AddOrReplaceReadGroups.jar  MAX_RECORDS_IN_RAM=5000000  TMP_DIR=/home/skhan/tmp/  I=$hash_files{'MarkDup'} O=$hash_files{'AddOrRep'} RGID=$gen_prefix LB=$gen_prefix PL=\"Illumina\" SM=$gen_prefix CN=BGI RGPU=$gen_prefix VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate CREATE_INDEX=TRUE;
"
);
system(@commands_argv) == 0
        or die "system @args failed: $?";

		print join("\n\n\n",@commands_argv),"\n";

delete_files($gen_prefix,('Sam','MarkDup','Sorted_sam'));

# $gen_prefix


`java -Xmx5g -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -I $hash_files{'AddOrRep'} -o $hash_files{'intervals'} 2> $gen_prefix.txt`;

my $out = read_file("$gen_prefix.txt");
if(($out=~/extremely\s+high\s+quality\s+score/)||($out=~/wrong\s+encoding/)){
# print "extremely high quality score detected\n";
# exit;
my @gatk_sys =
(
"
java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -I $hash_files{'AddOrRep'} -o $hash_files{'intervals'} --fix_misencoded_quality_scores;
java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -I $hash_files{'AddOrRep'} -targetIntervals $hash_files{'intervals'} -o $hash_files{'Indelrealign'} --fix_misencoded_quality_scores;
java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -I $hash_files{'Indelrealign'} -R $reference -run_without_dbsnp_potentially_ruining_quality -o $hash_files{'bqsr'} -plots $hash_files{'plot'} --intermediate_csv_file $hash_files{'csv_plot'};
java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -I $hash_files{'Indelrealign'} -R $reference -run_without_dbsnp_potentially_ruining_quality -BQSR $hash_files{'bqsr'} -o $hash_files{'recalib_bqsr'} -plots $hash_files{'recalib_plot'} --intermediate_csv_file $hash_files{'csv_recalib'};
java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T PrintReads -R $reference -I $hash_files{'Indelrealign'} -BQSR $hash_files{'bqsr'} -o $hash_files{'PrintReads'};
java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T ReduceReads -R $reference  -I $hash_files{'PrintReads'} -o $hash_files{'ReduceReads'};
java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reference -I $hash_files{'ReduceReads'} -o $hash_files{'UnifiedGenotyper'};
"
);
system(@gatk_sys) == 0
        or die "system @args failed: $?";
}
else{
# print "no fallacies detected\n";
# exit;
my @gatk_sys =
(
"
java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -I $hash_files{'AddOrRep'} -targetIntervals $hash_files{'intervals'} -o $hash_files{'Indelrealign'};
java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -I $hash_files{'Indelrealign'} -R $reference -run_without_dbsnp_potentially_ruining_quality -o $hash_files{'bqsr'} -plots $hash_files{'plot'} --intermediate_csv_file $hash_files{'csv_plot'};
java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -I $hash_files{'Indelrealign'} -R $reference -run_without_dbsnp_potentially_ruining_quality -BQSR $hash_files{'bqsr'} -o $hash_files{'recalib_bqsr'} -plots $hash_files{'recalib_plot'} --intermediate_csv_file $hash_files{'csv_recalib'};
java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T PrintReads -R $reference -I $hash_files{'Indelrealign'} -BQSR $hash_files{'bqsr'} -o $hash_files{'PrintReads'};
java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T ReduceReads -R $reference  -I $hash_files{'PrintReads'} -o $hash_files{'ReduceReads'};
java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reference -I $hash_files{'ReduceReads'} -o $hash_files{'UnifiedGenotyper'};
"
);
system(@gatk_sys) == 0
        or die "system @args failed: $?";
}


delete_files($gen_prefix,('PrintReads'));

unlink("$gen_prefix.txt");
# $pm->finish; # do the exit in the child process
		# exit;
# }
# $pm->wait_all_children;



sub delete_files{
my $gtype = shift;
my(@df) = @_;
foreach my $type(@df){
if($type=~/(PrintReads)/){
my $file = $type."_".$gtype.'.bam';
unlink "$file" or warn "Could not unlink $file: $!";
my $file = $type."_".$gtype.'.bai';
unlink "$file" or warn "Could not unlink $file: $!";
}
elsif($type=~/(Sam|Sorted_sam|MarkDup)/){
my $file = $type."_".$gtype.'.sam';
unlink "$file" or warn "Could not unlink $file: $!";
}
else{next;}
}
}
