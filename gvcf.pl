use strict;
use File::Find;
use Data::Dumper;
use Cwd;
use File::Path qw(make_path remove_tree);

my $num_args = $#ARGV + 1;
if ($num_args != 1) {
  print "\nUsage: gvcf.pl Soybean_line_type \n";
  exit;
}

my $condition_type = $ARGV[0];
chomp($condition_type);

my $outdir = getcwd."/$condition_type"."_GVCF/";
print $outdir,"\n";
my @all_files = ();
my $reference = "/home/skhan/bio/Gmax_assembly/Gmax_275_v2.0.fa";
my %hash;
find(\&print_name_if_dir, ".");
my @files=sort grep/HaplotypeCaller/,grep/\.vcf$/,@all_files;

my @file_types = qw(GVCF SNP_only INDEL_only SNP_filtered INDEL_filtered SNP_passed INDEL_passed);

my %hash_files = ();
%hash_files =  map{
my $output;
if($_=~/(GVCF|SNP_only|INDEL_only|SNP_filtered|INDEL_filtered|SNP_passed|INDEL_passed)/){
$output = $outdir.$_."_".$condition_type.".vcf";
($_,$output);
}
else{
();
}
} @file_types;

print Dumper(\%hash_files),"\n";

mkdir($outdir, 0770) unless(-d $outdir);

my $variant;
map{$variant.=" --variant $_ "}@files;
print $variant,"\n";
print "java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $reference $variant -o $hash_files{'GVCF'}\n";
system "java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $reference $variant -o $hash_files{'GVCF'}";
system "java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T SelectVariants -R $reference -V  $hash_files{'GVCF'} -selectType SNP -o $hash_files{'SNP_only'}";
system "java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T SelectVariants -R $reference -V  $hash_files{'GVCF'} -selectType INDEL -o $hash_files{'INDEL_only'}";
system "java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V $hash_files{'SNP_only'} --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0\" --filterName \"my_snp_filter\" -o $hash_files{'SNP_filtered'}";
system "java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V $hash_files{'INDEL_only'} --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filterName \"my_indel_filter\" -o $hash_files{'INDEL_filtered'}";

filtered_to_pass($hash_files{'INDEL_filtered'},$hash_files{'INDEL_passed'});

filtered_to_pass($hash_files{'SNP_filtered'},$hash_files{'SNP_passed'});

sub print_name_if_dir
{
    push(@all_files,$File::Find::name)if -f && $File::Find::name;
}

sub filtered_to_pass{
my $vcf = shift;
my $outvcf = shift;
my $temp_pass = "$vcf"."temp";
my $temp_header = "$vcf"."temp_header";
system "grep PASS $vcf >$temp_pass";
print "grep \# $vcf > $temp_header\n";
my $regex = '\\#';
`grep $regex $vcf > $temp_header`;
`cat $temp_header $temp_pass >$outvcf`;
my $dump = $outvcf.".dump";
system "vcf-stats $outvcf > $dump";
system "rm $temp_pass $temp_header";
}