use strict;
use File::Find;
use Data::Dumper;
use Cwd;
use File::Path qw(make_path remove_tree);

my $num_args = $#ARGV + 1;
if ($num_args != 3) {
  print "\nUsage: Unified_genotyper.pl Directory old_or_new_genome Soybean_line_type \n";
  exit;
}
my $condition_type = $ARGV[2];
my $gentyp = $ARGV[1];
my $reference;
if($gentyp=~/old/){$reference = "/home/santosj/data/Gmax_v9.0/Gmax_v_0.9.fa";}
else{$reference = "/home/skhan/bio/Gmax_assembly/Gmax_275_v2.0.fa";}

print $reference,"\n";

my $pref = $ARGV[0];
chomp($condition_type);
chomp $pref;
chop $pref;
my $outdir = getcwd."/$condition_type"."_UG/";
print "UG directory ",$outdir,"\n";
my($vol) = (split/\//,$pref)[-1];
my $dir = getcwd."/".$vol."/";
print "Indel Directory ",$dir,"\n";
my @all_files = ();
my %hash;
find(\&print_name_if_dir, $dir);
my @files=sort grep/Indelrealign/,grep/\.bam$/,@all_files;
my $indellist = $dir.'_'.'inderealign.list';
open(my $fh,">$indellist");
print $fh join "\n",@files;
print $fh "\n";

my @file_types = qw(UnifiedGenotyper SNP_only INDEL_only SNP_filtered INDEL_filtered SNP_passed INDEL_passed);

my %hash_files = ();
%hash_files =  map{
my $output;
if($_=~/(UnifiedGenotyper|SNP_only|INDEL_only|SNP_filtered|INDEL_filtered|SNP_passed|INDEL_passed)/){
$output = $outdir.$_."_".$condition_type."_multi.vcf";
($_,$output);
}
else{
();
}
} @file_types;

print Dumper(\%hash_files),"\n";

print "java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reference -I $indellist -o $hash_files{'UnifiedGenotyper'}\n";
mkdir($outdir, 0770) unless(-d $outdir);

system "java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reference -glm BOTH -I $indellist -o $hash_files{'UnifiedGenotyper'}";
system "java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T SelectVariants -R $reference -V  $hash_files{'UnifiedGenotyper'} -selectType SNP -o $hash_files{'SNP_only'}";
system "java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T SelectVariants -R $reference -V  $hash_files{'UnifiedGenotyper'} -selectType INDEL -o $hash_files{'INDEL_only'}";
system "java -Xmx10g -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V $hash_files{'SNP_only'} --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0\" --filterName \"my_snp_filter\" -o $hash_files{'SNP_filtered'}";
system "java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V $hash_files{'INDEL_only'} --filterExpression \"QD < 2.0 || FS > 200.0\" --filterName \"my_indel_filter\" -o $hash_files{'INDEL_filtered'}";

filtered_to_pass($hash_files{'INDEL_filtered'},$hash_files{'INDEL_passed'});

filtered_to_pass($hash_files{'SNP_filtered'},$hash_files{'SNP_passed'});

unlink "$indellist";

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
system "vcf-stats $hash_files{'SNP_passed'} > $dump";
system "rm $temp_pass $temp_header";
}