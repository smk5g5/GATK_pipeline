use strict;
use File::Find;
use Data::Dumper;
use Cwd;
use File::Path qw(make_path remove_tree);
use File::Basename;

my $num_args = $#ARGV + 1;
if ($num_args != 1) {
print "\nUsage: Haplotype_caller.pl genotype\n";
exit;
}
my $gtype = $ARGV[0];
chomp $gtype;
my $reference = "/home/skhan/bio/Gmax_assembly/Gmax_275_v2.0.fa";
my @all_files = ();

find(\&print_name_if_dir, ".");

# print join "\n",@all_files;

my($indel_realign)=sort grep/Indelrealign/,grep/\.bam$/,grep/$gtype/,@all_files;
print $indel_realign,"\n";
my $haplotype = $indel_realign;
$haplotype =~s/Indelrealign/HaplotypeCaller/g;
$haplotype =~s/bam$/vcf/g;
print "java -Xmx10g -XX:+UseSerialGC  -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -R $reference -I $indel_realign -o $haplotype\n";
system "java -Xmx10g -XX:+UseSerialGC  -Djava.io.tmpdir=/home/skhan/tmp -jar /home/skhan/bio/GATK_3.0/GenomeAnalysisTK.jar -T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -R $reference -I $indel_realign -o $haplotype";


sub print_name_if_dir
{
    push(@all_files,$File::Find::name)if -f && $File::Find::name;
}
