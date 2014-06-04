use lib qw(/home/skhan/local_lib/);
use Data::Dumper;
use File::Slurp;
use List::MoreUtils qw{uniq};
use Parallel::ForkManager;
 use Cwd;
my $num_args = $#ARGV + 1;
if ($num_args != 3) {
  print "\nUsage: newgatk_bsub.pl directory_containing_fastq_files old_or_new_genome Soybean_line_type\n";
  # exp_file
  exit;
}
my $dir = $ARGV[0];
my $genome_type =  $ARGV[1];
my $soybean_type = $ARGV[2];
chomp $soybean_type;
my $old_dir = getcwd;
chdir($dir);

my @HN = ();
my @files = grep/\.(fq|fastq)/,glob("*.*");
my @gt = map { (split/_/,$_)[0];} @files;

@HN = uniq(@gt);
print join "\n",@HN,"\n";
chdir($old_dir);
my @hap_jobs = ();
my @job_id;
foreach my $gtype(@HN){
print "bsub -R \"rusage[mem=15000] span[hosts=1]\" -J $gtype -o $gtype.o%J -e $gtype.e%J -q normal -n 1 \"perl gatk_new.pl $gtype $genome_type $dir\"\n\n";
my $out = `bsub -R \"rusage[mem=15000] span[hosts=1]\" -J $gtype -o $gtype.o%J -e $gtype.e%J -q normal -n 1 \"perl gatk_new.pl $gtype $genome_type $dir\"`;
my $id;
# Job <136145> is submitted to default queue <normal>.
if($out=~/Job\s+\<(\d{1,})\>\s+.*/g){$id = $1;}
print "jobid is ",$id,"\n";
push(@job_id,$id);
my $hap_job = "Haplotypecaller_".$gtype;
my $condition = "\"ended($id)\"";
my $hapout = `bsub -w $condition -R \"rusage[mem=15000] span[hosts=1]\" -J $hap_job -o $hap_job.o%J -e $hap_job.e%J -q normal -n 1 \"perl Haplotype_caller.pl $genome_type $gtype\"`;
my $hap_id;
if($hapout=~/Job\s+\<(\d{1,})\>\s+.*/g){$hap_id = $1;}
print "jobid is ",$hap_id,"\n";
push(@hap,$hap_id);
}
my $hap_cond = create_dependency(\@hap);
system  "bsub -w $hap_cond -R \"rusage[mem=15000] span[hosts=1]\" -J gvcf_$soybean_type -o gvcf_$soybean_type.o%J -e gvcf_$soybean_type.e%J -q normal -n 1 \"perl gvcf.pl $genome_type $soybean_type $dir\"";

# my $UG_cond = create_dependency(\@job_id);
# system  "bsub -w $UG_cond -R \"rusage[mem=15000] span[hosts=1]\" -J UG_$soybean_type -o UG_$soybean_type.o%J -e UG_$soybean_type.e%J -q normal -n 1 \"perl Unified_genotyper.pl $dir $genome_type $soybean_type\"";


sub create_dependency{
my $arr = shift;
my $counter = 0;
my $condition = '"';
my $and = ' && ';
map{
$counter++;
if($counter < scalar @$arr){
$condition.="ended($_)".$and;
}
else {$condition.="ended($_)".'"';}
}@$arr;
return($condition);
}