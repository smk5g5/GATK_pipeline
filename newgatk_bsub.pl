use lib qw(/home/skhan/local_lib/);
use Data::Dumper;
use File::Slurp;
use List::MoreUtils qw{uniq};
use Parallel::ForkManager;
 use Cwd;
my $num_args = $#ARGV + 1;
if ($num_args != 1) {
  print "\nUsage: newgatk_bsub.pl directory_containing_fastq_files \n";
  # exp_file
  exit;
}
my $dir = $ARGV[0];
my $old_dir = getcwd;
chdir($dir);

my @HN = ();
my $reference = "/home/skhan/bio/Gmax_assembly/Gmax_275_v2.0.fa";
my @files = grep/\.(fq|fastq)/,glob("*.*");
my @gt = map { (split/_/,$_)[0];} @files;

@HN = uniq(@gt);
print join "\n",@HN;
chdir($old_dir);

foreach my $gtype(@HN){
print "bsub -R \"rusage[mem=15000] span[hosts=1]\" -J $gtype -o $gtype.o%J -e $gtype.e%J -q normal -n 1 \"perl gatk_new.pl $gtype $dir\"\n\n";
system  "bsub -R \"rusage[mem=15000] span[hosts=1]\" -J $gtype -o $gtype.o%J -e $gtype.e%J -q normal -n 1 \"perl gatk_new.pl $gtype $dir\"";
}
