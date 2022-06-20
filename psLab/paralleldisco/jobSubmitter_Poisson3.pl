#!/usr/bin/perl
use strict;

# Includes poisson weighting 

# NOTE:
# Numerical Recipes ran1 routine is initialized 
# with a NEGATIVE, NON-ZERO integer.


#my @nSigs = (0);
#my @nSigs = (10,20,30,40,50,60,70,80,90,100);
#my @nSigs = (120,140,160,180,200,220,240);
#my @nSigs = (140,160,180,200,220,240);
#my @nSigs = (300,350,400,450,500,550,600,650,700,750,800);
my @nSigs = (0,40,60,80,100,120,140,160,180,200,220,240,260,280,300);

my $nJobs   = 10; # per nSig
my $nTrials = 100; # per nSig, per job

#Get unique seed modifier, depending on time
my @timeData = localtime(time);
#my $seed = join('', @timeData); 
#print $seed;
my $seed = join('', -1*$timeData[0], $timeData[1], $timeData[2], $timeData[3]);print $seed."\n";

# USE DIFFERENT RANDOM NUMBER SEED FOR EACH 


# Name the search to identify scripts for setting up Llh
#my $searchName = "GP_FermiGalDiffuse";
#my $searchName = "fermibubble";
my $searchName = "UHECR";
#my $searchName = "GP_IngelmanThunman96";

print "SUBMITTING JOBS FOR: " . $searchName . "\n";

#my $labmaindir = $ENV{ 'LAB_MAIN_DIR' };
my $labmaindir = "/net/user/naoko/psLab";

## Testlab environment
my $envSh     = "$labmaindir/labcore/start_naoko_getos.sh";
my $mainDir   = "$labmaindir/paralleldisco";
print "Env: " . $envSh . "\n";
print "mainDir: "."$mainDir" . "\n";


my $scriptDir = "$mainDir/scripts/$searchName";
print $scriptDir . "\n";
mkdir("$scriptDir") || print $!;

my $arrayindex=0;
foreach(@nSigs)
{
  for (my $index=$arrayindex*$nJobs; $index<($arrayindex+1)*$nJobs; $index++)
  {
    my $number = sprintf("%08d",abs($seed-$index));
    my $scriptFile = sprintf("%s/script_%s.sh",$scriptDir,$number);


    my $nSigMin = $_; # $_ is one entry 
    my $nSigMax = $_;


    if (-e $scriptFile) {
	print "Warning: $number exists:  Overwriting  $scriptFile\n"; 
    }

    open OUTPUT, ">$scriptFile";
    
    my $script = "$envSh  root -b -q '$mainDir/rootSubmitMacro_Poisson_UHECR.C($seed-$index,\"$searchName\",$nSigMin,$nSigMax,$nTrials)'\n";

    print OUTPUT $script;

    close OUTPUT;
    system "chmod 744 $scriptFile";
  }
  $arrayindex++;
}


print $nJobs*$arrayindex. " shell scripts written.\n";
print "Type 'yes' to submit jobs, or anything else to quit.\n";

my $choice = <STDIN>;


if ($choice eq "yes\n") {

  my $count = 0;
  my $arrayindex=0;
  foreach(@nSigs)
  {
    for (my $index=$arrayindex*$nJobs; $index<($arrayindex+1)*$nJobs; $index++)
    {
	++$count;
    my $number = sprintf("%08d",abs($seed-$index));
	my $scriptFile = sprintf("%s/script_%s.sh",$scriptDir,$number);

	print "Job being submitted: $count\n";
	system "./scriptSubmit_npx2-uwa.sh $number $scriptFile";
	
    }

    print $count . " jobs submitted.\n";
    $arrayindex++;
  }
} else {
    print "No jobs submitted.\n";
}
