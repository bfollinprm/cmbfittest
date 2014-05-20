#!/usr/bin/perl
#use Cwd;

#Use current directory as root
#$cosmomc = cwd;
$cosmomc = pwd;

$name =  $ARGV[0];
$numnodes = $ARGV[1];
$walltime = $ARGV[2];

#$dir = "/lustre/bbenson/params_042811_".$name."/";
#$ini = "/home/bbenson/paramfits/cosmomc.b11f/".$name."/params.clusters.hst.wcdm.ini";
$ini = "/home/bbenson/paramfits/cosmomc.lps12_wmap9/params.".$name.".ini";

#$dir = "/lustre/bbenson/chains/asz8_800_070511/";
$dir = "/lustre/bbenson/chains/lps12_wmap9_020213/";
#$dir = "/lustre/bbenson/chains/r11_asz8_800_070711/";

#if ($ini !~ m/\.ini/) {$ini= "$ini.ini"}

#system("mkdir", $dir);
system("cp", $ini, $dir);

$path = $cosmomc;
system("rm ./script_MPI");
open(Fout,">./script_MPI");
print Fout <<EMP;
#PBS -S /bin/bash
#PBS -N $name
#PBS -j oe
#PBS -o /lustre/bbenson/chains/logs/$name.out
#PBS -m n
#PBS -l nodes=$numnodes:ppn=8,walltime=$walltime,vmem=1300mb
#PBS -q normal
cd \${PBS_O_WORKDIR}

# print identifying info for this job
echo "Job \${PBS_JOBNAME} submitted from \${PBS_O_HOST} started `date` jobid \${PBS_JOBID}"

#export OMP_NUM_THREADS=8

application=\${PBS_O_WORKDIR}/cosmomc
/usr/local/mpich2-intel/bin/mpiexec -v -pernode \$application $ini

EMP
close(Fout);

#chdir("./scripts");
@args=("qsub","./script_MPI");
system(@args);
#chdir("../");
