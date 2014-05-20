# before running this script you should:
# mv /home/paramfits/wmap7 /lustre/rkeisler/wmap7

#chmod u+x runMPI.pl
#chmod u+x go_button

rm -r data
rm -r wmap9

ln -s ../wmap9 wmap9
ln -s ../data/ data
ln -s ../ptsrc ptsrc
ln -s ../spt_lowell spt_lowell

#cd ../wmap7
#ln -s ../wmap7/data data

cd ../cosmomc.lps12

rm camb/Makefile
rm source/Makefile
rm wmap9/Makefile

cd camb
ln -s Makefile_spt Makefile
cd ../source
ln -s Makefile_spt Makefile
cd ../wmap9
ln -s Makefile_spt Makefile


