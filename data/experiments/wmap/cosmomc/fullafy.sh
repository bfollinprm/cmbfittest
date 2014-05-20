# before running this script you should:
# mv /home/paramfits/wmap7 /lustre/rkeisler/wmap7

#chmod u+x runMPI.pl
#chmod u+x go_button

rm -r ptsrc
rm -r data
rm -r wmap9
rm -r spt_lowell

ln -s /lustre/bbenson/paramfits/wmap9/ wmap9
ln -s /lustre/bbenson/paramfits/data/ data
ln -s /lustre/bbenson/paramfits/ptsrc/ ptsrc
ln -s /lustre/bbenson/paramfits/spt_lowell/ spt_lowell

cd ../cosmomc.lps12_wmap9

rm camb/Makefile
rm source/Makefile
rm wmap9/Makefile

cd camb
ln -s Makefile_fulla Makefile
cd ../source
ln -s Makefile_fulla Makefile
cd ../wmap9
ln -s Makefile_fulla Makefile



