
#Notes
# These are the notes to make cosmomc.lps12_wmap9

#in WMAP9 do, I did this in the paramfits/wmap9 directory (BB, 02-Feb-2013)
The recommended configuration is to turn off all but  "WMAPFLAGS  =
-DOPTIMIZE" flag. Then, in "WMAP_9yr_options.F90", make the following changes:

    Logical :: use_gibbs = .False.
->
    Logical :: use_gibbs = .True.

    Integer :: lowl_max = 30 ! use low l TT code 2<l<lowl_max
->
    Integer :: lowl_max = 32 ! use low l TT code 2<l<lowl_max

#Now do this:

#run this program
./sptafy.sh

#now do these steps to make
cd ../wmap9
make clean
make

cd ../cosmomc.lps12_wmap9/source
make clean
make

#cosmomc should be made now

