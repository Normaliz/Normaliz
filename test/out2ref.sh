for a in test-*
do
pwd
echo $a
cd $a
mkdir -p copy
cd copy
cp ../*.out .
cp ../*.aut .
cp ../*.ind .
cp ../*.cst .
cp ../*.tri .
cp ../*.fac .
cp ../*.gen .
cp ../*.grb .
cp ../*.inc .
cp ../*.mrk .
rename 's/\.out$/.ref/' *.out
rename 's/\.aut$/.aut.ref/' *.aut
rename 's/\.ind$/.ind.ref/' *.ind
rename 's/\.cst$/.cst.ref/' *.cst
rename 's/\.tri$/.tri.ref/' *.tri
rename 's/\.fac$/.fac.ref/' *.fac
rename 's/\.gen$/.gen.ref/' *.gen
rename 's/\.grb$/.grb.ref/' *.grb
rename 's/\.inc$/.inc.ref/' *.inc
rename 's/\.mrk$/.mrk.ref/' *.mrk
cp *.ref ..
rm *.ref
cd ..
rmdir copy
cd ..
done
