for a in test-*
do
pwd
echo $a
cd $a
mkdir copy
cd copy
cp ../*.out .
cp ../*.aut
rename 's/\.out$/.ref/' *.out
rename 's/\.out$/aut.ref/' *.aut
cp *.ref ..
rm *.ref
cd ..
rmdir copy
cd ..
done