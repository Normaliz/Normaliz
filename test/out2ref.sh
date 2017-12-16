for a in test-*
do
pwd
echo $a
cd $a
mkdir copy
cd copy
cp ../*.out .
rename 's/\.out$/.ref/' *.out
cp *.ref ..
rm *.ref
cd ..
rmdir copy
cd ..
done