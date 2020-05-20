for a in test-*
do
pwd
echo $a
cd $a
mkdir -p copy
cd copy
cp ../*.out .
cp ../*.aut .
rename 's/\.out$/.ref/' *.out
rename 's/\.aut$/.aut.ref/' *.aut
cp *.ref ..
rm *.ref
cd ..
rmdir copy
cd ..
done
