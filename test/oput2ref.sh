for a in test*
do
cd a
mkdir copy
cd copy
cp ../*.out
cd ..
done