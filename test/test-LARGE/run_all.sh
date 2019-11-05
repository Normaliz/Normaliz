for a in $1*.in
do
    time ../../source/normaliz -c  -x=32 $a
    sleep 10
done
