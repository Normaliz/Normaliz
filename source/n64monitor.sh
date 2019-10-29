#!/bin/sh

CMDLIST="normaliz,normaliz-27,normaliz-28,normaliz-29,normaliz-2101,normaliz-master,norm64-22,count,hilbert,zsolve,totalpyr,nopyr"

if [ $1 ]
then
	WAITTIME=$1
else
	WAITTIME=5s
fi



ps -C $CMDLIST  o pid,user,nice,pcpu,cputime,etime,pmem,vsize=VIRT,rss=RES #,thcount=THR

while [ 1 ]
do
	sleep $WAITTIME
	ps -C $CMDLIST  o pid=,user=,nice=,pcpu=,cputime=,etime=,pmem=,vsize=,rss= #,thcount=
done
