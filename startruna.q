#!/bin/bash
#

Npatterns=2

i=1
echo 'i='$i;
 
while [  $i -lt 351 ]; 
	do
	echo "$i";
        cd sima"$i";
        qsub rs-"$i".q
        cd ..
	let "i+=1";
done
