#!/bin/bash
#move everything that is inside a folder to the outside
subs=`ls $1`
mkdir fastg
for dir in $subs; do
	mv $1/$dir/* /users/nanje/fastg
done
