#! /bin/bash

for sample in `cat Gambianids.txt`
do
	ln -s /cbio/datasets/human/1kg/b37/low-coverage/fastq/$sample /users/nanje/link_Gwddata
done

