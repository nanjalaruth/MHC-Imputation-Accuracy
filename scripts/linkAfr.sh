#! /bin/bash

for sample in `cat Africanids.txt`
do
	ln -s /cbio/datasets/human/1kg/b37/low-coverage/fastq/$sample /users/nanje/link_Afrdata
done

