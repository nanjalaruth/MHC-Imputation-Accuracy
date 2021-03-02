#! /bin/bash

for sample in `cat kgids.txt`
do
	ln -s /cbio/datasets/human/1kg/b37/low-coverage/fastq/$sample /users/nanje/data/link_kg_all_data
done

