indir=/scratch3/users/nanje/hlatyping/h3a_hlatypes/results/optitype

for file in $(ls ${indir})
do 
	new=${file%%.*}
	#echo ${new}
	mv ${indir}/${file}/${file}_result.tsv ${indir}/${file}/${new}_result.tsv
done
