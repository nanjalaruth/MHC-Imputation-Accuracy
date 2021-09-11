#!/usr/bin/env python3

ann_file="${ann}"
r2_file="${rsquared}"

ann = { data.strip().split('\\t')[2]:data.strip().split('\\t') for data in open(ann_file).readlines() }
r2_data = { data.strip().split('\\t')[0]:data.strip().split('\\t')[1] for data in open(r2_file).readlines() }

out = open("${info_out}", 'w')
out.writelines('\\t'.join(["SNP", "REF(0)", "ALT(1)", "ALT_Frq", "MAF", "Rsq", "Genotyped"])+"\\n")
for data in ann:
	chr = ann[data][0]
	pos = ann[data][1]
	id = ann[data][2]
	ref = ann[data][3]
	alt = ann[data][4]
	maf = ann[data][5]
	r2 = '0'
	type = '-'
	if not id.startswith('SNP') and not id.startswith('AA') and not id.startswith('HLA'):
		type = "Genotyped"
		if id in r2_data:
                        r2 = r2_data[id]
	if id.startswith('SNP'):
		#print(data)
		#time.sleep(4)
		type = 'Imputed'
		if id in r2_data:
			r2 = r2_data[id]
	new_id = f"{chr}:{pos}:{ref}:{alt}"
	new_data = [new_id, ref, alt, maf, maf, r2, type]
	out.writelines('\\t'.join(new_data)+"\\n")
out.close()
