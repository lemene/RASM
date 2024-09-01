chromosomes=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
rep='rep3'

# ground truth SV callset
f=open('simulation_'+rep+'/translocation-random','r')

vcf=f.read().split('\n')[:-1]
f.close()
for i in range(len(vcf)):
	if 'chr23' in vcf[i]:
		vcf[i]=vcf[i].split('chr23')[0]+'chrX'+vcf[i].split('chr23')[1]

# TRA calls reported by SV caller
f=open('simulation_'+rep+'/debreak_2021/debreak-allsv-merged-final','r')
pac=f.read().split('\n')[:-1]
new=[]
for c in pac:
	if 'Tra' in c and c.split('\t')[0] in chromosomes and c.split('\t')[2] in chromosomes : 
		new+=[c]
pac=new

f.close()

detected=[]
reported=len(pac)
window=1000
for cc in vcf:
	c=cc.split('\t')
	for dd in pac:
		d=dd.split('\t')
		if c[0]==d[0] and c[2]==d[2] and abs(int(c[1])-int(d[1]))<=window and abs(int(c[3])-int(d[3]))<=window:
			detected+=[cc]
			break
		if c[0]==d[2] and c[2]==d[0] and abs(int(c[1])-int(d[3]))<=window and abs(int(c[3])-int(d[1]))<=window:
			detected+=[cc]
			break

pac=new
correct=[]
for cc in pac:
	c=cc.split('\t')
	for dd in vcf:
		d=dd.split('\t')
		if c[0]==d[0] and c[2]==d[2] and abs(int(c[1])-int(d[1]))<=window and abs(int(c[3])-int(d[3]))<=window:
			correct+=[cc]
			break
		if c[0]==d[2] and c[2]==d[0] and abs(int(c[1])-int(d[3]))<=window and abs(int(c[3])-int(d[1]))<=window:
			correct+=[cc]
			break
print(str(len(detected))+'\t'+str(len(vcf))+'\t'+str(100*float(len(detected))/len(vcf))+'\n'+str(len(detected))+'\t'+str(reported)+'\t'+str(100*float(len(detected))/reported)+'\n')
try:
	print (2*(float(len(detected))/len(vcf))*(float(len(detected))/reported)/(float(len(detected))/len(vcf)+float(len(detected))/reported))
except:
	pass
