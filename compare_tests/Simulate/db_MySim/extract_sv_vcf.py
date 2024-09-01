
def run_extract(filename):
	f=open(filename,'r')
	allvcf=f.read().split('\n')[:-1]
	f.close()
	allins=[];alldel=[];alldup=[];allinv=[];alltra=[];allcnv=[]
	chromosomes=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

	for c in allvcf:
		if c[0]=='#' or c.split('\t')[0] not in chromosomes or 'PASS' not in c:
			continue
		if 'SVTYPE=DEL;' in c :
			alldel+=[c]
		if ('SVTYPE=INS;' in c and 'SVLEN=999999999;' not in c) :
			allins+=[c]
		if 'SVTYPE=DUP;' in c :
			alldup+=[c]
		if 'SVTYPE=INV;' in c:
			allinv+=[c]
		if 'SVTYPE=BND;' in c:
			alltra+=[c]
		if 'SVTYPE=cnv;' in c:
			allcnv+=[c]
	ins=[]

	for c in allins:
		c=c.split('\t')
		if int(c[7].split('SVLEN=')[1].split(';')[0])>=45  :
			ins+=[c[0]+'\t'+c[1]+'\t'+c[7].split('SVLEN=')[1].split(';')[0]+'\t'+c[9].split(':')[0]]

	f=open(filename[:-4]+'_ins.info','w')
	for c in ins:
		f.write(c+'\n')
	f.close()

	dels=[]
	for c in alldel:
		c=c.split('\t')
		if int(c[7].split('SVLEN=-')[1].split(';')[0])>=45:
			dels+=[c[0]+'\t'+c[1]+'\t'+c[7].split('SVLEN=-')[1].split(';')[0]+'\t'+c[9].split(':')[0]]
	

	f=open(filename[:-4]+'_del.info','w')
	for c in dels:
		f.write(c+'\n')
	f.close()
	
	dup=[]
	for c in alldup:
		c=c.split('\t')
		dup+=[c[0]+'\t'+c[1]+'\t'+c[7].split('SVLEN=')[1].split(';')[0]+'\t'+c[9].split(':')[0]]

	f=open(filename[:-4]+'_dup.info','w')
	for c in dup:
		f.write(c+'\n')
	f.close()
	inv=[]
	for c in allinv:
		c=c.split('\t')
		inv+=[c[0]+'\t'+c[1]+'\t'+str(int(c[7].split('END=')[1].split(';')[0])-int(c[1]))+'\t'+c[9].split(':')[0]]

	f=open(filename[:-4]+'_inv.info','w')
	for c in inv:
		f.write(c+'\n')
	f.close()



	f=open(filename[:-4]+'_insdup.info','w')
	for c in dup+ins:
		f.write(c+'\n')

	f.close()



	tra=[]
	for c in alltra:
		chr1=c.split('\t')[0]
		chr2=c.split('\t')[2].split('-')[1].split(':')[0]
		if chr1 in chromosomes and chr2 in chromosomes and chr1!=chr2:
			pos1=c.split('\t')[1]
			pos2=c.split('\t')[2].split('-')[1].split(':')[1]
			tra+=[chr1+'\t'+pos1+'\t'+chr2+'\t'+pos2+'\t'+c[9].split(':')[0]]
	f=open(filename[:-4]+'_tra.info','w')
	for c in tra:
		f.write(c+'\n')
	f.close()
	cnv=[]
	for c in allcnv:
		c=c.split('\t')
		cnv+=[c[0]+'\t'+c[1]+'\t'+c[7].split('SVLEN=')[1].split(';')[0]+'\t'+c[9].split(':')[0]]
	
	f=open(filename[:-4]+'_cnv.info','w')
	for c in cnv:
		f.write(c+'\n')
	f.close()	
	

	return True


if __name__ == "__main__":
	run_extract('sim_rep1_pbsv.vcf')

