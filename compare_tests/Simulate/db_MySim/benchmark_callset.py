import sys

def test(size1,size2,svty1,svty2,dep):
	chromosomes=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
	# ground truth sv set
	f=open('confident-'+svty1,'r')
	vcf=f.read().split('\n')[:-1]
	f.close()

	# callset of SV callers
	f=open('simulation_'+dep+'/debreak-allsv-merged-final','r')
	pac=f.read().split('\n')[:-1]
	f.close()

	pac=[c for c in pac if svty2 in c]

	#sensitivity / number of detected vcf
	groups=[]
	for c in pac:
		if groups==[]:
			groups=[[c]]
		else:
			i=0
			for group in groups:
				chrom=group[0].split('\t')[0]
				if c.split('\t')[0]==chrom:
					group+=[c]
					i=1
					break
			if i==0:
				groups+=[[c]]
	detected=[]
	detected2=[]
	no=[]
	bpshift=[]
	correct=[]
	for c in vcf:
		iii=0
		chrom=c.split('\t')[0]
		s1=int(c.split('\t')[1])
		s2=int(c.split('\t')[2])+s1
		length=int(c.split('\t')[2])
		window=1000
		sizeratio=0.5
		for group in groups:
			ch=group[0].split('\t')[0]
			if chrom!=ch:
				continue
			for d in group:
				d1=int(d.split('\t')[1])	
				if abs(d1-s1)<=window and sizeratio<=(float(d.split('\t')[2])/length)<=1/sizeratio:
					detected+=[c]
					correct+=[d]
					bpshift+=[d1-s1]
					group.remove(d)
					if group==[]:
						groups.remove([])	
					
					break
	#precision / number of correct pac
	correct=[]
	groups=[]
	for c in vcf:
		if groups==[]:
			groups=[[c]]
		else:
			i=0
			for group in groups:
				chrom=group[0].split('\t')[0]
				if c.split('\t')[0]==chrom:
					group+=[c]
					i=1
					break
			if i==0:
				groups+=[[c]]
	for c in pac:
		iii=0
		chrom=c.split('\t')[0]
		s1=int(c.split('\t')[1])
		s2=int(c.split('\t')[2])+s1
		length=abs(int(c.split('\t')[2]))
		window=1000
		sizeratio=0.5
		for group in groups:
			ch=group[0].split('\t')[0]
			if chrom!=ch:
				continue
			for d in group:
				d1=int(d.split('\t')[1])
				if abs(s1-d1)<=window and sizeratio<=(float(d.split('\t')[2])/length)<=1/sizeratio:
					correct+=[c]
					group.remove(d)
					if group==[]:
						groups.remove([])
					break
	#done




	false=[c for c in pac if c not in correct]
	missed=[c for c in vcf if c not in detected]
	if len(vcf)==0:
		sen=0.0
	else:
		sen=round(float(len(detected))/len(vcf)*10000)/100.0
	if len(pac)==0:
		pre=0.0
	else:
		pre=round(float(len(correct))/len(pac)*10000)/100.0
	if sen==0 and pre==0:
		f1=0
	else:
		f1=2*float(len(detected))/len(vcf)*len(correct)/len(pac)/(float(len(detected))/len(vcf)+float(len(correct))/len(pac))
	
	print(str(size1)+'-'+str(size2)+'\t'+str(len(detected))+'\t'+str(len(vcf))+'\t'+str(sen)+'\t'+str(len(correct))+'\t'+str(len(pac))+'\t'+str(pre))
	try:
		print(2*sen*pre/(sen+pre)/100)
	except:
		pass
	f=open('bpshift_debreak_'+dep,'a')
	for c in bpshift:
		f.write(svty2+'\t'+str(c)+'\n')
	f.close()
	return True


test(50,40000000,'deletion','Deletion','rep1')
