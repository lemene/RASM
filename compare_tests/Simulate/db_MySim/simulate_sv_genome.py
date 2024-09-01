import random
import os

class Random():
	def __init__(self) -> None:
		pass

def random_length(svtype):
	'''随机生成SV的大小'''
	if svtype in ['deletion']:
		i=random.randint(121,8131)
	elif  svtype in ['insertion']:
		i=random.randint(227,8131)
	elif  svtype in ['duplication']:
		i=random.randint(227,3069)
	else:
		i=random.randint(121,3069)
	if 120<i<=146:
		length=random.randint(1000000,10000000)
	if 146<i<=226:
		length=random.randint(50000,1000000)
	if 226<i<=382:
		length=random.randint(10000,50000)
	if 382<i<=450:
		length=random.randint(6200,10000)
	if 450<i<=510:
		length=random.randint(5900,6200)
	if 510<i<=802:
		length=random.randint(2500,5900)
	if 802<i<=1248:
		length=random.randint(1000,2500)
	if 1248<i<=1500:
		length=random.randint(500,1000)
	if 1500<i<=1871:
		length=random.randint(400,500)
	if 1871<i<=2009:
		length=random.randint(350,400)
	if 2009<i<=3069:
		length=random.randint(300,350)
	if 3069<i<=3675:
		length=random.randint(200,300)
	if 3675<i<=4199:
		length=random.randint(150,200)
	if 4199<i<=5177:
		length=random.randint(100,150)
	if 5177<i<=6118:
		length=random.randint(75,100)
	if 6118<i<=7077:
		length=random.randint(60,75)
	if 7077<i<=8131:
		length=random.randint(50,60)
	return length

def random_chr():
	i=random.randint(1,200)
	if 1<=i<=32:
			length=random.randint(1,2)
	if 32<i<=60:
			length=random.randint(3,4)
	if 60<i<=108:
			length=random.randint(5,8)
	if 108<i<=128:
			length=random.randint(9,10)
	if 128<i<=136:
			length=11
	if 136<i<=150:
			length=random.randint(12,13)
	if 150<i<=162:
			length=random.randint(14,15)
	if 162<i<=174:
			length=random.randint(16,18)
	if 174<i<=186:
			length=random.randint(19,22)
	if 186<i<=200:
			length=23
	return length

def random_seq(i):
	iii=0
	seq=''
	while iii<i:
		iii+=1
		r=random.randint(1,4)
		if r==1:
			seq+='A'
		if r==2:
			seq+='T'
		if r==3:
			seq+='C'
		if r==4:
			seq+='G'
	return seq

def invert_seq(a):
	aa=''
	for i in range(len(a)):
		if a[-1-i]=='T':
			aa+='A'
		if a[-1-i]=='A':
			aa+='T'
		if a[-1-i]=='C':
			aa+='G'
		if a[-1-i]=='G':
			aa+='C'
		if a[-1-i] not in 'ATCG':
			aa+=a[-1-i]
	if len(a)!=len(aa):
		print('length not correct: '+str(len(a))+'\t'+str(len(aa)))
	return aa


def Takepos(a):
	return int(a.split('\t')[1])


def makesv():
	print("----------------------------makeSV----------------------------")
	f=open('hg38.fa','r')
	refs=f.read().split('>')[1:24]
	f.close()
	ref=[]
	names=[]
	for c in refs:
		seq=''
		seqparts=c.split('\n')[1:-1]
		names+=[c.split('\n')[0].split(' ')[0]]
		for d in seqparts:
			seq+=d
		ref+=[seq]
	
	numbers1=[800,800,700,700,600,600,600,600,500,500,400,400,400,300,300,200,200,200,150,150,150,150,600]
	numbers2=[80,80,70,70,60,60,60,60,50,50,40,40,40,30,30,20,20,20,15,15,15,15,60]
		
	for i in range(23):
		allsv=[]
		if i<22:
			chrom='chr'+str(i+1)
		else:
			chrom='chrX'
		reference=ref[i]

		insertions=[]
		while len(insertions)<numbers1[i]:
			po=random.randint(200000,len(reference)-200000)
			length=random_length('insertion')
			testif=0
			for c in allsv:
				if min(int(c.split('\t')[1])+int(c.split('\t')[2]),po+length) - max(int(c.split('\t')[1]),po) >=-1000 or 'N' in reference[1000+po:po+length+1000]:
				#if int(c.split('\t')[1])-1000<po<int(c.split('\t')[1])+int(c.split('\t')[2])+1000 or int(c.split('\t')[1])-1000<po+length<int(c.split('\t')[1])+int(c.split('\t')[2])+500 or 'N' in reference[500+po:po+length+500] or po-1000<int(c.split('\t')[1])<po+length+1000 or po-1000<int(c.split('\t')[1])+int(c.split('\t')[2])<po+length+1000:
					testif=1;  break
			if testif==0:
				hetero=random.randint(1,3)
				print(chrom+'\t'+str(po)+'\t'+str(length)+'\t'+str(hetero)+'\tinsertion')
				insertions+=[chrom+'\t'+str(po)+'\t'+str(length)+'\t'+str(hetero)+'\tinsertion']
				allsv+=[chrom+'\t'+str(po)+'\t'+str(length)+'\t'+str(hetero)+'\tinsertion']
		print(chrom+'\t'+str(len(insertions)))
		deletions=[]
		while len(deletions)<numbers1[i]:
			po=random.randint(200000,len(reference)-200000)
			length=random_length('deletion')
			testif=0
			for c in allsv:
				if min(int(c.split('\t')[1])+int(c.split('\t')[2]),po+length) - max(int(c.split('\t')[1]),po) >=-1000 or 'N' in reference[1000+po:po+length+1000]:
					testif=1;  break
			if testif==0:
				hetero=random.randint(1,3)
				print(chrom+'\t'+str(po)+'\t'+str(length)+'\t'+str(hetero)+'\tdeletion')
				deletions+=[chrom+'\t'+str(po)+'\t'+str(length)+'\t'+str(hetero)+'\tdeletion']
				allsv+=[chrom+'\t'+str(po)+'\t'+str(length)+'\t'+str(hetero)+'\tdeletion']
		print(chrom+'\t'+str(len(deletions)))
	
		inversions=[]
		while len(inversions)<numbers2[i]:
			po=random.randint(200000,len(reference)-200000)
			length=random_length('inversion')
			testif=0
			for c in allsv:
				if min(int(c.split('\t')[1])+int(c.split('\t')[2]),po+length) - max(int(c.split('\t')[1]),po) >=-1000 or 'N' in reference[1000+po:po+length+1000]:
					testif=1;  break
			if testif==0:
				hetero=random.randint(1,3)
				print(chrom+'\t'+str(po)+'\t'+str(length)+'\t'+str(hetero)+'\tinversion')
				inversions+=[chrom+'\t'+str(po)+'\t'+str(length)+'\t'+str(hetero)+'\tinversion']
				allsv+=[chrom+'\t'+str(po)+'\t'+str(length)+'\t'+str(hetero)+'\tinversion']
		print(chrom+'\t'+str(len(inversions)))

		duplications=[]
		while len(duplications)<numbers2[i]:
			po=random.randint(200000,len(reference)-200000)
			length=random_length('duplication')		
			testif=0
			for c in allsv:
				if min(int(c.split('\t')[1])+int(c.split('\t')[2]),po+length) - max(int(c.split('\t')[1]),po) >=-1000 or 'N' in reference[1000+po:po+length+1000]:
					testif=1;  break
			if testif==0:
				hetero=random.randint(1,3)
				print(chrom+'\t'+str(po)+'\t'+str(length)+'\t'+str(hetero)+'\tduplication')
				duplications+=[chrom+'\t'+str(po)+'\t'+str(length)+'\t'+str(hetero)+'\tduplication']
				allsv+=[chrom+'\t'+str(po)+'\t'+str(length)+'\t'+str(hetero)+'\tduplication']
		print(chrom+'\t'+str(len(duplications)))

		f=open('random_allsv_hetero','a')
		for c in allsv:
			f.write(c+'\n')
		f.close()
	return 0


def maketra():
	print("----------------------------maketra----------------------------")
	f=open('random_allsv_hetero','r')
	allsv=f.read().split('\n')[:-1]
	f.close()


	f=open('hg38.fa','r')
	refs=f.read().split('>')[1:24]
	f.close()
	ref=[]
	names=[]
	for c in refs:
		seq=''
		seqparts=c.split('\n')[1:-1]
		names+=[c.split('\n')[0].split(' ')[0]]
		for d in seqparts:
			seq+=d
		ref+=[seq]

	
	numbers3=[16,16,14,14,12,12,12,12,10,10,8,8,6,6,6,6,4,4,4,4,2,2,12]
	alltra=[]
	for i in range(23):
		translocations=[]
		if i<22:
			chr1='chr'+str(i+1)
		else:
			chr1='chrX'
		while len(translocations)<numbers3[i]:
			chr22=random_chr()
			pos1=random.randint(200000,len(ref[i])-200000)
			pos2=random.randint(200000,len(ref[chr22-1])-200000)

			testif=0
			
			if chr22<=22:
				chr2='chr'+str(chr22)
			else:
				chr2='chrX'

			if chr1==chr2:
				testif=1
			
			if testif==0:
				bp1region=ref[i][pos1-1000:pos1+1000]
				bp2region=ref[chr22-1][pos2-1000:pos2+1000]
				if 'N' in bp1region or 'N' in bp2region:
					testif=1
			if testif==0:
				for c in allsv:
					if (chr1==c.split('\t')[0] and int(c.split('\t')[1])-1000<=pos1<=int(c.split('\t')[1])+int(c.split('\t')[2])+1000) or (chr2==c.split('\t')[0] and int(c.split('\t')[1])-1000<=pos2<=int(c.split('\t')[1])+int(c.split('\t')[2])+1000):
						testif=1;  break
			if testif==0:
				for c in alltra:
					if (chr1==c.split('\t')[0] and abs(int(c.split('\t')[1])-pos1)<=5000) or (chr1==c.split('\t')[2] and abs(int(c.split('\t')[3])-pos1)<=5000) or (chr2==c.split('\t')[0] and abs(int(c.split('\t')[1])-pos2)<=5000) or (chr2==c.split('\t')[2] and abs(int(c.split('\t')[3])-pos2)<=5000):
						testif=1;  break
					chr11=c.split('\t')[0]
					chr22=c.split('\t')[2]
					pos11=int(c.split('\t')[1])
					pos22=int(c.split('\t')[3])
					if chr1==chr11 and chr2==chr22:
						if ( pos1<pos11 and pos2>pos22) or (pos1>pos11 and pos2<pos22):
							testif=1;  break
					if chr1==chr22 and chr2==chr11:
						if (pos1<pos22 and pos2>pos11) or (pos1>pos22 and pos2<pos11):					
							testif=1;  break

			if testif==0:
				hetero=random.randint(1,3)
				translocations+=[chr1+'\t'+str(pos1)+'\t'+chr2+'\t'+str(pos2)+'\t'+str(hetero)+'\ttranslocation']
				alltra+=[chr1+'\t'+str(pos1)+'\t'+chr2+'\t'+str(pos2)+'\t'+str(hetero)+'\ttranslocation']
		f=open('translocation-random','a')
		for c in translocations:
			f.write(c+'\n')
		f.close()

	return 0

def modif():
	print("----------------------------modif----------------------------")
	f=open('random_allsv_hetero','r')
	allsv=f.read().split('\n')[:-1]
	f.close()

	f=open('hg38.fa','r')
	refs=f.read().split('>')[1:24]
	f.close()
	ref=[]; names=[]
	for c in refs:
		seq=''
		seqparts=c.split('\n')[1:-1]
		names+=[c.split('\n')[0].split(' ')[0]]
		for d in seqparts:
			seq+=d
		ref+=[seq]



	print ('Finish Reading') 

	for i in range(23):
		chrom=names[i]
		g=open('output','a')
		f=open('modifiedref-deletion-'+chrom,'w')
		seq1=ref[i]
		seq2=ref[i]
		originlen=len(seq1)
		insertpart=[]
		length1=0
		length2=0
		for c in allsv:
			if c.split('\t')[0]==chrom:
				insertpart+=[c]
		insertpart.sort(key=Takepos,reverse=True)
		h=open('allsv-random-sorted-rev','a')
		for c in insertpart:
			h.write(c+'\n')
		h.close() 
		print ('Modify Ref....')
		for c in insertpart:
			c=c.split('\t')
			po=int(c[1])
			if c[-1]=='deletion':
				stop=po+int(c[2])
				if c[3] in '134':
					length1-=int(c[2])
					seq1=seq1[:po]+seq1[stop:]
				if c[3] in '234':
					length2-=int(c[2])
					seq2=seq2[:po]+seq2[stop:]
				continue
			if c[-1]=='insertion':
				insertseq=random_seq(int(c[2]))
				if c[3] in '134':
					length1+=int(c[2])
					seq1=seq1[:po]+insertseq+seq1[po:]
				if c[3] in '234':
					length2+=int(c[2])
					seq2=seq2[:po]+insertseq+seq2[po:]
				continue
			if c[-1]=='inversion':
				stop=po+int(c[2])
				invertseq=invert_seq(seq1[po:stop])
				if c[3] in '134':
					seq1=seq1[:po]+invertseq+seq1[stop:]
				if c[3] in '234':
					seq2=seq2[:po]+invertseq+seq2[stop:]
				continue
			if c[-1]=='duplication':
				stop=po+int(c[2])
				if c[3] in '134':
					length1+=int(c[2])
					seq1=seq1[:stop]+seq1[po:]
				if c[3] in '234':
					length2+=int(c[2])
					seq2=seq2[:stop]+seq2[po:]
				continue
			print ('Wrong SV found:  '+str(c))
		g.write('after modify,halploid 1 is corret :  '+str(len(seq1)-originlen==length1)+'\n')
		g.write('after modify,2 sequence length is :  '+str(len(seq2)-originlen==length2)+'\n\n\n')
		print('modify length correct? ')
		print(length1-len(seq1)==0-len(ref[i]))
		print(length2-len(seq2)==0-len(ref[i]))
		
		f.write('>'+chrom+'_1\n')
		iii=0
		while iii<len(seq1)/70:
			c=seq1[70*iii:70*(iii+1)]
			f.write(c+'\n')
			iii+=1
		c=seq1[70*iii:]
		f.write(c+'\n')
		f.write('>'+chrom+'_2\n')
		iii=0
		while iii<len(seq2)/70:
			c=seq2[70*iii:70*(iii+1)]
			f.write(c+'\n')
			iii+=1
		c=seq2[70*iii:]
		f.write(c+'\n')
		f.close()
	return 0	
	
def makrref():
	print("----------------------------makrref----------------------------")
	#os.system('cat modifiedref-deletion-chr* > modifiedref-allsv-beforetra.fa')

	reflist=open('filelist','r').read().split('\n')[:-1]
	f=open('modifiedref-allsv-beforetra.fa','w')
	for nnn in reflist:
		a=open(nnn,'r').read()
		f.write(a)
	f.close()


	f=open('modifiedref-allsv-beforetra.fa','r')
	allref=f.read().split('>')[1:]
	ref1=[]
	ref2=[]
	for i in range(23):
		ref1+=[clean_ref(allref[2*i])]
		ref2+=[clean_ref(allref[2*i+1])]


	f=open('translocation-modified-position','r')
	alltra=f.read().split('\n')[:-1]
	f.close()

	tra1=[]
	tra2=[]
	for c in alltra:
		if c.split('\t')[4] in ['1','3-1']:
			tra1+=[c];  continue
		if c.split('\t')[4] in ['2','3-2']:
			tra2+=[c];  continue
		print ('wrong translocation: '+c)

	print ('Finish reading')	
	print ('length of tra1:  '+str(len(tra1)))
	print ('length of tra2:  '+str(len(tra2)))
	f=open('modifiedref-allsv-aftertra.fa','w')
	g=open('output','w')
	donetra=0

	for i in range(23):
		chrom='chr'+str(i+1)
		if chrom=='chr23':
			chrom='chrX'
		seq1=''
		newpos=[i,0]
		nextbp=findnext(chrom,1,tra1)
		while nextbp!=[]:
			donetra+=1
			g.write(str(nextbp)+'\n')
			seq1+=ref1[newpos[0]][newpos[1]:nextbp[1]]
			try:
				newpos=[int(nextbp[2].split('hr')[1])-1,nextbp[3]]
			except:
				newpos=[22,nextbp[3]]
			nextbp=findnext(nextbp[2],nextbp[3]+100,tra1)
		seq1+=ref1[newpos[0]][newpos[1]:]
		f.write('>'+chrom+'_1\n')
		iii=0
		while iii<len(seq1)/70:
			c=seq1[70*iii:70*(iii+1)]
			f.write(c+'\n')
			iii+=1
		c=seq1[70*iii:]
		f.write(c+'\n')
	
	print ('done translocation:   '+str(donetra))
	print (donetra==len(tra1)*2)

	donetra=0
	for i in range(23):
		chrom='chr'+str(i+1)
		if chrom=='chr23':
			chrom='chrX'
		seq1=''
		newpos=[i,0]
		nextbp=findnext(chrom,1,tra2)
		while nextbp!=[]:
			
			donetra+=1
			g.write(str(nextbp)+'\n')
			seq1+=ref2[newpos[0]][newpos[1]:nextbp[1]]
			try:
				newpos=[int(nextbp[2].split('hr')[1])-1,nextbp[3]]
			except:
				newpos=[22,nextbp[3]]
			nextbp=findnext(nextbp[2],nextbp[3]+100,tra2)
		seq1+=ref2[newpos[0]][newpos[1]:]
		f.write('>'+chrom+'_2\n')
		iii=0
		while iii<len(seq1)/70:
			c=seq1[70*iii:70*(iii+1)]
			f.write(c+'\n')
			iii+=1
		c=seq1[70*iii:]
		f.write(c+'\n')

	g.close()
	print (donetra==len(tra2)*2)
	f.close()
	print (donetra)
	
	return 0


def findnext(chrom,pos,translocations):
	nextbp=[]
	for cc in translocations:
		c=cc.split('\t')
		if c[0]==chrom and int(c[1])>pos:
			if nextbp==[]:
				nextbp=[c[0],int(c[1]),c[2],int(c[3]),c[4]]
			else:
				if int(c[1])<nextbp[1]:
					nextbp=[c[0],int(c[1]),c[2],int(c[3]),c[4]]
		if c[2]==chrom and int(c[3])>pos:
			if nextbp==[]:
				nextbp=[c[2],int(c[3]),c[0],int(c[1]),c[4]]
			else:
				if int(c[3])<nextbp[1]:
					nextbp=[c[2],int(c[3]),c[0],int(c[1]),c[4]]
	return nextbp

		
def clean_ref(refseq):
	refseq=refseq.split('\n')[1:]
	refseq1=''
	for c in refseq:
		refseq1+=c
	return refseq1



if __name__ == "__main__":
	makesv()
	maketra()
	import sim_find_pos
	modif()
	makrref()
