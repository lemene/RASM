f=open('translocation-random','r')
alltra=f.read().split('\n')[:-1]
f.close()

f=open('random_allsv_hetero','r')
allsv=f.read().split('\n')[:-1]
f.close()

allsv=[c for c in allsv if 'inversion' not in c]

allsv1=[c for c in allsv if c.split('\t')[3] in '134']
allsv2=[c for c in allsv if c.split('\t')[3] in '234']




modiftra=[]
while True:
	for cc in alltra:
		chr1=cc.split('\t')[0]
		pos1=int(cc.split('\t')[1])
		chr2=cc.split('\t')[2]
		pos2=int(cc.split('\t')[3])
		hetero=cc.split('\t')[4]

		#pos1:
		sv1=[c for c in allsv1 if c.split('\t')[0]==chr1 and int(c.split('\t')[1])<pos1]
		sv2=[c for c in allsv2 if c.split('\t')[0]==chr1 and int(c.split('\t')[1])<pos1]
		shift1=0
		shift2=0
		for c in sv1:
			if 'insertion' in c or 'duplication' in c:
				shift1+=int(c.split('\t')[2]);   continue
			if 'deletion' in c:
				shift1-=int(c.split('\t')[2]);   continue
			print ('One event wrong:   '+c)

		for c in sv2:
			if 'insertion' in c or 'duplication' in c:
				shift2+=int(c.split('\t')[2]);   continue
			if 'deletion' in c:
				shift2-=int(c.split('\t')[2]);   continue
			print ('One event wrong:   '+c)

		modified1pos1=pos1+shift1
		modified2pos1=pos1+shift2
		#pos2:
		sv1=[c for c in allsv1 if c.split('\t')[0]==chr2 and int(c.split('\t')[1])<pos2]
		sv2=[c for c in allsv2 if c.split('\t')[0]==chr2 and int(c.split('\t')[1])<pos2]
		shift1=0
		shift2=0
		for c in sv1:
			if 'insertion' in c or 'duplication' in c:
				shift1+=int(c.split('\t')[2]);   continue
			if 'deletion' in c:
				shift1-=int(c.split('\t')[2]);   continue
			print ('One event wrong:   '+c)
		for c in sv2:
			if 'insertion' in c or 'duplication' in c:
				shift2+=int(c.split('\t')[2]);   continue
			if 'deletion' in c:
				shift2-=int(c.split('\t')[2]);   continue
			print ('One event wrong:   '+c)
		modified1pos2=pos2+shift1
		modified2pos2=pos2+shift2

		if hetero=='1':
			modiftra+=[chr1+'\t'+str(modified1pos1)+'\t'+chr2+'\t'+str(modified1pos2)+'\t'+hetero+'\ttranslocation']
		if hetero=='2':
			modiftra+=[chr1+'\t'+str(modified2pos1)+'\t'+chr2+'\t'+str(modified2pos2)+'\t'+hetero+'\ttranslocation']
		if hetero in '34':
			modiftra+=[chr1+'\t'+str(modified1pos1)+'\t'+chr2+'\t'+str(modified1pos2)+'\t'+hetero+'-1\ttranslocation']
			modiftra+=[chr1+'\t'+str(modified2pos1)+'\t'+chr2+'\t'+str(modified2pos2)+'\t'+hetero+'-2\ttranslocation']
	break
f=open('translocation-modified-position','w')
for c in modiftra:
	f.write(c+'\n')
f.close()

