import csv
 f= open("target1.faa","w+")
 e= open("target0.faa","w+")
 with open( 'training_set.csv', 'r' ) as theFile:
     reader = csv.DictReader(theFile)
     for line in reader:
         if line['label']=='1':
              f.write('>'+line.get('id')+'\n'+line.get('seq')+'\n')
         else:
             e.write('>'+line.get('id')+'\n'+line.get('seq')+'\n') 
 f.close()
 e.close()
