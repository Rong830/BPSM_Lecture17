'''
Here's a list of made-up gene accession numbers:
xkn59438, yhdck2, eihd39d9, chdsye847, hedle3455, xjhd53e, 45da, de37dp
Write a program that will print only the accessions that satisfy the following criteria individually (i.e. treat each criterion separately):
contain the number 5
contain the letter d or e
contain the letters d and e in that order
contain the letters d and e in that order with a single letter between them
contain both the letters d and e in any order
start with x or y
start with x or y and end with e
contains any 3 numbers in any order
contains 3 different numbers in the accession
contain three or more numbers in a row
end with d followed by either a, r or p
'''
import re

gene = [ 'xkn59438', 'yhdck2', 'eihd39d9', 'chdsye847', 'hedle3455', 'xjhd53e', '45da', 'de37dp']

# contain the number 5
contain_5 = []
contain_de = []
contain_de_oder = []
contain_dnume = []
de_anyoder = []
xy = []
xy_end_e = []
num = []
diffnum = []
num_in_row = []
end_darp = []

for i in gene:
    if re.search('5', i):
        contain_5.append(i)
# contain the letter d or e   
    if re.search('[de]',i):
        contain_de.append(i)
# contain the letters d and e in that order       
    if re.search('de',i):
        contain_de_oder.append(i)
# contain the letters d and e in that order with a single letter between them
    if re.search('[de][0-9][de]',i):
        contain_dnume.append(i)
# contain both the letters d and e in any order
    if re.search('d',i) and re.search('e',i):
        de_anyoder.append(i)
# start with x or y
    if re.search('^[xy]',i):
        xy.append(i)
# start with x or y and end with e
    if re.search('^[cy]',i) and re.search('e$',i):
        xy_end_e.append(i)
# contains any 3 numbers in any order
    if len(re.findall('[0-9]',i)) >= 3:
        num.append(i)
# contains 3 different numbers in the accession
    if len(set(re.findall('[0-9]',i))) >= 3:
        diffnum.append(i)
# contain three or more numbers in a row
    if re.search('[0-9]{3,}',i):
        num_in_row.append(i)
# end with d followed by either a, r or p
    if re.search('d[arp]$',i):
        end_darp.append(i)

print('5',contain_5)
print('de',contain_de)
print('de order', contain_de_oder)




'''
DNA sequence: a double digest

A file called long_dna.txt has been put in the following directory:
/localdisk/data/BPSM/Lecture17

It contains a made-up DNA sequence.
What fragment lengths will we get if we digest the sequence with a novel restriction enzyme BpsmI, whose recognition site is ANT*AAT, where * indicates the position of the cut site.
What will the fragment lengths be if we do a double digest with both BpsmI and BpsmII (whose recognition site is GCRW*TG)?
What are the sequences of the fragments themselves?
'''
import os
dna = open('long_dna.txt').read().rstrip()

# What fragment lengths will we get if we digest the sequence with a novel restriction enzyme BpsmI, whose recognition site is ANT*AAT, where * indicates the position of the cut site.
num_frag = len(re.findall('A[ATCG]TAAT',dna)) + 1
# What will the fragment lengths be if we do a double digest with both BpsmI and BpsmII (whose recognition site is GCRW*TG)?
num_frag2 = len(re.findall('A[ATCG]TAAT',dna) + re.findall('GC[AG][AT]TG',dna)) + 1

# What are the sequences of the fragments themselves?
find = []
find1 = re.finditer('A[ATCG]TAAT',dna)
find2 = re.finditer('GC[AG][AT]TG',dna)
counter=0
for i in find1:
    counter += 1 
    start = i.start()
    end = i.end()


BpsmI='A[GATC]TAAT'
print('BpsmI cuts at:',BpsmI) 
# Find the sites, incrementing by three to account for where the enzyme cuts in the recognition sequence
for matching in re.finditer(BpsmI, dna): 
    print(matching.start()+3) 


# Start: open and read the file
dna = open('/localdisk/data/BPSM/Lecture17/long_dna.txt').read().rstrip('\n') 
last_cut = 0
findnum=0
for matching in re.finditer(BpsmI, dna):
    findnum += 1
    cut_position = matching.start() + 3
# Distance from the current cut site to the previous one
    fragment_size = cut_position - last_cut
    print('Fragment size is ' + str(fragment_size))
    last_cut = cut_position
# We also have to remember the last fragment, from the last cut to the end:
    if findnum == len(list(re.finditer(BpsmI, dna))) :
       fragment_size = len(dna) - last_cut
       print('Fragment size is ' + str(fragment_size))



# What will the fragment lengths be if we do a double digest with both BpsmI and BpsmII (whose recognition site is GCRW*TG)?
BpsmI='A[GATC]TAAT'
BpsmII='GC[AG][AT]TG'
# Make a list to store the cut positions for both enzymes
all_cuts = []
# Add cut positions for BpsmI
for match in re.finditer(BpsmI, dna): 
    all_cuts.append(match.start() + 3) 

# Add cut positions for BpsmII
for match in re.finditer(BpsmII, dna): 
    all_cuts.append(match.start() + 4)
all_cuts.sort()
print(all_cuts)

# Double digest run
last_cut = 0
counter = 0
for cut_position in all_cuts:
    counter +=1
    fragment_size = cut_position - last_cut
    print('Fragment '+str(counter)+' size is ' + \
       str(fragment_size) +': '+ str(last_cut)+ ' to ' +str(cut_position) )
    last_cut = cut_position
# Now the last fragment
fragment_size = len(dna) - last_cut
counter +=1
print('Fragment '+str(counter)+' size is ' + \
  str(fragment_size) +': '+ str(last_cut)+ ' to ' +str(len(dna)) )

  # Let's use a dictionary to store the fragment sequences
fragment_sequences = {}


# Double digest run
last_cut = 0
counter = 0
for cut_position in all_cuts:
    counter +=1
    fragment_size = cut_position - last_cut
    print('Fragment '+str(counter)+' size is ' + \
       str(fragment_size) +': '+ str(last_cut)+ ' to ' +str(cut_position) )
# Get the sequence substring
    fragment_sequences['Fragment'+str(counter)] = dna[last_cut:cut_position]
    print(fragment_sequences['Fragment'+str(counter)])
# Get the fragment start and end
    fragends = dna[last_cut:cut_position][0:6] + '...' + dna[last_cut:cut_position][-6:]
    print('Fragment '+str(counter)+ ' has ends: '+fragends+'\n')
    last_cut = cut_position
# Now the last fragment
fragment_size = len(dna) - last_cut
counter +=1
print('Fragment '+str(counter)+' size is ' + \
  str(fragment_size) +': '+ str(last_cut)+ ' to ' +str(len(dna)) )
fragment_sequences['Fragment'+str(counter)] = dna[last_cut:]
print(fragment_sequences['Fragment'+str(counter)])
fragends = dna[last_cut:][0:6] + '...' + dna[last_cut:][-6:]
print('Fragment has ends: '+fragends)
