# This script modifies the complete mitochondrial nt database to be relavant for species in California

spp = []

for line in open("ca_birds.txt"):
        i=line.strip().split(',')
        if len(i)<4:
                continue
        if len(i[1].split(":"))==2:
                continue
        if "Rare" in i[3]:
		continue
	bird="_".join(i[2].split())
        spp.append(bird)

for line in open("ca_mammal.txt"):
        i=line.strip().split()
        if len(i)<2:
                continue
        mam="_".join([i[0],i[1]])
        spp.append(mam)

for line in open("ca_reptile.txt"):
        if line.strip()=="Species":
		continue
	line = line.strip().strip("[edit]")
	line = line.strip("[e]")
	line = line.strip("[i]")
	i=line.strip().split()
        if i[0]=="Family":
		continue
	if len(i)<2:
                continue
        rep="_".join([i[0],i[1]])
        spp.append(rep)
	#print rep

add=['Aedes_aegypti', 'Aedes_albopictus', 'Aedes_notoscriptus', 'Culex_quinquefasciatus', 'Culex_pipiens', 'Equus_caballus','Gallus_gallus','Capra_hircus','Bos_taurus','Homo_sapiens','Canis_lupus']

# remove unlikely species. e.g. bison and raven
bad=['Pinicola_enucleator','Nucifraga_columbiana','Corvus_corax', 'Bison_bison', 'Sitta_carolinensis', ]
#Sitta_carolinensis 2 sightings in 2018

for item in add:
	spp.append(item)

for b in bad:
	spp.remove(b)


switch=0
for line in open("mito.nt"):
	if line[0]==">":
		i=line.strip().split()
		name = "_".join([i[1],i[2]])
		if name not in spp:
			switch=1
			continue
                if "neanderthalensis" in line.strip(): # skip prehistoric homonids
			switch=1
			continue

                else:	
			switch=0
			print line.strip()
			continue
	if not switch:
		print line.strip()


