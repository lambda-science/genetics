import re

def isSequenceValid(sequenceADN):
	sequenceADN = set(list(sequenceADN.lower()))
	a = True
	for i in sequenceADN:
		if not (i in ['a', 'c', 't', 'g']):
			a = False
	return a

def composition(sequenceADN):
	composition={}
	for i in sequenceADN.lower():
		composition[i] = composition.get(i, 0) + 1
	return composition

def pourcentGC(sequenceADN):
	compositionADN = composition(sequenceADN)
	pourcentGC = 100*(compositionADN["g"] + compositionADN["c"])/(len(sequenceADN))
	return pourcentGC

def tempFusionHowley(sequenceADN):
	GC = pourcentGC(sequenceADN)
	Tm = 67.5 + (0.34*GC)-(395/len(sequenceADN))
	return Tm

def estComplementaire(brinA, brinB):
	codeComplementaire = {"a":"t", "t":"a", "g":"c", "c":"g"}
	brinB = brinB[::-1]
	complementaireBrinB = [ codeComplementaire[i] for i in brinB ]
	complementaireBrinB = ''.join(complementaireBrinB)

	if brinA == complementaireBrinB:
		return True
	else:
		return False

def ADN2ARN(sequenceADN):
	sequenceARN = []
	sequenceADN = sequenceADN.lower()
	for i in sequenceADN:
		if i == "t":
			sequenceARN.append("a")
		elif i == "a":
			sequenceARN.append("u")
		elif i == "c":
			sequenceARN.append("g")
		elif i == "g":
			sequenceARN.append("c")

	sequenceARN = ''.join(sequenceARN)
	return sequenceARN


def traduction(sequenceARN, cadre):
	code_genetique = {
	'uuu': 'F', 'ucu': 'S', 'uau': 'Y', 'ugu': 'C',
	'uuc': 'F', 'ucc': 'S', 'uac': 'Y', 'ugc': 'C',
	'uua': 'L', 'uca': 'S', 'uaa': '*', 'uga': '*',
	'uug': 'L', 'ucg': 'S', 'uag': '*', 'ugg': 'W',
	'cuu': 'L', 'ccu': 'P', 'cau': 'H', 'cgu': 'R',
	'cuc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
	'cua': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
	'cug': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
	'auu': 'I', 'acu': 'T', 'aau': 'N', 'agu': 'S',
	'auc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
	'aua': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
	'aug': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
	'guu': 'V', 'gcu': 'A', 'gau': 'D', 'ggu': 'G',
	'guc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
	'gua': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
	'gug': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G',
	}	


	if cadre in range(0,3):
		sequenceCodon = [ sequenceARN[j:j+3] for j in range(cadre, len(sequenceARN), 3) ]
		sequenceProteine = [ code_genetique[i] for i in sequenceCodon if i in code_genetique ]
		sequenceProteine = '-'.join(sequenceProteine)
		return sequenceProteine

	else:
		return "Erreur dans le choix du cadre"

def localiserMotifSimple(sequenceADN, motif, positionRecherche):
	motif = motif.lower()
	sequenceADN = sequenceADN.lower()
	sequenceADN = sequenceADN[positionRecherche::]
	pos = 0
	listePosMotif = []
	pos = sequenceADN.find(motif)
	while pos != -1:
		listePosMotif.append(pos)
		pos = sequenceADN.find(motif, pos+1)
	return listePosMotif

def localiserMotifRE(sequenceADN, motif, positionRecherche):
	sequenceADN = sequenceADN.upper()
	pattern = re.compile(motif)
	listPosMotif = []
	m = pattern.search(sequenceADN, positionRecherche)
	while m != None:
		listPosMotif.append(m.start(0))
		m = pattern.search(sequenceADN, m.start(0)+1)
	return listPosMotif

def signature(sequenceADN, m):
	sequenceADN = sequenceADN.lower()
	dictSignature = {}	
	if m in range(0,11):
		for i in range(0, len(sequenceADN), 1):
			dictSignature[sequenceADN[i:i+m]] = dictSignature.get(sequenceADN[i:i+m], 0) + 1
		cleInvalide = [ j for j in dictSignature if len(j) != m ]
		for k in cleInvalide:
			del dictSignature[k]
		return dictSignature
	else:
		return "Erreur: motif trop grand"


if __name__ == '__main__':
	s = "cctagctagctacgtATGatgcatgctacgatgTAGgtcatcgatcatgTAGCTAGC"
	s2 = "gtctactagctctga"
	resultats = {
	'ADN': s, 'composition': composition(s), 'GC': pourcentGC(s),
	'fusion': tempFusionHowley(s), 'complementaire': estComplementaire(s, s2), 'ARN': ADN2ARN(s),
	'motifsimple': localiserMotifSimple(s, "AATTGC", 7),
	'motifRE': localiserMotifRE(s, "[AT][GC]..AT*", 4), 'signature': signature(s, 3)
	}

	print("Sequence de test: ", resultats['ADN'],  '\n',
	"Composition:", resultats['composition'], '\n',
	"Pourcent de GC: %d " % resultats['GC'], "% \n",
	"Température de fusion: %d" % resultats['fusion'], "°C \n",
	"Test de complémentarité avec la sequence :", s2, ": ", resultats['complementaire'], '\n',
	"Sequence ARN correspondante: ", resultats['ARN'], '\n',
	"Motif AATTGC a partir de la position 7 aux positions: ", resultats['motifsimple'], '\n',
	"Motif regulier [AT][GC]..AT* a partir de la position 4 aux positions:", resultats['motifRE'], '\n',
	"Nombre de motif de taille 3:", len(resultats['signature']), '\n',
	"Liste et nombre de ces motif:", resultats['signature'])