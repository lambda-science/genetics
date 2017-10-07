from flask import Flask, flash, render_template, redirect, request, session
from flask_wtf import FlaskForm
from wtforms import TextAreaField, BooleanField, IntegerField, StringField
from wtforms.validators import DataRequired, Optional
from genetics import *

app = Flask(__name__)
app.secret_key = 'myverylongsecretkey'

class SequenceForm(FlaskForm):
	sequenceADN = TextAreaField('sequenceADN')
	motifsimpleBool = BooleanField('motifsimpleBool', default=False)
	motifsimple = StringField('motifsimple', validators=[Optional()])
	motifcomplexeBool = BooleanField('motifcomplexeBool', default=False)
	motifcomplexe = StringField('motifcomplexe', validators=[Optional()])
	signatureBool = BooleanField('signature', default=False)
	signature = IntegerField('signature', validators=[Optional()])

@app.route("/", methods=['GET', 'POST'])
def index():
	form = SequenceForm()
	if form.validate_on_submit():
		if isSequenceValid(fastaConvert(form.sequenceADN.data)[1]):
			session['SequenceADN'] = form.sequenceADN.data
			session['simpleBool'] = form.motifsimpleBool.data
			session['motifSimple'] = form.motifsimple.data
			session['complexeBool'] = form.motifcomplexeBool.data
			session['motifComplexe'] = form.motifcomplexe.data
			session['signaBool'] = form.signatureBool.data
			session['signatureSize'] = form.signature.data
			return redirect('/resultats')
		else:
			flash(repr(fastaConvert(form.sequenceADN.data)[1]))
	else:
		flash('FORMULAIRE VIDE')
	return render_template('index.html', title ='Acceuil:', form=form)

@app.route('/resultats')
def resultat():
	seqTitle = fastaConvert(session['SequenceADN'])[0]
	sequence = fastaConvert(session['SequenceADN'])[1]
	compo = composition(sequence)
	gc = int(pourcentGC(sequence))
	fusion = int(tempFusionHowley(sequence))
	arn = ADN2ARN(sequence)
	prot = traduction(arn)
	posMotifSimple = ["N/A"]
	posMotifComplexe = ["N/A"]
	sortedSigna = ["N/A"]
	if session['simpleBool']:
		posMotifSimple = localiserMotifSimple(sequence,session['motifSimple'], 0)

	if session['complexeBool']:
		posMotifComplexe = localiserMotifSimple(sequence,session['motifComplexe'], 0)

	if session['signaBool']:
		signa = signature(sequence, session['signatureSize'])
		sortedSigna = sorted(signa.items(), key=lambda x: x[1], reverse=True)
	return render_template('resultats.html', title ='Resultats', seqTitle = seqTitle, composition=compo, gc=gc, fusion=fusion, arn=arn, prot=prot,
	length=len(sequence),posMotifSimple= posMotifSimple, posMotifComplexe=posMotifComplexe, signa=sortedSigna, session=session)

if __name__ == "__main__":
	# For server developpement
	# app.run(debug=False, host="0.0.0.0", port=5010)
	# For local testing
	app.run(debug=True, host="127.0.0.1", port=5010)