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
		if isSequenceValid(form.sequenceADN.data):
			session['SequenceADN'] = form.sequenceADN.data
			session['simpleBool'] = form.motifsimpleBool.data
			session['motifSimple'] = form.motifsimple.data
			session['complexeBool'] = form.motifcomplexeBool.data
			session['motifComplexe'] = form.motifcomplexe.data
			session['signaBool'] = form.signatureBool.data
			session['signatureSize'] = form.signature.data
			return redirect('/resultats')
		else:
			flash('SEQUENCE INVALIDE')
	else:
		flash('FORMULAIRE VIDE')
	return render_template('index.html', title ='Acceuil:', form=form)

@app.route('/resultats')
def resultat():
	compo = composition(session['SequenceADN'])
	gc = int(pourcentGC(session['SequenceADN']))
	fusion = int(tempFusionHowley(session['SequenceADN']))
	arn = ADN2ARN(session['SequenceADN'])
	prot = traduction(arn)
	posMotifSimple = ["N/A"]
	posMotifComplexe = ["N/A"]
	sortedSigna = ["N/A"]
	if session['simpleBool']:
		posMotifSimple = localiserMotifSimple(session['SequenceADN'],session['motifSimple'], 0)

	if session['complexeBool']:
		posMotifComplexe = localiserMotifSimple(session['SequenceADN'],session['motifComplexe'], 0)

	if session['signaBool']:
		signa = signature(session['SequenceADN'], session['signatureSize'])
		sortedSigna = sorted(signa.items(), key=lambda x: x[1], reverse=True)
	return render_template('resultats.html', title ='Resultats', composition=compo, gc=gc, fusion=fusion, arn=arn, prot=prot,
	length=len(session['SequenceADN']),posMotifSimple= posMotifSimple, posMotifComplexe=posMotifComplexe, signa=sortedSigna, session=session)

if __name__ == "__main__":
	app.run(debug=True)