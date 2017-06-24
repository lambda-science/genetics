from flask import Flask, flash, render_template, redirect, request
from flask_wtf import FlaskForm
from wtforms import TextAreaField
from wtforms.validators import DataRequired
from genetics import *

app = Flask(__name__)
app.secret_key = 'myverylongsecretkey'

class SequenceForm(FlaskForm):
	sequenceADN = TextAreaField('sequenceADN')

globalSequenceADN = ""

@app.route("/", methods=['GET', 'POST'])
def index():
	form = SequenceForm()
	global globalSequenceADN
	if form.validate_on_submit():
		if isSequenceValid(form.sequenceADN.data):
			globalSequenceADN = form.sequenceADN.data
			return redirect('/resultats')
		else:
			flash('SEQUENCE INVALIDE')
	else:
		flash('FORMULAIRE VIDE')
	return render_template('index.html', title ='Acceuil:', form=form)

@app.route('/resultats')
def resultat():
	compo = composition(globalSequenceADN)
	gc = int(pourcentGC(globalSequenceADN))
	fusion = int(tempFusionHowley(globalSequenceADN))
	arn = ADN2ARN(globalSequenceADN)
	prot = traduction(arn)
	return render_template('resultats.html', title ='Resultats', composition=compo, gc=gc, fusion=fusion, arn=arn, prot=prot, length=len(globalSequenceADN))

if __name__ == "__main__":
	app.run(debug=True)