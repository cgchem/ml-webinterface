import os

from flask import Blueprint, flash, g, redirect, render_template, request, session, url_for, jsonify, current_app, send_from_directory

from . import predict_2

bp = Blueprint('lookup_page', __name__, url_prefix='')


def get_models():

    return [os.path.splitext(x)[0] for x in os.listdir(r".\main\models")]



@bp.route('/', methods=('GET', 'POST'))
def begin():

    error = None

    if request.method == 'POST':

        query = request.form['qry']
        smiles_list = query.rstrip("\n").split("\n")
        result, file_location = predict_2.predict_from_smiles("solubility", smiles_list)
        file_link = "/result/" + file_location
        return render_template('page1/page1.html', results=result, model_list=get_models(), file_loc=file_link)

    if request.method == 'GET':

        get_models()
        pass

    if error:
        flash(error)
    return render_template('page1/page1.html', model_list=get_models())


@bp.route('/result/<fname>', methods=('GET',))
def download_results(fname):

    return send_from_directory(os.path.join("static", "output"), fname, as_attachment=True)


@bp.route('/predict/<model>/<smiles>', methods=('GET',))
def url_prediction(model, smiles):

    result = predict_2.predict_single(model, smiles)

    return jsonify(result)
