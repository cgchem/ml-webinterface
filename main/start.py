import os

from flask import Blueprint, flash, g, redirect, render_template, request, session, url_for, jsonify, current_app, send_from_directory

from . import predict_2

bp = Blueprint('lookup_page', __name__, url_prefix='')


def get_models():

    return [os.path.splitext(x)[0] for x in os.listdir(r"./main/models")]


@bp.route('/', methods=('GET', 'POST'))
def begin():

    error = None

    if request.method == 'POST':

        query = request.form['qry']
        mod = request.form['model']
        print(f"Model is {mod}")
        smiles_list = query.rstrip("\r\n").split(os.linesep)
        smiles_list = [x.strip("\r\n") for x in smiles_list]  # get rid of the damn extra newlines!!
        result, file_location = predict_2.predict_from_smiles(mod, smiles_list)
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


@bp.route('/predict/<model>', methods=('POST',))
def post_prediction(model):

    """Expects a POST of comma-delimited SMILES strings, like PubChem"""

    try:
        smiles_list = request.data.decode("UTF-8").split(",")  # data arrives in binary format, convert to list of strs
        print(smiles_list)
    except:
        return jsonify({"Fault": "The text of the request could not be decoded."})

    result, _ = predict_2.predict_from_smiles(model, smiles_list)  # don't care about the file name, not using it

    outdc = {}  # where we will store the predictions from the results zip object
    outer = {model:outdc}

    for smi, pred in result:  # a zip of tuples (smiles, prediction)
        outdc[smi] = pred  # rather than tuples we will use them to give key:value in a dict

    return jsonify(outer)
