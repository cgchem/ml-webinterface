import os

from flask import Blueprint, flash, g, redirect, render_template, request, session, url_for, jsonify, current_app, send_from_directory

from . import predict_2

bp = Blueprint('lookup_page', __name__, url_prefix='')


def get_models():

    """Gets a list of the names of all pickled model files in the model directory, to present in dropddown list"""

    return [os.path.splitext(x)[0] for x in os.listdir(r"./main/models")]


@bp.route('/', methods=('GET', 'POST'))
def begin():

    """Main page method for interactive prediction via textentry box"""

    error = None

    if request.method == 'POST':

        query = request.form['qry']
        mod = request.form['model']
        print(f"Predicting values using {mod} model")
        smiles_list = query.rstrip("\r\n").split(os.linesep)
        smiles_list = [x.strip("\r\n") for x in smiles_list]  # get rid of the damn extra newlines!!
        result, file_location = predict_2.predict_multiple(mod, smiles_list)
        file_link = "/result/" + file_location
        return render_template('page1/page1.html', results=result, model_list=get_models(), file_loc=file_link)

    if request.method == 'GET':

        pass  # nothing extra to do beyond the standard render template

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

    """Expects a POST of comma-delimited SMILES strings, like PubChem. Returns a JSON dict of smiles: value
    inside another dict that gives the name of the prediction model that was passed in the URL.
    e.g {"solubility":{"CCN":1.0799283637300687,"CNC":1.1766386397282302}}"""

    try:
        smiles_list = request.data.decode("UTF-8").split(",")  # data arrives in binary format, convert to list of strs
    except:
        return jsonify({"Fault": "The text of the request could not be decoded."})

    result, _ = predict_2.predict_multiple(model, smiles_list, write_file=False)

    outdc = {}  # where we will store the predictions from the results zip object
    outer = {model:outdc}

    for smi, pred in result:  # a zip of tuples (smiles, prediction)
        outdc[smi] = pred  # rather than tuples we will use them to give key:value in a dict

    return jsonify(outer)
