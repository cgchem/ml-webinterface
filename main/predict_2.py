import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import numpy as np
import pandas as pd
import datetime

PROPERTY_NAMES = list(rdMolDescriptors.Properties.GetAvailableProperties())
PROPERTY_GETTER = rdMolDescriptors.Properties(PROPERTY_NAMES)


def smi2props(smi):

    mol = Chem.MolFromSmiles(smi)
    props = None
    if mol:
        props = np.array(PROPERTY_GETTER.ComputeProperties(mol))
    return props


def unique_name():

    return datetime.datetime.now().strftime("%Y%m%d%H%M%S")


def predict_from_smiles(model_name, smiles_list):

    """Takes a SMILES list and returns a DataFrame with the SMILES, their properties, and the prediction value.
     Predicts based on the model name which must match a pkl file in the models directory, this will always be true
     from the web interface."""

    with open("main/models/" + model_name + ".pkl", "rb") as f:
        model = pickle.load(f)

    df = pd.DataFrame({"SMILES": smiles_list})  # model expects a DataFrame
    df["props"] = df["SMILES"].apply(smi2props, 1)
    df.dropna(inplace=True)  # throw out rows where we got None for property calculation because the SMILES was bad
    df[PROPERTY_NAMES] = df['props'].to_list()  # PROPERTY_NAMES is a list of column headings and we are expanding
    # the list that is in the "props" cell

    vals = df[PROPERTY_NAMES].astype(float)

    result_list = model.predict(vals)
    df["prediction"] = result_list

    return df  # various functions that call this function can do what they want with the dataframe


def predict_single(model_name, smiles):

    """For the URL API."""

    pred = predict_from_smiles(model_name, [smiles])  # pass our single SMILES as a list of length 1
    prediction = pred["prediction"][0]  # there is just one result in a dataframe with one row

    return prediction


def predict_multiple(model_name, smiles_list, write_file=True):

    df = predict_from_smiles(model_name, smiles_list)
    prediction = zip(list(df["SMILES"]), df["prediction"])
    # we are getting the SMILES back out of the df because we may have discarded some invalid ones

    if write_file:  # set up a new df with minimal columns to write to a csv
        df2 = pd.DataFrame()
        df2["SMILES"] = df["SMILES"]
        df2[model_name] = df["prediction"]
        fname = model_name + "_" + unique_name() + ".csv"
        df2.to_csv("main/static/output/" + fname, index=False)  # put here so they can be served as static files
    else:
        fname = None  # may as well be consistent and always return two things

    return prediction, fname  # provide the name where file is saved so the user can download it
