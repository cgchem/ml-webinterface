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

    with open("main\\models\\" + model_name + ".pkl", "rb") as f:
        model = pickle.load(f)

    df = pd.DataFrame({"SMILES": smiles_list})  # model expects a DataFrame
    df["props"] = df["SMILES"].apply(smi2props, 1)
    df.dropna(inplace=True)  # throw out rows where we got None for property calculation because the SMILES was bad
    df[PROPERTY_NAMES] = df['props'].to_list()  # PROPERTY_NAMES is a list of column headings and we are expanding
    # the list that is in the "props" cell

    vals = df[PROPERTY_NAMES].astype(float)

    result_list = model.predict(vals)

    prediction = zip(list(df["SMILES"]), result_list)
    # we are getting the SMILES back out of the df because we may have discarded some invalid ones

    df2 = pd.DataFrame()
    df2["SMILES"] = df["SMILES"]
    df2[model_name] = result_list
    fname = model_name + "_" + unique_name() + ".csv"
    df2.to_csv("main\\static\\output\\" + fname, index=False)

    return prediction, fname  # provide the name where file is saved so the user can download it


def predict_single(model_name, smiles):

    """For the URL API. Some code duplication for now."""

    #TODO: error messages like model name not recognised, etc
    with open("main\\models\\" + model_name + ".pkl", "rb") as f:
        model = pickle.load(f)

    df = pd.DataFrame({"SMILES": [smiles]})  # model expects a DataFrame
    df["props"] = df["SMILES"].apply(smi2props, 1)
    df.dropna(inplace=True)  # throw out rows where we got None for property calculation because the SMILES was bad
    df[PROPERTY_NAMES] = df['props'].to_list()  # PROPERTY_NAMES is a list of column headings and we are expanding
    # the list that is in the "props" cell

    vals = df[PROPERTY_NAMES].astype(float)

    result_list = model.predict(vals)

    prediction = result_list[0]  # just one result
    # we are getting the SMILES back out of the df because we may have discarded some invalid ones

    return prediction
