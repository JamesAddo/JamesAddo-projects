# Python script to generate ML regression QSAR models using Morgan Fingerprints and two machine-learning methods: random forest and naive Bayes using python, RDKit and scikit-learn machine learning functions
#
# Example: QSAR modeling of 5-HT1A Receptor Binding: Predictions of 5-HT1A Ligands Using Machine Learning Regression/Classification Models.
# This is the code for ML regression QSAR/From an sdf file/File contains SDF and activities information:

#
# Created by James Addo, July 2014
import cPickle, gzip, numpy, copy, math
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, rdmolops
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import rdMolDescriptors
from matplotlib import cm
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import BernoulliNB
from rdkit.Chem.Draw import IPythonConsole
from sklearn import cross_validation
from sklearn import metrics
from sklearn.cross_validation import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
from sklearn import preprocessing
import numpy as np
import matplotlib.pyplot as plt
import os,sys,csv
from time import time

# prepare command-line argument parser
parser = argparse.ArgumentParser()
parser.add_argument('-i', "-- sdf_input_file", nargs='?', help="input sdf file")
parser.add_argument('-t', "-- test_sdf_input_file ", nargs='?', help="input sdf file")
parser.add_argument('-o', "--output", type=str, default=sys.stdout)

# read the mols
# generate dict to store molecules/mols and activities
mols_dict = {}
suppl = Chem.SDMolSupplier('sdf_input_file')
for mol in suppl:
    # mols.key = mol
    # mols_dic[k]) = act_value
    act_value =  mol.GetProp('_Name')
    mols_dic = dict(zip(mol, act_value)) # mols_dic [k] = act_value             

# calculate training fingerprints    	
fps = []
acts = []
for k,v in mols_dic.items():       
    fp = AllChem.GetMorganFingerprintAsBitVect(k, 2)
    arr = np.zeros((1,))
    fps.append(fp)

print "morgan2 generated:", len(fps)
    acts.append(mols_dic[k])

acts_array = np.array(acts)
fps_array   = np.array(fps)

    rf = RandomForestClassifier(n_estimators=100, max_depth=10, min_samples_split=2, min_samples_leaf=1)
    rf.fit(acts_array, fps_array)
    # write the model
    cPickle.dump(rf, outfile, 2)
    i += 1
outfile.close()
print "done"
# get test fingerprints
# precalculate fingerprints for test molecules
testmols_dict = {}
suppl = Chem.SDMolSupplier('test_sdf_input_file')
for testmol in suppl:
    # testmols.key = test_mol
    # testmols_dic[k]) = test_act_value
    test_act_value =  mol.GetProp('_Name')
    test_mols_dic = dict(zip(test_mol, test_act_value)) # testmols_dic [k] = test_act_value             

test_fps = []
test_acts = []
for k,v in testmols_dic.items():       
    test_fp = AllChem.GetMorganFingerprintAsBitVect(k, 2)
    test_fps.append(test_fp)
    test_fps_array = np.array(test_fps)
    test_acts.append(testmols_dic [k])
    test_acts_array = np.array(test_acts)
  
# make predictions for test compounds
test_SVMpredictions = rf.predict(test_fps_array)

# cross-validate
scores = cross_validation.cross_val_score(rf, fps_array, acts_array, cv=5)
# scores = cross_validation.cross_val_score(rf, X_train,y_train, cv=cv_counter,score_func=metrics.zero_one_score)

######################## MAIN PART ###########################
if __name__=='__main__':
# read in command line options
args = parser.parse_args()
# required arguments
if args.test_sdf_input_file:
    sdf_input_file = args.sdf_input_file
else:
    raise RuntimeError('input sdf file missing')
print "RF training for input sdf file:", sdf_input_file 

sdf_input_file = args.sdf_input_file
test_sdf_input_file 
test_sdf_input_file  = args.test_sdf_input_file
