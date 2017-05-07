# Python script to generate ML regression QSAR models using Morgan Fingerprints and two machine-learning methods: random forest and naive Bayes using python, RDKit and scikit-learn machine learning functions
#
# Example: QSAR modeling of 5-HT1A Receptor Binding: Predictions of 5-HT1A Ligands Using Machine Learning Regression/Classification Models.
# This is the code for ML regression QSAR/From an sdf file/File contains SDF and activities information: using cross_validation.train_test_split

## 02_train_models.py
from rdkit import Chem
from scipy import interp
import numpy as np
from sklearn import cross_validation
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.cross_validation import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_score,recall_score
from sklearn import preprocessing
import os,sys,csv,cPickle,gzip
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
    DataStructs.ConvertToNumpyArray(fp, arr)
    fps.append(arr)

print "morgan2 generated:", len(fps)
    acts.append(mols_dic[k])
    acts = []
    
fps_array = np.asarray(fp)
acts_array   = np.array(acts)

fps_train,fps_test,act_train,act_test = cross_validation.train_test_split(fps_array, acts_array ,test_size=.4,random_state=randomseedcounter)

    rf = RandomForestClassifier(n_estimators=100, max_depth=10, min_samples_split=2, min_samples_leaf=1)
    rf.fit(fps_train, act_train)
    # write the model
    cPickle.dump(rf, outfile, 2)
    i += 1
outfile.close()
print "done"
# make predictions for test compounds
test_SVMpredictions = rf.predict(fps_test)

# cross-validate
scores = cross_validation.cross_val_score(rf, fps_train, act_train, cv=5)
# scores = cross_validation.cross_val_score(clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.zero_one_score)
