# Python script to generate ML classification QSAR models using Morgan Fingerprints and two machine-learning methods: random forest and naive Bayes using python, RDKit and scikit-learn machine learning functions
#
# Example: QSAR modeling of 5-HT1A Receptor Binding: Predictions of 5-HT1A Ligands Using Machine Learning Regression/Classification Models.
# This is the code for ML classification QSAR/File contains SDF and activities information.
import cPickle, gzip, numpy, copy, math
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, rdmolops
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import rdMolDescriptors
from matplotlib import cm
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import BernoulliNB
from rdkit.Chem.Draw import IPythonConsole
from sklearn.cross_validation import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
import numpy as np
import matplotlib.pyplot as plt
import os

# prepare command-line argument parser
parser = argparse.ArgumentParser()
parser.add_argument('-i', "-- sdf_input_file", nargs='?', help="input sdf file")
parser.add_argument('-t', "-- test_sdf_input_file ", nargs='?', help="input sdf file")
parser.add_argument('-o', "--output", type=str, default=sys.stdout)

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
# read the mols
# generate dict to store molecules/mols and activities
mols_dict = {}
suppl = Chem.SDMolSupplier('sdf_input_file')
for mol in suppl:
    # mols.key = mol
    # mols_dic[k] = act_value
    act_value =  mol.GetProp('_Name')
    mols_dic = dict(zip(mol, act_value)) # mols_dic [k] = act_value             

# if file contains: [sdf, activities]
# generate outcomes
act_threshold = ""
mols_act_label = []
mols_act_class_dict = {}

for k,v in mols_dic.items():
if mols_dic[k]) < act_threshold, mols_act_class_dict[k] == 'Active' and mols_act_class[] = 1
elif
mols_dic[k]) > act_threshold, mols_act_class_dict[k] == 'Active' and mols_act_label[] = 0
# calculate fingerprints    	
fps = []
ys_fit = []
for k,v in mols_dic.items():       
    fp = AllChem.GetMorganFingerprintAsBitVect(k, 2)
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    fps.append(arr)

print "morgan2 generated:", len(fps)
    ys_fit.append(mols_act_label)
    rf = RandomForestClassifier(n_estimators=100, max_depth=10, min_samples_split=2, min_samples_leaf=1)
    rf.fit(fps, ys_fit)
    # write the model
    cPickle.dump(rf, outfile, 2)
    i += 1
outfile.close()
print "done"
# get test fps
# precalculate fingerprints for test molecules
testmols_dict = {}
suppl = Chem.SDMolSupplier('test_sdf_input_file')
for testmol in suppl:
    # testmols.key = test_mol
    # testmols_dic[k]) = test_act_value
    test_act_value =  mol.GetProp('_Name')
    test_mols_dic = dict(zip(test_mol, test_act_value)) # testmols_dic [k] = test_act_value             


test_mols_act_label = []
test_mols_act_class_dict = {}

for k,v in test_mols_dic.items():
if test_mols_dic[k]) < act_threshold, test_mols_act_class_dict[k] == 'Active' and mols_act_class[] = 1
elif
test_mols_dic[k]) > act_threshold, test_mols_act_class_dict[k] == 'Active' and test_mols_act_label[] = 0

test_fps = []
for k,v in testmols_dic.items():       
    test_fp = AllChem.GetMorganFingerprintAsBitVect(k, 2)
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(test_fp, arr)
    test_fps.append(arr)

# make predictions for test compounds
test_SVMpredictions = rf.predict(test_fps)

# output confusion matrix and percentage accuracy on test sets    
    print metrics.confusionmatrix(test_mols_act_label, test_SVMpredictions) 
    accuracy = rf.score(test_fps, test_mols_act_label)
           print accuracy 

# calculate probabilities for each test molecules
test_SVMprobabilities = rf.predict_proba(test_fps)

# compute AUC metric for this CV fold
        fpr, tpr, thresholds = metrics.roc_curve(test_mols_act_label, test_SVMprobabilities)
        roc_auc = metrics.auc(fpr, tpr)
        print "AUC (fold %d/%d): %f" % (i + 1, n, roc_auc)
        mean_auc += roc_auc

    print "Mean AUC: %f" % (mean_auc/n)
    plt.plot(fpr, tpr, label="Model#%d (AUC=%.2f)" % (model_id + 1, roc_auc))

# baseline, axes, labels, etc
plt.plot([0, 1], [0, 1], "k--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc="lower right")
plt.show()

