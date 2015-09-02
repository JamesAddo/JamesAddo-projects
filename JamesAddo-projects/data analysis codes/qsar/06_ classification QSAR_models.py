# Python script to generate ML classification QSAR models using Morgan Fingerprints and two machine-learning methods: random forest and naive Bayes using python, RDKit and scikit-learn machine learning functions
#
# Example: QSAR modeling of 5-HT1A Receptor Binding: Predictions of 5-HT1A Ligands Using Machine Learning Regression/Classification Models.
# This is the code for ML classification QSAR/File contains SMILES and active/inactive information
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

# read in command line options
(options, args) = parser.parse_args()
# required arguments
if options.assay:
    aid = options.assay
else:
    raise RuntimeError('test assay missing')
print "RF training for test assay:", aid

mols = {}
with open(path+'test_assays/'+aid+'.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        # contains: [SID, CID, outcome, score, url, comment, value, qc]
        if row[2] == 'Active': 
            mols[row[1]] = 1
        elif row[2] == 'Inactive': 
            mols[row[1]] = 0
num = len(mols)
print "number of molecules:", num

# if file contains: [smiles, activities]
# generate outcomes
threshold = ""
for line in sub_strct:
if line[1] < act_threshold, row[1] == 'Active' and dic[line[0]] = 1
elif
if line[1] > act_threshold, row[1] == 'Inactive' and dic[line[0]] = 0
act = dic[line[0]
# generate dict to store molecules/smarts and activities and calculate fingerprints
inf = open("pains.txt", "r")
sub_strct = [ line.rstrip().split(" ") for line in inf ]# contains: [smiles, activities]
smiles = [ line[0] for line in sub_strct]
act = [ line[1] for line in sub_strct] # act = 1 or act = 0
dic = dict(zip(smiles, act)) # dic[line[0]] = float(line[1])
        
    for k,v in dic.items():
        mols = Chem.MolFromSmarts( k )
fps = []
ys_fit = []
for mol in mols:
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    fps.append(arr)

print "morgan2 generated:", len(fps)

            ys_fit.append(dic[k])
    rf = RandomForestClassifier(n_estimators=100, max_depth=10, min_samples_split=2, min_samples_leaf=1)
    rf.fit(fps, ys_fit)
    # write the model
    cPickle.dump(rf, outfile, 2)
    i += 1
outfile.close()
print "done"
