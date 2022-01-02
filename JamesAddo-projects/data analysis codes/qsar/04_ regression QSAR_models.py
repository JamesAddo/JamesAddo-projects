# Python script to generate ML regression QSAR models using Morgan Fingerprints and two machine-learning methods: random forest and naive Bayes using python, RDKit and scikit-learn machine learning functions
#
# Example: QSAR modeling of 5-HT1A Receptor Binding: Predictions of 5-HT1A Ligands Using Machine Learning Regression/Classification Models.
# This is the code for ML regression QSAR/File contains SMILES and activities information using an all dictionary data structure.
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

# generate dict to store molecules/smarts and activities and calculate fingerprints

inf = open("pains.txt", "r")
sub_strct = [ line.rstrip().split(" ") for line in inf ]# contains: [smiles, activities]
smiles = [ line[0] for line in sub_strct]
act = [ line[1] for line in sub_strct]
mols_dic = dict(zip(smiles, act)) # mols_dic [line[0]] = float(line[1])
        
    for k,v in dic.items():
        mols = Chem.MolFromSmiles( k )
fps = []
ys_fit = []
morgan2_fps = {}
for k in mols.keys():          
m = Chem.MolFromSmiles(k)
        fp = AllChem.GetMorganFingerprint(m, 2)
	   morgan2_fps[k] = fp
	   if morgan2_fps.has_key(k):
        fps.append(morgan2_fps[k])
print "morgan2 generated:", len(fps)
# transform fps to scipy/sklearn compatible sparse format
fps_scipy = [dict([(str(k),v) for (k,v) in fp.GetNonzeroElements().items()]) for fp in fps]
dv = DictVectorizer().fit(fps_scipy)
            ys_fit.append((mols_dic[k])
    rf = RandomForestRegressor(n_estimators=100, max_depth=10, min_samples_split=2, min_samples_leaf=1)
    rf.fit(fps, ys_fit)
    # write the model
    cPickle.dump(rf, outfile, 2)
    i += 1
outfile.close()
print "done"




