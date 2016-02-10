/*C++ version of MedChemCompanion MMPMaker's Molecule Fragmenter and Matched Molecular Pairs Identifier based on C++ implementation of Jameed Hussain and Ceara Rea algorithm-Computationally Efficient Algorithm to Identify Matched Molecular Pairs (MMPs) in Large Data Sets, J. Chem. Inf. Model., 2010, 50 (3), pp 339–348.
/
/*
 * Copyright (C)2015, James Addo.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided
 *   with the distribution.
 *
 * - Neither the name of James Addo
 *   nor the names of its contributors may be used to endorse or promote
 *   products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//
// Created by James Addo, November 2015
//

#include "MMPMaker.h"
#include "utilities.h"

#include <RDKitBase.h>
#include <SmilesParse/SmilesParse.h>
#include <BadFileException.h>
#include <boost/foreach.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>
#include <Exceptions.h>

using namespace RDKit;
using namespace std;
//Global smarts used by the program
String acyc_smarts = ("[*]!@!=!#[*]");
String cyc_smarts = ("[R1,R2]@[r;!R1]");

//smarts used to find appropriate fragment for
//would use SMARTS: [$([#0][r].[r][#0]),$([#0][r][#0])]
//but rdkit doesn't support component SMARTS in recursive one - $([#0][r].[r][#0])
//hence split into two
String cSma1 = ("[#0][r].[r][#0]");
String cSma2 = ("[#0][r][#0]");

//get RDKit molecule from smarts
ROMol *acyc = SmartsToMol(acyc_smarts);
ROMol *cyc = SmartsToMol(cyc_smarts);
ROMol *Sma1 = SmartsToMol(cSma1);
ROMol *Sma2 = SmartsToMol(cSma2);


//initiate and manipulate RDKit molecule
const ROMol &mol;
RWMol *newMol = static_cast<RWMol *>(new ROMol(mol));
//query mol heavy atom count
int hac = newMol->getNumAtoms();

void deletebonds(const ROMol &mol, String ftype, int hac){
RWMol *newMol = static_cast<RWMol*>(new ROMol(mol,false));
int total_acyclic = 0;
int total_cyclic = 0;
//find the relevant bonds to break

// Single acyclic Cuts
if(ftype == acyc_smarts){
acyclic_matching_atoms = mol.getSubstructMatches(acyc);
total_acyclic = acyclic_matching_atoms.size();
MatchVectType acyclic_matching_atoms;
std::vector<int> bonds_selected;
Match_Vect bonds_selected;
SubstructMatch(mol, acyc, acyclic_matching_atoms);
// if we didn't find any matches, there's nothing to be done here
    // simply return a list with a copy of the starting molecule
    if (acyclic_matching_atoms.size() == 0) {
      newMol.push_back(ROMOL_SPTR(new ROMol(mol,false)));
      newMol[0]->clearComputedProps(false);
      return newMol;
for(MatchVectType::const_iterator mvit=acyclic_matching_atoms.begin();
        mvit!=acyclic_matching_atoms.end(); mvit++){
      bonds_selected.push_back(mvit->first);
      bonds_selected.push_back(mvit->second);
	bonds_selected[0] = mvit->first;
	bonds_selected[1] = mvit->second;

Atom *at1 = newMol->getAtomWithIdx(bonds_selected[0]);
Atom *at2 = newMol->getAtomWithIdx(bonds_selected[1]);
Atom *atom0 = newMol->getAtomWithIdx(0);
Atom *atom1 = newMol->getAtomWithIdx(0);
newMol->removeBond(bonds_selected[0], bonds_selected[1]);// Break the bond with idx=at1, at2.

// Introduce two dummy atoms in the molecule
newMol->addAtom(atom0);
newMol->addAtom(atom1);
// Bond the dummy atoms to the new terminal atoms
Bond *oBond=newMol->getBondBetweenAtoms(bonds_selected[0],atom0->getIdx());
          CHECK_INVARIANT(oBond,"required bond not found");
          newMol->addBond(at1,
                          *atom0,Bond::SINGLE);

Bond *oBond=newMol->getBondBetweenAtoms(bonds_selected[1],atom1->getIdx());
          CHECK_INVARIANT(oBond,"required bond not found");
          newMol->addBond(at2,
                          *atom1,Bond::SINGLE);

break;
}

//Now get the modified fragment in smiles(i.e. smi2);

String smi2 = MolToSmiles(*newMol);

//printf (smi2);
	}
}

//cyclic Cuts
if(ftype == cyc_smarts){
cyclic_matching_atoms = mol.getSubstructMatches(cyc);
int total_cyclic = cyclic_matching_atoms.size();
MatchVectType cyclic_matching_atoms;
std::vector<int> bonds_selected;
Match_Vect bonds_selected;
SubstructMatch(mol, cyc, cyclic_matching_atoms);
// if we didn't find any matches, there's nothing to be done here
    // simply return a list with a copy of the starting molecule
    if (cyclic_matching_atoms.size() == 0) {
      newMol.push_back(ROMOL_SPTR(new ROMol(mol,false)));
      newMol[0]->clearComputedProps(false);
      return newMol;
for(MatchVectType::const_iterator mvit=cyclic_matching_atoms.begin();
        mvit!=cyclic_matching_atoms.end(); mvit++){
      bonds_selected.push_back(mvit->first);
      bonds_selected.push_back(mvit->second);
	bonds_selected[0] = mvit->first;
	bonds_selected[1] = mvit->second;

Atom *at1 = newMol->getAtomWithIdx(bonds_selected[0]);
Atom *at2 = newMol->getAtomWithIdx(bonds_selected[1]);
Atom *atom0 = newMol->getAtomWithIdx(0);
Atom *atom1 = newMol->getAtomWithIdx(0);
newMol->removeBond(bonds_selected[0], bonds_selected[1]);// Break the bond with idx=at1, at2.

// Introduce two dummy atoms in the molecule
newMol->addAtom(atom0);
newMol->addAtom(atom1);
// Bond the dummy atoms to the new terminal atoms
Bond *oBond=newMol->getBondBetweenAtoms(bonds_selected[0],atom0->getIdx());
          CHECK_INVARIANT(oBond,"required bond not found");
          newMol->addBond(at1,
                          *atom0,Bond::SINGLE);

Bond *oBond=newMol->getBondBetweenAtoms(bonds_selected[1],atom1->getIdx());
          CHECK_INVARIANT(oBond,"required bond not found");
          newMol->addBond(at2,
                          *atom1,Bond::SINGLE);

continue;
}
//Now get the modified fragment in smiles(i.e. smi2);

String smi2 = MolToSmiles(*newMol);
//printf (smi2);

//now do an acyclic cut with the successful cyclic cut on the mol
for(MatchVectType::const_iterator mvit=acyclic_matching_atoms.begin();
        mvit!=acyclic_matching_atoms.end(); mvit++){
      bonds_selected.push_back(mvit->first);
      bonds_selected.push_back(mvit->second);
	bonds_selected[0] = mvit->first;
	bonds_selected[1] = mvit->second;

Atom *at1 = newMol->getAtomWithIdx(bonds_selected[0]);
Atom *at2 = newMol->getAtomWithIdx(bonds_selected[1]);
Atom *atom0 = newMol->getAtomWithIdx(0);
Atom *atom1 = newMol->getAtomWithIdx(0);
newMol->removeBond(bonds_selected[0], bonds_selected[1]);// Break the bond with idx=at1, at2.

// Introduce two dummy atoms in the molecule
newMol->addAtom(atom0);
newMol->addAtom(atom1);
// Bond the dummy atoms to the new terminal atoms
Bond *oBond=newMol->getBondBetweenAtoms(bonds_selected[0],atom0->getIdx());
          CHECK_INVARIANT(oBond,"required bond not found");
          newMol->addBond(at1,
                          *nbrIdx,oBond->getBondType());

Bond *oBond=newMol->getBondBetweenAtoms(bonds_selected[1],atom1->getIdx());
          CHECK_INVARIANT(oBond,"required bond not found");
          newMol->addBond(at2,
                          *atom1,Bond::SINGLE);

continue;
}

//Now get the modified fragment in smiles(i.e. smi2);

String smi2 = MolToSmiles(*newMol);
//printf (smi2);
	}
}

//determine whether ring cut is valid
String cSma1 = ("[#0][r].[r][#0]");
static ROMol *Sma1 = SmartsToMol("[#0][r].[r][#0]");
String cSma2 = ("[#0][r][#0]");
static ROMol *Sma2 = SmartsToMol("[#0][r][#0]");
void is_ring_cut_valid(ROMol *fMol, ROMol *Sma1, ROMol *Sma2){
//to check is a fragment is a valid ring cut, it needs to match the
//smarts: [$([#0][r].[r][#0]),$([#0][r][#0])]

	boolean valid = false;
	ROMol *m = new RWMol();
//if m is not None:
	if (m != NULL){
    //use global smarts
    if(m->hasSubstructMatch(Sma1) || m->hasSubstructMatch(Sma2)){
    	int atom_count = m->getNumAtoms();
    	valid = true;
    		}
	}
}

//determine number of attachmentPoints given a fragmented molecule,
void attachmentPoints(String smi2){
//Creates an Integer array for storing count results
	int attachments;
	int asterisks = 0;
	boolean match = false;

	//loop //if (smi2[i] == '*')
	for (int i = 0; i < smi2.length(); i ++) {
	    if (smi2.at(i) == '*') {
	        asterisks++;
	        match = true;
	        attachments = asterisks++;
	        break;
	   }
    }
}

void select_fragments(String smi2, int attachments,String ftype,int hac){


	int result_hcount = 0;
      std::vector<ROMOL_SPTR> &molfrag;
      std::vector<std::string> molfragsmi;
      vector<string> molfragsmi = split(smi2, '.');
	if(ftype == "acyclic"){
//count each split fragments number of attachments ("*")
//for each fragment
for (int i = 0; i < molfragsmi.size(); i ++) {
        attachments = attachmentPoints(molfragsmi[i]);
        	// check if terminal fragment
            if(attachments == 1){
                //fMol = Chem.MolFromSmiles(f)
		    RMol *fMol=SmilesToMol(molfragsmi[i]);
                fhac = fMol->GetNumAtoms();
	//check fragment heavy atom count
    	//if the fragment is 2 atoms (or less - includes attachment) it 	  is too small
    	//to be interesting. This check has the additional benefit
    	//of pulling out the relevant single cuts as it discards
    	//fragments where we only chop off a small part of the input cmpd
    	//iterate through fragments

    	if(fhac > 3){
    	String result = molfragsmi[i];
    	result_hcount = result_hcount + fhac;   	
    	//needs to be greater than 60% of parent mol
    if ((result != NULL) && (result_hcount > 0.6*hac)){
    		//remove first character as it will always be "."	
    		 result = result.substr(1);
    		 result = molfragsmi[i];
    	//to output the smiles
	vector<string>results;
	results.push_back(molfragsmi[i]);
    	molfrag.push_back(fMol);
		  }   	
        }
    }
    	else {
    	molfrag = NULL;
    	vector<string>results = NULL;
    		 }   	  	
       }
  }

   
else if(ftype == "cyclic"){
molfrag = NULL;
int result_hcount = 0;
//make sure it is 2 components and valid ring cut
	for (int i = 0; i < molfragsmi.size; i++){
	 	if(molfragsmi.size == 2){
       		//check if a valid cut
	 		ROMol *fMol=SmilesToMol(molfragsmi[i]);
			int atom_count = 0;
			Boolean valid = true;
			//check if a valid cut
	 		valid, result_hcount = is_ring_cut_valid(fMol, Sma1, Sma2);
		      //int hac = fMol->getNumAtoms();
	 		//boolean valid = true;
	 		if(valid){
				//needs to be greater 3 heavy atoms and greater than 40% of parent mol
		     		if((result_hcount  > 3) && (result_hcount  > 0.4	*hac)){
		     			//to output the smiles
	vector<string>results;
	results.push_back(molfragsmi[i]);
    	molfrag.push_back(fMol);
				}
			}
	 	}
	 		else{
	 		    	molfrag = NULL;
	 		    	vector<string>results = NULL;		     	 	
	 	}
	}
}
           	   
else if(ftype == "cyclic_and_acyclic"){
	//printf (f_smi);
	molfrag  = Null;
	int result_hcount = 0;
	String result = NULL;
	boolean valid = true;
	for (int i = 0; i < molfragsmi.size; i++){
		attachments = attachmentPoints(molfragsmi[i]);
                ROMol *fMol=SmilesToMol(molfragsmi[i]);
				
				//int atom_count = fMol->getNumAtoms();			
				//attachments = f.count("*")
				//for (int i = 0; i < molfragsmi.size(); i++) {
        			//attachments = attachmentPoints(molfragsmi[i])
	//need to find the fragments which are valid which means they must be:
				//  Terminal (one attachment point) or valid ring cut

            if(attachments >= 3){
                continue;
}
            if(attachments == 2){
            	//check if a valid cut
            	valid,result_hcount = is_ring_cut_valid	(fMol,Sma1, Sma2);
            	// needs to be greater 3 heavy atoms
            	if(valid){
                    if(result_hcount > 3){
                       String result = molfragsmi[i];
                       vector<string>results;
                       	results.push_back(molfragsmi[i]);
                           	molfrag.push_back(fMol);
                       				}
                       			}
                       	 	}
                    }
            }
else if(attachments == 1){
                //needs to be greater 3 heavy atoms
	int fhac = fMol->getNumAtoms();
	fhac = fMol->GetNumAtoms();
	//needs to be greater 3 heavy atoms
                if(fhac > 3){
                    String result = molfragsmi[i];
                    vector<string>results;
                    results.push_back(molfragsmi[i]);
                    molfrag.push_back(fMol);
                    result_hcount = result_hcount + fhac;

                    		//printf(result)
                }
		   }

		}
	}

            //appropriate fragmentations must have 2 components
            //result will always start with . because of the way it is constructed
            //hence 2 component result wil contain 2 dots

int count_dot(string molfragsmi[i]){
int count = 0;
for (int i = 0; i < molfragsmi[i].size(); i++){
if ((molfragsmi[i]) == '.'){
count++;
	}
return count;
}
if( (result != NULL) && (count_dot(result) == 2)) {
	//take off the starting dot when building smiles
ROMol *fMol=SmilesToMol(result.substr(1));

				//needs to be greater 3 heavy atoms and greater than 60% of parent mol
if((result_hcount  > 3) && (result_hcount  > 0.6*hac)){
//take off the starting dot
                result = result.substr(1);
                String result = molfragsmi[i];
            	std::vector<string>results;
            	results.push_back(molfragsmi[i]);
                	molfrag.push_back(fMol);
}

            else {
                result = NULL;
	}

}
    return result;
		}

//needs to be greater 3 heavy atoms //to output the smiles
//if(result_hcount > 3){
//if(fhac > 3){
//if(result_hcount > 3){
		
		
	
	//result_hcount = result_hcount + fhac;
			
		
	

//generate_molecule_fragments
void generate_molecule_fragments (RWMol *mol){
int hac = (Integer) NULL;
String ftype = NULL;
	get_bonds_to_break(mol);
deletebonds(mol, ftype,hac);
//int attachments = 0;
//String smi2 = NULL;
String smi2 = MolToSmiles(*mol);
int attachments = attachmentPoints(smi2);
select_fragments(smi2, attachments,ftype,hac);
}

int main (int argc, char* argv[]){

	  if (!RDKit.activate()) {
	      throw new UnsatisfiedLinkError("RDKit library could not be loaded.");
	    }
	String Fragtype = NULL;
		// Check the number of parameters
    if (argc >=  2) {
			String filename1 = args[0];
			Fragtype = args[2];}
    
    {
     	// create new file
  SDMolSupplier msuppl(filename1);
  BOOST_LOG(rdInfoLog)<<"loading mols: "<<std::endl;
  while(!msuppl.atEnd()){
  ROMol *mol=msuppl.next();
  RWMol *currentmol = new RWMol(mol);
  generate_molecule_fragments(currentmol);
  		}
      }
  }

//std::vector<ROMOL_SPTR> &mols
//generate_MMPs
void generate_MMPs (std::vector<ROMOL_SPTR> molfrag, vector<string>results){

  		if (argc >= 2) {
  			String filename2 = args[1];
  		// create new file
		//std::ifstream MMPmolFile(filename2.c_str());
		std::ifstream MMPmolFile;
		MMPmolFile.open("filename2");
		if (!MMPmolFile || (MMPmolFile.bad()) ) {
       	std::ostringstream errout;
        	errout << "Bad input file " << filename2;
         	throw BadFileException(errout.str());
  	 	SDMolSupplier *suppl = new SDMolSupplier(filename2);
		BOOST_LOG(rdInfoLog)<<"loading mols: "<<std::endl;
  	 			while (!suppl.atEnd()){
  	  				RWMol *MMPmol = suppl.next();
  	  				RWMol *currentMMPmol = new RWMol(MMPmol);
  	  			//* create a substructure search pattern object from molecule fragments*/
  	    for (int i = 0; i < molfrag.size(); i++){
  	    ROMol *ss = SmilesToMol("molfrag[i]");
  	    ROMol *ssh = SmilesToMol("results[i]");
  	        if (currentMMPmol.hasSubstructMatch(ss));
		std::ofstream fw;
		String filename3 = args[2];
		fw.open("filename3");
  		try {
  			fw.write(MolToSmiles(*MMPmol));
  		} catch (IOException& e) {
  			// TODO Auto-generated catch block
  			e.printStackTrace();
  		}
  		fw.close();

      	   }
     	     }
        }
     }
}
