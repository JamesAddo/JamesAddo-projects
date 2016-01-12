/*MedChemCompanion MMPMaker's Molecule Fragmenter and Matched Molecular Pairs Identifier based on Java implementation of Jameed Hussain and Ceara Rea algorithm-Computationally Efficient Algorithm to Identify Matched Molecular Pairs (MMPs) in Large Data Sets, J. Chem. Inf. Model., 2010, 50 (3), pp 339–348.
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

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.io.*;

import org.RDKit.Atom;
import org.RDKit.Bond;
import org.RDKit.Match_Vect;
import org.RDKit.Match_Vect_Vect;
import org.RDKit.ROMol;
import org.RDKit.RWMol;
import org.RDKit.SDMolSupplier;
import org.RDKit.MolSupplier;
import org.RDKit.Bond_Vect;
import org.RDKit.ROMol_Vect;
import org.RDKit.ROMol_Vect_Vect;
import org.RDKit.Atom_Vect;
public class MMPMaker{
   private static final RWMol RWMol = null;
	//get_bonds_to_break
   public static void get_bonds_to_break(RWMol mol){
		  //Global smarts used by the program
		  //acyclic bond smarts 
	   // String acyc_smarts = ("[*]!@!=!#[*]");
		// ROMol acyc;
		// acyc = RWMol.MolFromSmarts("[*]!@!=!#[*]");
	
		//exocyclic/fused exocyclic bond smarts
		// String cyc_smarts = ("[R1,R2]@[r;!R1]");
		// ROMol cyc;
		//cyc = RWMol.MolFromSmarts("[R1,R2]@[r;!R1]");
		  //smarts used to find appropriate fragment for
		  //would use SMARTS: [$([#0][r].[r][#0]),$([#0][r][#0])]
		  //but rdkit doesn't support component SMARTS in recursive one - $([#0][r].[r][#0])
		  //hence split into two
		// String cSma1 = ("[#0][r].[r][#0]");
		  ROMol Sma1;
		  Sma1 = RWMol.MolFromSmarts("[#0][r].[r][#0]");
		  // String cSma2 = ("[#0][r][#0]");
		  ROMol Sma2;
		  Sma2 = RWMol.MolFromSmarts("[#0][r][#0]");
	      //query mol heavy atom count
	      int hac = (int) mol.getNumAtoms();
	      //different cuts can give the same fragments
	      //to use out_fragments to remove them
	      
	      // Single acyclic Cuts
	  //find the relevant bonds to break
	      //Match_Vect_Vect acyclic_matching_atoms = mol.getSubstructMatches(acyc);
	  //System.out.println (Matching Atoms:);
	  //System.out.println (acyclic matching atoms: ",acyclic_matching_atoms);
	      // return acyclic_matching_atoms ;
	  	
}

protected static void deletebonds(RWMol mol, String ftype, int hac){
			String acyc_smarts = ("[*]!@!=!#[*]");
			ROMol acyc;
			acyc = RWMol.MolFromSmarts("[*]!@!=!#[*]");
			String cyc_smarts = ("[R1,R2]@[r;!R1]");
			ROMol cyc;
			cyc = RWMol.MolFromSmarts("[R1,R2]@[r;!R1]");
			Match_Vect_Vect acyclic_matching_atoms = null;
			Match_Vect_Vect cyclic_matching_atoms = mol.getSubstructMatches(cyc);
			int total_acyclic = 0;
		if(ftype == acyc_smarts){	
			acyclic_matching_atoms = mol.getSubstructMatches(acyc);
			total_acyclic = (int) acyclic_matching_atoms.size(); 
				for (int i = 0; i < total_acyclic; i++){
			  		Match_Vect bonds_selected = acyclic_matching_atoms.get(i); //define bonds_selected
			  			int bond_to_break_atom1 = bonds_selected.get(0).getFirst();
			  			int bond_to_break_atom2 = bonds_selected.get(0).getSecond();
			  		  Atom at1 = mol.getAtomWithIdx(bond_to_break_atom1);
					  Atom at2 = mol.getAtomWithIdx(bond_to_break_atom2);
					  Atom atom0 = mol.getAtomWithIdx(0);
					  Atom atom1 = mol.getAtomWithIdx(0);
				mol.removeBond(bond_to_break_atom1, bond_to_break_atom2);// Break the bond with idx=at1, at2. 
				// Introduce two dummy atoms in the molecule
				mol.addAtom(atom0);
				mol.addAtom(atom1);
			    // Bond the dummy atoms to the new terminal atoms	
				mol.addBond(at1, atom0, Bond.BondType.SINGLE); 
				mol.addBond(at2, atom1, Bond.BondType.SINGLE);
			    continue;
		}
				//Now get the modified fragment in smiles(i.e. smi2);
			    RWMol newmol = new RWMol(mol);
			    String smi2 = ((RWMol) newmol).MolToSmiles();
			    //System.out.println (smi2);
   }		   

		//Fused/Spiro exocyclic bond Cuts
		//find the relevant cyclic bonds to break
		else if(ftype == cyc_smarts){
		cyclic_matching_atoms = mol.getSubstructMatches(cyc);
		int total_cyclic = (int) cyclic_matching_atoms.size(); 
		for (int i = 0; i < total_acyclic; i++){
				Match_Vect bonds_selected = cyclic_matching_atoms.get(i); //define bonds_selected
					int bond_to_break_atom1 = bonds_selected.get(0).getFirst();
					int bond_to_break_atom2 = bonds_selected.get(0).getSecond();
				  Atom at1 = mol.getAtomWithIdx(bond_to_break_atom1);
			  Atom at2 = mol.getAtomWithIdx(bond_to_break_atom2);
			  Atom atom0 = mol.getAtomWithIdx(0);
			  Atom atom1 = mol.getAtomWithIdx(0);
		mol.removeBond(bond_to_break_atom1, bond_to_break_atom2);// Break the bond with idx=at1, at2. 
		// Introduce two dummy atoms in the molecule
		mol.addAtom(atom0);
		mol.addAtom(atom1);
		// Bond the dummy atoms to the new terminal atoms	
		mol.addBond(at1, atom0, Bond.BondType.SINGLE); 
		mol.addBond(at2, atom1, Bond.BondType.SINGLE);
		continue;
		}
		//now do an acyclic cut with the successful cyclic cut on the mol
		for (int i = 0; i < total_acyclic; i++){
	  		Match_Vect bonds_selected = acyclic_matching_atoms.get(i); //define bonds_selected
	  			int bond_to_break_atom1 = bonds_selected.get(0).getFirst();
	  			int bond_to_break_atom2 = bonds_selected.get(0).getSecond();
	  		  Atom at1 = mol.getAtomWithIdx(bond_to_break_atom1);
			  Atom at2 = mol.getAtomWithIdx(bond_to_break_atom2);
			  Atom atom0 = mol.getAtomWithIdx(0);
			  Atom atom1 = mol.getAtomWithIdx(0);
		mol.removeBond(bond_to_break_atom1, bond_to_break_atom2);// Break the bond with idx=at1, at2. 
		// Introduce two dummy atoms in the molecule
		mol.addAtom(atom0);
		mol.addAtom(atom1);
	    // Bond the dummy atoms to the new terminal atoms	
		mol.addBond(at1, atom0, Bond.BondType.SINGLE); 
		mol.addBond(at2, atom1, Bond.BondType.SINGLE);
	    continue;
		}
		//Now get the modified fragment in smiles(i.e. smi2);
		RWMol newmol = new RWMol(mol);
		String smi2 = ((RWMol) newmol).MolToSmiles();
	}
}
   
public void testBasicInstantiation_ConstBondIterator() {
ROMol m = RWMol.MolFromSmiles("CS");
	}
	//Chem.SanitizeMol(newmol)
	//Use itertools.combinations at_pairs, N) for N-cuts
//determine whether ring cut is valid
String cSma1 = ("[#0][r].[r][#0]");
static ROMol Sma1 = org.RDKit.RWMol.MolFromSmarts("[#0][r].[r][#0]");
String cSma2 = ("[#0][r][#0]");
static ROMol Sma2 = RWMol.MolFromSmarts("[#0][r][#0]");
protected static void is_ring_cut_valid(ROMol fMol, ROMol Sma1, ROMol Sma2){
//to check is a fragment is a valid ring cut, it needs to match the
//smarts: [$([#0][r].[r][#0]),$([#0][r][#0])]
	
	boolean valid = false;
	ROMol m = new RWMol();
//if m is not None:
	if (m != null){
    //use global smarts
    if(m.hasSubstructMatch(Sma1) || m.hasSubstructMatch(Sma2)){
    	int atom_count = (int) m.getNumAtoms();
    	valid = true;
    	}  
	}
}
//determine number of attachmentPoints given a fragmented molecule, 
static void attachmentPoints(String smi2){
//Creates an Integer array for storing count results 
	int attachments; 
	int asterisks = 0;
	//loop
	for (int i = 0; i < smi2.length() - 1; i ++) {
	    if (smi2.charAt(i) == '*') {
	        asterisks++;
	        boolean valid = smi2.matches("*");
	        attachments = asterisks++;
	        break;
	   }	    

    }
}

protected static void select_fragments(String smi2, int attachments, StringBuffer f_smi,String ftype,int hac){
    int result_hcount = 0;
    ROMol_Vect molfrag = new ROMol_Vect(attachments);
    String[] parts = smi2.split("*");
    if(ftype == "acyclic"){
    for (int i = 0; i < attachments; i ++) {
    	attachmentPoints(parts[i]);   
    }
    for (int i = 0; i < attachments; i ++) {
    	//check if terminal fragment
    	if(attachments == 1){
    	ROMol fMol = RWMol.MolFromSmiles(parts[i]);
    	//check fragment heavy atom count
    	//if the fragment is 2 atoms (or less - includes attachement) it is too small
    	//to be interesting. This check has the additional benefit
    	//of pulling out the relevant single cuts as it discards
    	//fragments where we only chop off a small part of the input cmpd
    	//iterate through fragments
    	int fhac = (int) fMol.getNumAtoms();
    	if(fhac > 3){
    	result_hcount = result_hcount + fhac;
    	//needs to be greater than 60% of parent mol
    	if(result_hcount > 0.6*hac){
    	//to output the smiles instead remove first character as it will always be "."
    	StringBuilder sb = new StringBuilder(parts[i]);
	    sb.deleteCharAt(0);
	    String resultString = sb.toString();
	    ArrayList<String>results = new ArrayList<String>();
	    results.add(resultString);
    	molfrag.add(fMol);
    	}
    	else{
    	molfrag = null;
    	ArrayList<String>results = null;
    	}    	
      }
    }
  }
 }

else if(ftype == "cyclic"){
molfrag = null;
//make sure it is 2 components
	for (int i = 0; i < parts.length; i++){
	 	if(parts.length == 2){
       		//check if a valid cut
	 		ROMol fMol = RWMol.MolFromSmiles(parts[i]);
	 		is_ring_cut_valid(fMol, Sma1, Sma2);
	 		int atom_count = (int) fMol.getNumAtoms();
	 		boolean valid = true;
	 		result_hcount = atom_count; 
	 		if(valid){
				//needs to be greater 3 heavy atoms and greater than 40% of parent mol
		     		if((atom_count > 3) && (atom_count > 0.4	*hac)){
		     			StringBuilder sb = new StringBuilder(parts[i]);
		     		    sb.deleteCharAt(0);
		     		    String resultString = sb.toString();
		     		    ArrayList<String>results = new ArrayList<String>();
		     		    results.add(resultString);
		     	    	molfrag.add(fMol);
		     	    	}
		     	    	else{
		     	    	molfrag = null;
		     	    	ArrayList<String>results = null;
		     	    	}    	
				}
			}
		}
	}

else if(ftype == "cyclic_and_acyclic"){
	//System.out.println (f_smi);
	molfrag  = null;
	//need to find the fragments which are valid which means they must be:
	//  Terminal (one attachment point) or valid ring cut
		for (int i = 0; i < parts.length; i++) {
			attachmentPoints(parts[i]);
				if(attachments >= 3){
					continue;
				}
				if(attachments == 2){
		 			//check if a valid cut
				ROMol fMol = RWMol.MolFromSmiles(parts[i]);
				is_ring_cut_valid(fMol,Sma1, Sma2);
				int atom_count = (int) fMol.getNumAtoms();
		 		boolean valid = true;
						result_hcount = atom_count; 
						if(valid){
							//needs to be greater 3 heavy atoms
							if(result_hcount > 3){																
							molfrag.add(fMol);
							}
						}
					}
				else if(attachments == 1){
				
					result_hcount = 0;
					ROMol fMol = RWMol.MolFromSmiles(parts[i]);
					int fhac = (int) fMol.getNumAtoms();
					//needs to be greater 3 heavy atoms
					if(fhac > 3){
					molfrag.add(fMol);
					result_hcount = result_hcount + fhac;
					}
					
						
					}
				}
		}
	}

//generate_molecule_fragments
protected static void generate_molecule_fragments (RWMol mol){
int hac = (Integer) null;
String ftype = null;
	get_bonds_to_break(mol);				
deletebonds(mol, ftype,hac);
	StringBuffer StringBuffer;
int attachments = 0;
String smi2 = null;
StringBuffer f_smi = null;
select_fragments(smi2, attachments,f_smi,ftype,hac);
}

public static void main (String[] args) throws Exception {
	
	  if (!RDKit.activate()) {
	      throw new UnsatisfiedLinkError("RDKit library could not be loaded.");
	    }
	String Fragtype = null;
		if (args.length >= 2) {
			String filename1 = args[0];
			Fragtype = args[2];}
    boolean bool = false;
    
    {
     	// create new file
     	File molFile = new File("filename1");
     	SDMolSupplier suppl = new SDMolSupplier(molFile.getPath());
     	// true if the file path exists
     	bool = molFile.exists();
     	// if file exists
     	if(bool){
       	ROMol mol;
for (int i = 0; i < 10; i++){
     			while (!suppl.atEnd()){
      				mol = suppl.next();
      				RWMol currentmol = new RWMol(mol);    				
  	generate_molecule_fragments(currentmol);
     			}
            }
     	 }
      }

   }
  //generate_MMPs
  	protected static void generate_MMPs (ROMol_Vect molfrag, Match_Vect results) throws IOException{
  		String[] args = null;
  		if (args.length >= 2) {
  			String filename2 = args[1];
  			boolean bool = false;
  		// create new file
  	 	File MMPmolFile = new File("filename2");
  	 	SDMolSupplier suppl = new SDMolSupplier(MMPmolFile.getPath());
  	 	// true if the file path exists
  	 	bool = MMPmolFile.exists();
  	 	// if file exists
  	 	if(bool){
  	   	ROMol MMPmol;
  	 			while (!suppl.atEnd()){
  	  				MMPmol = suppl.next();
  	  				RWMol currentMMPmol = new RWMol(MMPmol);
  	  			//* create a substructure search pattern object from molecule fragments*/
  	    for (int i = 0; i < molfrag.size(); i++){
  	    ROMol ss = RWMol.MolFromSmiles("molfrag[i]");
  	    ROMol ssh = RWMol.MolFromSmiles("results[i]");
  	        if (currentMMPmol.hasSubstructMatch(ss));
  		File f = File.createTempFile("smi", ".csv");
  		FileWriter fw = new FileWriter(f);
  		try {
  			fw.write(MMPmol.MolToSmiles());
  		} catch (IOException e) {
  			// TODO Auto-generated catch block
  			e.printStackTrace();
  		}
  		fw.close();
  	
     					
     			}
      		}
     	}
     }

   }

}
