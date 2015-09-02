/*
 * =====================================================================
 *
 *        Filename:  BioactiveConstituentsManager.java
 *
 *     Description:  
 *
 *         Version:  1.0
 *         Created:  Jul 15, 2014 5:44:04 PM
 *        Revision:  none
 *
 *          Author:  James Addo (UMKC)
 *         Affiliations:  University of Missouri-Kansas City
 *
 * =====================================================================
 */

package managers;

import java.util.List;

import org.hibernate.SessionFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Repository;
import org.springframework.transaction.annotation.Transactional;

import entities.BioactiveConstituents;

@Repository
@Transactional
public class BioactiveConstituentsManager {
	
	@Autowired
	private SessionFactory sessionFactory;
	
	public BioactiveConstituentsManager() { };
	
	public List<BioactiveConstituents> getEntityByName(String keyword) {
		String query = "FROM BioactiveConstituents WHERE name LIKE '%" + keyword +"%'";
		
		return sessionFactory.getCurrentSession().createQuery(query).list();
	}
	
	public List<BioactiveConstituents> getEntityBySmiles(String keyword) {
		String query = "FROM BioactiveConstituents WHERE smiles LIKE '%" + keyword +"%' ORDER BY id ASC";
		
		return sessionFactory.getCurrentSession().createQuery(query).list();
	}
}
