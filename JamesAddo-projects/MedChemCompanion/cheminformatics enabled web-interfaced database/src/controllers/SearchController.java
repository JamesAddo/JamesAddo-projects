/*
 * =====================================================================
 *
 *        Filename:  SearchController.java
 *
 *     Description:  
 *
 *         Version:  1.0
 *         Created:  Jun 24, 2014 3:35:25 PM
 *        Revision:  none
 *
 *          Author:  James Addo(UMKC)
 *         Affiliations:  University of Missouri-Kansas City
 *
 * =====================================================================
 */

package controllers;


import javax.servlet.http.HttpServletRequest;

import managers.BioactiveConstituentsManager;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Controller;
import org.springframework.ui.ModelMap;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RequestParam;
import org.springframework.web.bind.annotation.ResponseBody;

@Controller
public class SearchController {
	
	@Autowired
	private BioactiveConstituentsManager bioactiveConstituentsManager;
	
	@RequestMapping(value="/search.html", method=RequestMethod.GET)
	public String fetchSearch(ModelMap modelJ) {
		return "search";
	}
	
	@RequestMapping(value="/bioactive_search.do", method=RequestMethod.POST)
	public @ResponseBody String handleBioactiveSearch(@RequestParam("keyword") String keyword, HttpServletRequest request) {
		return String.valueOf(bioactiveConstituentsManager.getEntityBySmiles(keyword).size());
	}
	
	@RequestMapping(value="/medicinal_search.do", method=RequestMethod.POST)
	public @ResponseBody String handleMedicinalSearch(@RequestParam("keyword") String keyword, HttpServletRequest request) {
		
		return "Implement Medicinal search! Please";
	}
}
