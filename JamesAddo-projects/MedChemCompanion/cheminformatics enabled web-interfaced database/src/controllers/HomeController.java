/*
 * =====================================================================
 *
 *        Filename:  HomeController.java
 *
 *     Description:  Home page controller class.
 *
 *         Version:  1.0
 *         Created:  Jun 10, 2014 4:37:34 PM
 *        Revision:  none
 *
 *          Author:  James Addo (UMKC)
 *         Affiliations:  University of Missouri-Kansas City
 *
 * =====================================================================
 */

package controllers;

import org.springframework.stereotype.Controller;
import org.springframework.ui.ModelMap;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;

@Controller
public class HomeController {
	
	@RequestMapping(value="/index.html", method=RequestMethod.GET)
	public String fetchHome(ModelMap model) {
		return "index";
	}
}
