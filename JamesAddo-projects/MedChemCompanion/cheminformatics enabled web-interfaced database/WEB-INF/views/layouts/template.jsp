<%@ page language="java" contentType="text/html; charset=ISO-8859-1" pageEncoding="ISO-8859-1"%>
<%@ taglib uri="http://tiles.apache.org/tags-tiles" prefix="tiles"%>

<!DOCTYPE html>
<html>
	<head>
		<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1">
	
		<title>Advanced Biomedical</title>

		<link rel="stylesheet" type="text/css" href="assets/css/style.css"/>
		
		<script type="text/javascript" src="assets/js/plugins/jquery-1.7.2.min.js"></script>
	</head>

	<body>
		<div id="hdr" class="ctr_cls">
			<a href="index.html">
				<img alt="header image" src="assets/images/advancedbiomedical_header.jpg"/>
			</a>
		</div>
	
		<div id="ctr" class="ctr_cls">
			<div id="ctn">
				<div style= "text-align: left; margin: 15px 218px 0px 230px; text-align: justify;">
					<tiles:insertAttribute name="content"/>
				</div>
			</div>
		
			<div id="lpnl">
				<ul id="navl">
					<li>
						<img src="assets/images/li_icon.png" style="padding-right: 1em;"><a href="index.html">Home</a>
					</li>
					<li>
						<img src="assets/images/li_icon.png" style="padding-right: 1em;"><a href="#">About Us</a>
					</li>
					<li>
						<img src="assets/images/li_icon.png" style="padding-right: 1em;"><a href="#">Browse</a>
					</li>
					<li>
						<img src="assets/images/li_icon.png" style="padding-right: 1em;"><a href="search.html">Search</a>
					</li>
					<li>
						<img src="assets/images/li_icon.png" style="padding-right: 1em;"><a href="#">Sign up</a>
					</li>
					<li>
						<img src="assets/images/li_icon.png" style="padding-right: 1em;"><a href="#">Login</a>
					</li>
				</ul>
			
			</div>

			<div id="rpnl">
				<!-- Products' information -->
			</div>
		</div>
	
		<div id="ftr" class="ctr_cls">
			<!-- Footer -->
		</div>
	</body>
</html>