<script type="text/javascript">
	$(function() {
		$('#search_form').submit(function(){
			$.ajax({
				url: $(':radio[name=searchOpt]:checked').val() + 'search.do',
				type: 'POST',
				data: {keyword: $('#search_box').val()},
				success: function(data) {
					if(console) {
						console.log(data);
					}
				}
			})
			
			return false;
		});
		
		$('#mps').click(function() {
			$('#search_form span.hint').text('Enter medicinal plant name or name segments');
		})
		
		$('#bcs').click(function() {
			$('#search_form span.hint').text('Enter compound name or name segments, CAS Registry Number, or Smiles');
		})
		
	});
</script>

<label style="display: block; font-weight: bold;">Search Center:</label>

<form id="search_form">
	<ul style="list-style: none; padding: 0px;">
		<li>
			<input id="mps" type="radio" name="searchOpt" value="medicinal_"/>
			<label for="mps">Medicinal Plants</label>
		</li>
		<li>
			<input id="bcs" type="radio" name="searchOpt" value="bioactive_"/>
			<label for="bcs">Bioactive Constituents</label>
		</li>
	</ul>
	
	<input type="text" id="search_box"/> <input type="submit" value="Search"/> </br>
	<span class="hint">Select a catalog search</span>
</form>
