from . import seqviz


START1= '''<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd"><html>
<head><meta http-equiv="content-type" content="text/html; charset=utf-8" />
<title>clusters information</title>
<style type="text/css" title="currentStyle">
@import "../css/info_css.css";	
</style>
<script type="text/javascript">

function showstuff(boxid){
   document.getElementById(boxid).style.visibility="visible";
}
 
function hidestuff(boxid){
   document.getElementById(boxid).style.visibility="hidden";
}

</script>
<script src="../js/Scribl.1.1.4.min.js" ></script>
<script src="../js/jquery.min.js"></script>
<script src="../js/jquery-ui.min.js"></script>
<link href="../js/jquery-ui.css" rel="stylesheet" type="text/css"/>
<script type="text/javascript" src="../js/dragscrollable.js"></script> 
		<style>
		   #scribl-zoom-slider {
		      width: 4px;
		   }
		</style>
'''
START2='''
</head>
<body id="cluster">
'''

JSSTART= '''<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd"><html>
<head><meta http-equiv="content-type" content="text/html; charset=utf-8" />
<title>clusters information</title>
<style type="text/css" title="currentStyle">@import "css/demo_page.css";	@import "css/jquery.dataTables.css";</style><script type="text/javascript" language="javascript" src="js/jquery.js"></script>
<script type="text/javascript" language="javascript" src="js/jquery.dataTables.js"></script><script type="text/javascript" charset="utf-8">$(document).ready(function() {				$('#table').dataTable();} );</script>

</head>
<body id="dt_example">
'''

END="</body></html>"
STABLE='<table cellpadding="0" cellspacing="0" border="0" class="info" id="table" width="500px">' 
SHR="<tr>"
ETABLE="</tbody></table>"
EHR="</tr>\n"


HEADER_SEQ=["sequence","size"]
HEADER_L=["position","annotation (element:strand)-(distance to 5,distance to 3)"]

def tab2table(file,output):
	f = open(file, 'r')
	header=f.readline().strip()
	for line in f:
		line=line.strip()
		cols = line.split("\t")
	f.close()
	return 0


def make_cell(cell):
	return  "<td>%s</td>" % str(cell)

def make_hs_link(id):
	return "<a id=\"myLink\" href=\"javascript:showstuff(%s);\">show</a> <a id=\"myLink\" href=\"javascript:hidestuff(%s);\">hide</a>" % (id,id)

def make_cell_link(cell,link):
	return  "<td><a href=%s>%s</a></td>" % (link,str(cell))

def make_line(line):
	return SHR+line+EHR

def make_table(c,name):
	tmp=STABLE
	return "%s %s %s " % (tmp,c,ETABLE)

def make_div(c,name,css):
	return "<div class=%s name=%s id=%s>%s</div>" % (css,name,name,c)

def make_jshtml(c,name):
	tmp=JSSTART
	return "%s %s %s "  % (tmp,c,END )

def make_html(c,showseq,name):
	header=START1+showseq+START2
	return "%s %s %s "  % (header,c,END)

def make_header(c):
	return "<thead><tr>%s</tr></thead><tbody>" % c

def make_cell_header(c):
	return "<th>%s</th>\n" % c

def make_a(c,name):
	return "<a name=%s>%s</a>\n" % (name,c)

def make_link(c,link):
	return  "<a href=%s>%s</a>" % (link,c)