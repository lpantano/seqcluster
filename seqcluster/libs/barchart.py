

HEAD='''
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
    
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
        <title>DB summary</title>
        <link rel="stylesheet" href="css/style.css" type="text/css">
        <script src="js/amcharts.js" type="text/javascript"></script>        
        <script type="text/javascript">
            var chart;

'''

INIFUNC='''

            AmCharts.ready(function () {
                // SERIAL CHART
                chart = new AmCharts.AmSerialChart();
                chart.dataProvider = chartData;
                chart.categoryField = "args";
                chart.plotAreaBorderAlpha = 0.2;

                // AXES
                // category
                var categoryAxis = chart.categoryAxis;
                categoryAxis.labelRotation=90;
                categoryAxis.gridAlpha = 0.1;
                categoryAxis.axisAlpha = 0;
                categoryAxis.gridPosition = "start";

                // value
                var valueAxis = new AmCharts.ValueAxis();
                valueAxis.stackType = "regular";
                valueAxis.gridAlpha = 0.1;
                valueAxis.axisAlpha = 0;
                chart.addValueAxis(valueAxis);

'''

ENDFUNC='''
              // LEGEND                  
                var legend = new AmCharts.AmLegend();
                legend.borderAlpha = 0.2;
                legend.horizontalGap = 10;
                chart.addLegend(legend);

                // WRITE
                chart.write("chartdiv");
            });

 
'''
TAIL='''
            function setDepth() {
                if (document.getElementById("rb1").checked) {
                    chart.depth3D = 0;
                    chart.angle = 0;
                } else {
                    chart.depth3D = 25;
                    chart.angle = 30;
                }
                chart.validateNow();
            }

        </script>
    </head>
    
    <body>
        <div id="chartdiv" style="width: 900px; height: 400px;"></div>
        <div style="margin-left:30px;">
	        <input type="radio" checked="true" name="group" id="rb1" onclick="setDepth()">2D
	        <input type="radio" name="group" id="rb2" onclick="setDepth()">3D
		</div>
     <div class="footer">Created by armcharts.com<div>
    </body>
    
</html>
'''

GRAPH='''
                var graph = new AmCharts.AmGraph();
                graph.title = "%(title)s";
                graph.labelText = "[[value]]";
                graph.valueField = "%(field)s";
                graph.type = "column";
                graph.lineAlpha = 0;
                graph.fillAlphas = 1;
                graph.lineColor = "%(color)s";
                chart.addGraph(graph);
'''


COL=["black","orange","red","yellow","green","blue","grey"]

def createdata(info):
	chartData="var chartData = ["
	for db in info:
		infodb=db
		#print "%s %s \n" % (infodb["DB"],infodb["uni"])
		jsdb="{"
		for k in infodb.keys():
			jsdb+='"%s":"%s",' % (k,infodb[k])
		#jsdb='{"args":%s, "uni":%s , "mul":%s, "nocon":%s\n},' % (infodb["DB"],infodb["uni"],infodb["mul"],infodb["nocon"])
		chartData+=jsdb[:-1]+"},"
	return chartData[:-1]+"];"

def addgraph(title,field,color):
	graph= GRAPH % {'title':title,'field':field,'color':color}
	return graph

def createchart(keys):
	chart=""
	i=1
	for k in keys:
		chart+=addgraph(k,k,COL[i])
		i+=1
	return chart

def createhtml(info,keys):

	data=createdata(info)
	chart=createchart(keys)

	return HEAD+data+INIFUNC+chart+ENDFUNC+TAIL

#info=[{"DB":"Repeat","uni":34,"mul":45,"nocon":3},
#{"DB":"RNA","uni":15,"mul":4,"nocon":5}]
#bar=["uni","mul","nocon"]

#print createdata(info)
#print addgraph("unique","uni","green")
#print createchart(bar)




