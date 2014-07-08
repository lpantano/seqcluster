JSVIZ1='''<script type="text/javascript">
		  function go() {
            var zoomCanvas = document.getElementById('canvas');
            origZoomChart = new Scribl(zoomCanvas, 100);
            //origZoomChart.scale.min = 0;
	   // origZoomChart.scale.max = 12000;
'''

JSVIZ2='''
            origZoomChart.scrollable = true;
            origZoomChart.scrollValues = [10, 250];
            origZoomChart.draw();
		  }

		</script>
'''


CANVAS='''		<div id="container">
		            <canvas id="canvas" width="940px" height="400px"  style="margin-left:auto; margin-right:auto"></canvas>  
		</div>		
'''

SEQ='''origZoomChart.addFeature( new Seq('human', %s, %s, "%s") );'''

def addseq(pos,len,seq):
	return(SEQ % (pos,len,seq))
