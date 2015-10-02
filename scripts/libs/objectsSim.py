class pos:
	def __init__(self,n,s,e):
		self.s=s
		self.e=e
		self.id=n

class mirna:
	def __init__(self):
		self.p5=0
		self.p3=0
	def addp5(self,n,s,e):
		#print "adding p5"
		self.p5=pos(n,s,e)
		#print self.p5
	def addp3(self,n,s,e):
		self.p3=pos(n,s,e)
