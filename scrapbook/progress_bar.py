# Progress bar to visually show script progress
class ProgressBar:
	def __init__(self):
		self.percentComplete = 0
		self.prevOut = ""
	def addPercent(self, percent):
		if self.percentComplete + percent < 100:
			self.percentComplete += percent
	def draw(self):
		for i in range(len(self.prevOut)):
			sys.stdout.write("\b")
		self.prevOut = "["
		for i in range(int(self.percentComplete)/2):
			self.prevOut += "-"
		for i in range(50 - (int(self.percentComplete)/2)):
			self.prevOut += " "
		self.prevOut += "] " + str(int(self.percentComplete)) + "%"
		sys.stdout.write(self.prevOut)
		sys.stdout.flush()
	def close(self):
		self.percentComplete = 100
		for i in range(len(self.prevOut)):
			sys.stdout.write("\b")
		self.prevOut = "["
		for i in range(int(self.percentComplete)/2):
			self.prevOut += "-"
		for i in range(50 - (int(self.percentComplete)/2)):
			self.prevOut += " "
		self.prevOut += "] " + str(int(self.percentComplete)) + "%"
		sys.stdout.write(self.prevOut)
		sys.stdout.write('\n')
		sys.stdout.flush()