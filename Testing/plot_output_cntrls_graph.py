import matplotlib.pyplot as plt
import csv
plt.figure(dpi=150)
state_var=["netT", "Mx", "My", "Mz"]
for stt in state_var:
	x = []
	y = []

	with open(('../Outputs/'+stt+'.txt'),'r') as csvfile:
	    plots = csv.reader(csvfile, delimiter=',')
	    for row in plots:
	        x.append(float(row[0]))
	        y.append(float(row[1]))
	#print(x)
	#print(y)
	
	plt.plot(x,y,label = ("Pertubrations in "+stt))

	plt.scatter(x, y,label= "CGL Nodes", color= "red", marker= "*", s=10) 
	plt.xlabel('Time')
	plt.ylabel('Amplitude')
	plt.title('IPOPT Controls')
	#plt.legend()
plt.savefig('IPOPT output cntrls.png')
plt.show()