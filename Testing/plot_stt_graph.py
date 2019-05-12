import matplotlib.pyplot as plt
import csv

state_var=["p", "q" , "r", "phi", "theta", "psi", "z", "Vz", "y", "Vy", "x", "Vx"]
#stt = "x"
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
	
	# plt.plot(x,y,label = ("Pertubrations in "+stt))
	# plt.scatter(x, y,label= "CGL Nodes", color= "red", marker= "*", s=10) 
	# plt.xlabel('Time')
	# plt.ylabel('Amplitude')
	# plt.title('Quad Rotor Dyamics')
	# plt.legend()
	# # plt.savefig('Pertubrations.png')
	# # plt.show()
	
	
	a = []
	b = []
	
	with open('../Outputs/'+stt+'_inp.txt','r') as csvfile:
	   plots = csv.reader(csvfile, delimiter=',')
	   for row in plots:
	       a.append(float(row[0]))
	       b.append(float(row[1]))
	#print(x)
	#print(y)
	plt.figure(dpi=150)
	plt.plot(x,y,label = "IPOPT")
	# plt.hold(True);
	plt.plot(a,b,label = "LQR")
	plt.scatter(a, b,label= "CGL Nodes", color= "green", marker= "*", s=10)
	plt.scatter(x, y,label= "CGL Nodes", color= "red", marker= "*", s=10)  
	plt.xlabel('Time')
	plt.ylabel('Amplitude')
	plt.title(stt)
	plt.legend()
	
	plt.savefig('Pertubrations in '+stt+'.png')
	plt.show()