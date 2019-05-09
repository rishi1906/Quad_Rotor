import matplotlib.pyplot as plt
import csv
plt.figure(dpi=150)
state_var=["netT", "Mx", "My", "Mz"]
for stt in state_var:
	x = []
	y = []

	with open(('../Outputs/'+stt+'_inp.txt'),'r') as csvfile:
	    plots = csv.reader(csvfile, delimiter=',')
	    for row in plots:
	        x.append(float(row[0]))
	        y.append(float(row[1]))
	#print(x)
	#print(y)
	
	plt.plot(x,y)

	plt.scatter(x, y, color= "red", marker= ".", s=10) 
	plt.xlabel('Time')
	plt.ylabel('Amplitude')
	plt.title('LQR input controls')
	# plt.legend()
plt.savefig('LQR_input_controls.png')
plt.show()


#a = []
#b = []

#with open('test_output_2.txt','r') as csvfile:
#    plots = csv.reader(csvfile, delimiter=',')
#    for row in plots:
#        a.append(float(row[0]))
#        b.append(float(row[1]))
#print(x)
#print(y)
#plt.figure(dpi=150)
#plt.plot(x,y,label = "Function Curve")
## plt.hold(True);
#plt.plot(a,b,label = "Derivative Curve")
#plt.scatter(a, b,label= "CGL Nodes", color= "green", marker= "*", s=10)
#plt.scatter(x, y,label= "CGL Nodes", color= "red", marker= "*", s=10)  
#plt.xlabel('Time')
#plt.ylabel('Amplitude')
#plt.title('Integral Weights Validation')
#plt.legend()

#plt.savefig('Integral_Weights_Validation.png')
#plt.show()