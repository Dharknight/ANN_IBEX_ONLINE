"""import os
import sys
import subprocess
import pandas as pd
from keras.models import load_model


route = os.path.dirname(os.path.abspath(__file__))
route += "/"

FNULL = open(os.devnull, 'w') 

filename = sys.argv[1]

# Imprimir el contenido del archivo recibido
#with open(filename, 'r') as file:
#    print("Contenido del archivo recibido:\n")
#    print(file.read())
    
if (len(sys.argv) >= 3):
  seeds = sys.argv[2]
else:
  seeds = ""
if (len(sys.argv) >= 4):
  timeout = sys.argv[3]
else: 
  timeout = ""

#Lectura y formato del problema y sus datos
######################
######################



query = str(route) + "foo.exe " + str(filename)
output = subprocess.check_output(query,shell=True,)

header = "Var,Res,JG1,JG2,JG3,JG4,JG5,JG6,JG7,JG8,JG9,JG10,JG11,JG12,JG13,JG14,JG15,JG16,JG17,JG18,JG19,JG20,JG21,JG22,JG23,JG24,JG25,JG26,JG27,JG28,JG29,JG30,JG31,JG32,JG33,JG34,JG35,JG36,JG37,JG38\n"

file = open("/home/abel/TESIS_ANN_ONLINE/env_py2/ibex-lib/temp.txt","r")
Lines = file.readlines()
file.close()

data = Lines[0]
data = data.replace(" ; ", ",")
data = data.strip("((")
data = data.strip("))")

file = open("/home/abel/TESIS_ANN_ONLINE/env_py3/temp.csv","w")
file.write(header)
file.write(data)
file.close()


#Red Neuronal
######################
######################
#print("NN")

model = load_model('/home/abel/TESIS_ANN_ONLINE/env_py3/MFV3.h5')

Test = pd.read_csv("/home/abel/TESIS_ANN_ONLINE/env_py3/temp.csv", index_col=False)
y_predtest = model.predict(Test)


#Formateo del output
######################
######################
#print("FO")

cont = 0
for i in y_predtest:
  if(i[0]>i[1] and i[0]>i[2]):
    y_predtest[cont] = [1, 0, 0]
  if(i[1]>i[0] and i[1]>i[2]):
    y_predtest[cont] = [0, 1, 0]
  if(i[2]>i[1] and i[2]>i[0]):
    y_predtest[cont] = [0, 0, 1]
  cont+=1

#print(y_predtest)

#Solve
######################
######################
  
query = str(route) + "/"

#LSmear
if (y_predtest[0][0] == 1):
  #Resolucion
  solver = "bis:LSmear"
  query = query + "ibexoptLSmear " + filename + " " + seeds + " " + timeout
  #output = subprocess.check_output(query,shell=True,)
  #print("LSmear")

#RRobin
elif (y_predtest[0][1] == 1):
  #Resolucion
  solver = "bis:RoundRobin"
  query = query + "ibexoptRRobin " + filename + " " + seeds + " " + timeout
  #output = subprocess.check_output(query,shell=True,)
  #print("RoundRobin")

#LFO
elif (y_predtest[0][2] == 1):
  #Resolucion
  solver = "bis:LargestFirst"
  query = query + "ibexoptLFO " + filename + " " + seeds + " " + timeout
  #output = subprocess.check_output(query,shell=True,)
  #print("LargestFirst")

response = f"{solver}\n"
print(response, end='')"""

#Reporte
######################
######################
  
#limpieza
#output = str(output).split("\\n")
#output = output[-4]
  
#output = str(filename) + " " + str(solver) + " " + str(output)
#output += "\n"

#output = output.replace(" c:", " boxes:")
#output = output.replace(" p:", " sol:")
#output = output.replace(" v:", " val:")

#file = open("results.txt","a")
#file.write(output)
#file.close()

#file = open("rtemp.txt","w")
#file.write(output)
#file.close()
#print(output)
#escritura

import socket
import os
import subprocess
import sys
import pandas as pd
from keras.models import load_model

HOST = '127.0.0.1'
PORT = 8080
BUFFER_SIZE = 1024

route = os.path.dirname(os.path.abspath(__file__)) + "/"
header = "Var,Res,JG1,JG2,JG3,JG4,JG5,JG6,JG7,JG8,JG9,JG10,JG11,JG12,JG13,JG14,JG15,JG16,JG17,JG18,JG19,JG20,JG21,JG22,JG23,JG24,JG25,JG26,JG27,JG28,JG29,JG30,JG31,JG32,JG33,JG34,JG35,JG36,JG37,JG38\n"

model = load_model('/home/abel/TESIS_ANN_ONLINE/env_py3/MFV3.h5')

with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
    s.bind((HOST, PORT))
    s.listen()
    print(f"Servidor escuchando en {HOST}:{PORT}")

    conn, addr = s.accept()
    with conn:
        print(f"Conectado por {addr}")
        while True:
            data = conn.recv(BUFFER_SIZE)
            if not data:
                break

            filename = data.decode().strip()
            if not os.path.isfile(filename):
                error_msg = f"Error: File {filename} does not exist."
                print(error_msg, file=sys.stderr)
                conn.sendall(error_msg.encode())
                continue

            try:
                query = f"{route}foo.exe {filename}"
                output = subprocess.check_output(query, shell=True)

                with open("temp.txt", "r") as file:
                    lines = file.readlines()

                if not lines:
                    error_msg = "Error: No data found in the temp.txt file."
                    print(error_msg, file=sys.stderr)
                    conn.sendall(error_msg.encode())
                    continue

                data = lines[0].replace(" ; ", ",").strip("((").strip("))")

                with open("temp.csv", "w") as file:
                    file.write(header)
                    file.write(data)

                test_data = pd.read_csv("temp.csv", index_col=False)
                y_predtest = model.predict(test_data)

                cont = 0
                for i in y_predtest:
                    if i[0] > i[1] and i[0] > i[2]:
                        y_predtest[cont] = [1, 0, 0]
                    elif i[1] > i[0] and i[1] > i[2]:
                        y_predtest[cont] = [0, 1, 0]
                    elif i[2] > i[1] and i[2] > i[0]:
                        y_predtest[cont] = [0, 0, 1]
                    cont += 1

                if y_predtest[0][0] == 1:
                    solver = "bis:LSmear"
                elif y_predtest[0][1] == 1:
                    solver = "bis:RoundRobin"
                elif y_predtest[0][2] == 1:
                    solver = "bis:LargestFirst"

                response = f"{solver}\n"
                print(f"Sending response: {response}", file=sys.stderr)
                conn.sendall(response.encode())

            except Exception as e:
                error_msg = f"Error: {e}"
                print(error_msg, file=sys.stderr)
                conn.sendall(error_msg.encode())
        



