# ANN_IBEX_ONLINE

![License](https://img.shields.io/github/license/tu_usuario/ANN_IBEX_ONLINE)
![Python Version](https://img.shields.io/badge/python-3.x-blue.svg)

## Descripción

**ANN_IBEX_ONLINE** es un proyecto orientado a la creación de una metodología online para la selección de bisectores para la resolución de problemas NCOP (Non-Convex Optimization Problems). Utiliza la librería **Ibex**, que es un potente conjunto de herramientas para la optimización y análisis de intervalos.

Este repositorio incluye todos los archivos que se usaron para modificar la libreria ibex y crear las diferentes propuestas para mi proyecto de defensa de titulo. Se explicara como se deben usar cada uno de los archivos modificados a partir de la libreria ibex original.

## Tabla de Contenidos

- [Instalación](#instalación)
- [Estructura del Repositorio](#estructura-del-repositorio)
- [Uso](#uso)
- [Contribuciones](#contribuciones)
- [Licencia](#licencia)
- [Agradecimientos](#agradecimientos)

## Instalación

Antes de comenzar, se decidio crear un entorno virtual con python para la instalacion de IBEX y la ANN. 

1. **Creacion de entornos virtuales**:
   1.1. Entorno virtual con python3:
         ```bash
            python3 -m venv env_py3
      
   1.2. Entorno virtual con python2:
      ```bash
         pip install virtualenv
      ```bash
         virtualenv -p python2.7 env_py27

3. **Instalacion IBEX entorno virtual con python2**
   Seguir indicaciones de documentacion IBEX (http://ibex-team.github.io/ibex-lib/).

4. **Clonar carpeta env_ann en entorno virtual con python3**

## Estructura del Repositorio
1. **Carpeta env_ann**:
   1.1 **Explicacion archivos**:
      - main.py: Script que contiene el modelo de la red neuronal artificial y crea el lado servidor con socket para la comunicacion entre IBEX y ANN.
      - MFV3.h5: Modelo de la red neuronal artificial.
      - temp.csv: archivo csv que contiene las caracteristicas del problema a resolver. Cantidad de variables, restricciones, jacobianos calculados, etc.
      - temp.txt: mismo contenido que temp.csv pero en formato .txt .
      - foo.exe: ejecutable para limpiar archivo que contiene el problema y calcular sus caracteristicas.
      - updated_problem.bch: archivo que se crea cuando se crea una instancia nueva del problema original.

3. **Carpeta env_lib_ibex**
   2.1 **Explicacion archivos**:
      -  carpeta benchs: igual a la original de la libreria IBEX, dentro de optim esta la carpeta easy, medium, hard, cada una de estas contiene los problemas de la libreria.
      -  carpeta optim: contiene los archivos modificados para cada propuestas. El original se encuentra en /src/optim.

## Uso
   Ya instalada la libreria ibex en su entorno virtual y la ANN igual, entonces se procede a explicar como ejecutar problemas usando ambas partes.

   1. **Ejecutar script python main.py que tambien trabaja como servidor**
      ```bash
         python3 main.py
   2. **Ejecutar problema con libreria IBEX**
      Antes de la ejecucion, la libreria ibex crea un ejecutable de la siguiente manera.
      ```bash
         waf ./install
      generando el archivo ibexopt ubicado en ./__build__/src/ibexopt
      
      Y para ejecutar un caso de prueba se usa el siguiente comando:
      ```bash
         ./__build__/src/ibexopt benchs/optim/medium/ex6_2_9.bch --random-seed=1

