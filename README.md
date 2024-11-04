# Visor-Protein-Metamorphisms
## 1.Introducción y descripción
Este proyecto tiene el objetivo de crear un visor que integre los resultados de la base de datos del sistema https://github.com/CBBIO/protein-metamorphisms-is con el visualizador de estructuras PyMOL[1]. Este visor facilita el análisis manual de los alineamientos con el objetivo de seleccionar los posibles metamorfismos.

## 2.Instalación
Una vez clonado el repositorio, para instalar las dependencias necesarias, se puede usar el archivo environment.yml para crear un entorno virtual en Conda [2] usando “conda env create -f environment.yml”; o bien directamente usar el fichero requirements.txt

## 3.Consideraciones importantes
La aplicación usa los datos a partir de un fichero csv, cuyo formato es el obtenido a partir del sistema https://github.com/CBBIO/protein-metamorphisms-is y se puede ver un ejemplo del mismo en el fichero example.csv . Para ejecutar el visor con tus datos es necesario escribir la ruta de tu csv en la variable self.path del script visor.py

## 4.Guía de uso
Al ejecutar el script visor.py aparecerá una ventana con diferentes botones, para empezar el análisis e iniciar pymol es necesario pulsar el botón “Iniciar”. Una vez iniciada la ventana de pymol y cargados las primeras estructuras, se puede mover entre los diferentes subclústers hacia el siguiente alineamiento (con el botón “siguiente alineamiento” o la flecha “Right” del teclado) o el alinemiento anterior (con el botón “anterior alineamiento” o la flecha “Left” del teclado); o puede moverse entre los diferentes clúster, usando el botón “Siguiente clúster” (o la flecha “Down” del teclado) para el siguiente clúster del archivo csv, o usando el botón “Anterior clúster” (o la flecha “Up” del teclado) para el clúster anterior.

Para facilitar la anotación de los polimorfismos en los alineamientos se encuentran dos checkbox a la derecha de la aplicación, donde puede seleccionar si hay o no polimorfismo para modificar dicha anotación en el fichero csv, una vez terminado el análisis puede guardar los cambios en el fichero usando el botón “Guardar fichero” para sobreescribir su fichero csv

## 5.Referencias
[1] Schrödinger, L., & DeLano, W. (2020). PyMOL. Retrieved from http://www.pymol.org/pymol

[2] Anaconda Software Distribution. (2020). Anaconda Documentation. Anaconda Inc. Retrieved from https://docs.anaconda.com/
