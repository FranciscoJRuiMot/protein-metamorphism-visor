# Visor-Protein-Metamorphisms
## 1.Introducción y descripción
Este proyecto tiene el objetivo de crear un visor que integre los resultados de la base de datos del sistema https://github.com/CBBIO/protein-metamorphisms-is con el visualizador de estructuras PyMOL[1]. Este visor facilita el análisis manual de los alineamientos con el objetivo de seleccionar los posibles metamorfismos.

## 2.Instalación
Una vez clonado el repositorio, para instalar las dependencias necesarias, se puede usar el archivo environment.yml para crear un entorno virtual en Conda [2] usando:

Instalar herramienta de grafos en el sistema.

sudo apt-get install graphviz

```bash
conda env create -f environment.yml
```

## 3.Consideraciones importantes
La aplicación extrae la información a partir de una base de datos, cuya estructura proviene del sistema https://github.com/CBBIO/protein-metamorphisms-is. Se puede ver un ejemplo de la consulta utilizada en **visor/db_consults**.

Para ejecutar el visor con tus propios datos puedes modificar las credenciales y la configuración de la base de datos en **visor/config/config.yaml**

## 4.Guía de uso
### 4.1 Inicio
El visor se ejecuta de la siguiente forma:

```bash
cd visor
python main.py
```

Al ejecutar el fichero aparecerá una ventana con diferentes botones, para empezar el análisis e iniciar PyMOL es necesario pulsar el botón **Iniciar**. 

### 4.2 Navegación entre los alineamientos

Una vez iniciada la ventana de PyMOL y cargados las primeras estructuras, puedes navegar por los diferentes clústers y alineamientos de la siguiente manera:

- Moverse entre alineamientos:

    - Siguiente alineamiento/subclúster: Botón **Siguiente alineamiento** o tecla ↓ (Down).

    - Alineamiento/subclúster anterior: Botón **Anterior alineamiento** o tecla ↑ (Up).

- Moverse entre clústeres:

    - Siguiente clúster: Botón **Siguiente clúster** o tecla → (Right).

    - Clúster anterior: Botón **Anterior clúster** o tecla ← (Left).

### 4.3 Alineamiento de las estructuras

Las diferentes proteínas se cargan sin alinear para permitir al usuario que elija el algoritmo de alineamiento que desee en cada caso. El visor te permite alinear las estructuras que aparecen por pantalla en cualquier momento usando los tres diferentes algoritmos que utiliza *Protein-metamorphisms-is*.

- Combinatorial Extension (CE) algorithm: Botón **CE-Align** o tecla **C**

- Universal Structure (US) alignment algorithm: Botón **US-Align** o tecla **U**

- Flexible structure AlignmenT by Chaining Aligned fragment pairs allowing Twists (FATCAT) algorithm: Botón **FATCAT-Align** o tecla **F**

### 4.4 Visualización de las métricas de interés

El visor muestra diversas tablas con las métricas obtenidas, y almacenadas en la base de datos, para el clúster que se esté visualizando en ese momento.

Dichas tablas están generadas utilizando *Treeviews* de tkinter, estos elementos son personalizables mediante el fichero de configuración **visor/config/config.yaml**. Donde se pueden crear subdataframes de la consulta, seleccionando las columnas de inteŕes y asignándolas a los *Treeviews*, que adicionalmente, pueden añadirse o eliminarse siguiente la misma estructura definida.

### 4.5 Anotación de metamorfismos y comentarios

Para facilitar la anotación de los metamorfismos en los alineamientos se encuentran dos checkbox en la zona izquierda de la aplicación, donde puedes seleccionar si hay o no polimorfismo para modificar dicha anotación en la base de datos.

Además, el visor te permite escribir comentarios sobre los alineamientos, para guardar dichos comentarios en el alineamiento actual es necesario usar el botón **Guardar comentarios actuales**

Siempre que desees guardar los cambios en la base de datos, tanto la anotación como los comentarios, puedes usar el botón **Actualizar BD** para actualizarla.

## 5.Referencias
[1] Schrödinger, L., & DeLano, W. (2020). PyMOL. Retrieved from http://www.pymol.org/pymol

[2] Anaconda Software Distribution. (2020). Anaconda Documentation. Anaconda Inc. Retrieved from https://docs.anaconda.com/
