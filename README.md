# Proyecto AutoFemurLoc

El proyecto AutoFemurLoc proporciona herramientas para la detección y medición automáticas del fémur fetal en imágenes de ultrasonido 2D. Utilizando técnicas de procesamiento de imágenes, este proyecto tiene como objetivo agilizar el proceso de mediciones fetales cruciales para la estimación de la edad gestacional y el monitoreo del crecimiento.

## Antecedentes

Basado en el método presentado en "Fully automatic segmentation and measurement of the fetal femur," published in the Proceedings 
Volume 10975 of the 14th International Symposium on Medical Information Processing and Analysis. 
DOI: 10.1117/12.2511534., AutoFemurLoc identifica los puntos finales del fémur en imágenes de ultrasonido. Cuando se dispone de metadatos de resolución espacial, se utilizan para proporcionar mediciones en milímetros.

## Instalación

Para configurar tu entorno para el proyecto AutoFemurLoc:

1. Asegúrate de que Python 3.6+ esté instalado en tu sistema.
2. Clona este repositorio en tu máquina local.
3. Navega al directorio clonado y crea un entorno virtual:

   ```bash
   python3 -m venv env

4. Activa el entorno virtual::

   ```bash
   source env/bin/activate

4. Instala las dependencias necesarias con el archivo **requirements.txt** incluido en este repositorio:

   ```bash
   pip install -r requirements.txt

## Uso

 Para utilizar la funcion sigue estos pasos :

 1. Importa la clase **AutoFemurLoc** de **AutoObst**:

    ```python
    from AutoObst import AutoFemurLoc

 2. Crea una instancia de **AutoFemurLoc** :

    ```python
    myObst = AutoFemurLoc()

 3. Llama al método **FindFemur** con la imagen de US 2D:
    
    ```python
    end_points = myObst.findFemur(image_data)

 4. Si se requiere se puede visaulizar los puntos con la función **set_markers**

    ```python
    myObst.set_markers(image_data,end_point)

 5. Si se tiene los metadatos con la resolución espacial se puede utilizar la función **meassure_length**:

    ```python
    meassure_mm = myObst.meassure_length(row_spacing,end_point)


## Ejemplo de uso 

Dentro del repositorio se provee un ejemplo de uso en unnotebook, también hay una imagen para probar el algoritmo propuesto. 
