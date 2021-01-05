Repositorio con los archivos del TFM-Bioinformática y Bioestadística Area 1 aula 1 de título:
Mejora metodológica de modelos de transmisión de COVID19 mediante la inclusión de las características demográficas y automatización de la búsqueda de fechas clave para la propagación de la pandemia

En primer lugar, debemos llevar a cabo el curado de los datos. En este repositorio, en la carpeta "data/", encontramos datos para varios paises en la primera ola, y para algunos en la segunda (Alemania, Francia y España). Una vez tenemos los datos, podemos llevar a cabo el análisis deseado.

Si se pretende estudiar la repercusión de una medida en particular, se introducirá en el archivo "interventions_template.csv" para crear nuestro marco temporal de la medida. 
Si se pretende encontrar una fecha significativa, se especificará el intervalo de tiempo en el que buscar una fecha en el script "date_grid.R". Después, se introducirán estas fechas en el archivo "interventions_template.csv"
Para llevar a cabo la predicción, se empleará el script "base.R", especificando el nombre del pais que queremos analizar (debe coincidir con el nombre del archivo de los datos), así como el nombre de nuestro archivo de intervenciones.

Una vez entrenado el modelo, se obtendrán los resultados en forma de gráficas y workspace de R mediante el script "pred.R"
