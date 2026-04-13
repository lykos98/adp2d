# ADP for image segmentation
Using Advanced Density Peak algorithm for satellite image deblending.

Set variables in ```fits_segmentation.py``` to point to the correct paths to the image. 
Supports either a background rejection based on the [ASTERISM](https://academic.oup.com/mnras/article/463/3/2939/2646553) paper

## Compile 
```make``` :)

## Usage 
``` python3 fits_segmentation.py ```

## API Documentation
Detailed API documentation for the adp2d module is available in [API.md](API.md)

## Suggestions
1) Test convoluzione preliminare dell'immagine con un kernel gaussiano con sigma 0.5 (Simil Asterism)

2) Ripetere run di ADP e vedere se il problema del core delle stelle viene mitigato

3) Creare un tile test con sorgenti simulate a mano nel seguente modo:

- Solo stelle: prima due stelle lontane ed iniziare ad avvicinarle. 
Ripetere ADP ad ogni step e vedere nel caso di due solo sorgenti quando ADP si rompe nel deblendare (sempre se si rompe). 
Raggiunto un overlap inferiore alla dimensione tipica della PSF capire se non deblenda più (giusto) oppure continua a deblendare (in teoria sbagliato)

- Solo stelle: aumentare il numero di sorgenti progressivamente e ripetere i test del punto precedente. 
Ovviamente le sorgenti dovranno avere dimensioni random (in raggio multipli della PSF). Per le stelle non ha senso parlare di orientamento

- Solo galassie: ripetere quanto fatto per due stelle iniziando da dei profili di brillanza semplici (tipo Sersic)

- Solo galassie: aumentare il numero di galassie e ripetere (in questo caso con orientazioni completamente random (eventualmente anche con profile più 
complicati tipo Sersic + bulge)

- Ibrido: stella e galassia (semplice). Stessa procedura del caso due stelle o due galassie

- Ibrido: campo reale con stelle e galassie.

NOTE: L'immagine simulata non deve essere grande per fare i test e il numero massimo di sorgenti da trattare sarà al massimo 20/30.

