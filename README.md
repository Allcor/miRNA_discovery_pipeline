# miRNA_discovery_pipeline

NL:
Als je database of script* al op een andere plek (op genseq) staat en je wilt 'm ook hier gebruiken, gebruik dan liever een 'symbolic link' naar het bestand dan dat je het hele bestand kopieert. Dit kun je doen met 

user@pc:$ cp -l <bestand> <plek waar de link moet komen>

Op deze manier kun je (een beetje) ruimte besparen op de harde schijf!

*: of wat voor bestand dan ook

============================================

EN:
If your database or script* already exists (on genseq) and you want to use it here as well, please use a 'symbolic link' rather than a copy of the file. You can do this by issuing the command 

user@pc:$ cp -l <file> <directory where you want to have it>

This way, you can preserve (a little) disk space!

*: or whatever file you want to use



hint: if you want to run multiple task at the same time, (for all samples at the same time) run snakemake with the -j # command where # is a number that corresponds to the amount of tasks being run. 
