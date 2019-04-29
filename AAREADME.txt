FOLDERS:

"ACEDATA" -> ACE rawdata
"getace.m" generates event data files in the folder "DATA"  


Reconstruction Pipeline:

1. use "getace.m"

type "mygui" to start a graphical user interface

2. choose event

use buttons to...
3. Plot data (calls "dataplot.m")

4. dht and MVUB analysis ("prefr.m")  

5. GS axis  ("optcloud.m", calls several other functions)

6. GS solver ("hu12n.m", calls "hu34n.m")


figures and results are saved in the event directories
e.g. ACE_303_306_05

(means event is from ACE spacecraft, full data interval selected from 
day of year 303 to 306 in year 2005)


The files
"ALLCASE.par" and "CASE1.par" (in Event - Directory)
control the reduction process!



ALLCASE.PAR looks like this: 

"
ACE_303_306_05   %event directory
0			% 1 to write data file
3			% porder - order of fitting polynomial
21 131 16		% nx ny (-)dmid
1 0.1 0.15		% get_Ab (0:yes; -1: sum for < Ab; +1: sum for > Ab) (+)dAl0 (-)dAr0
"

These parameters can be controlled using "mygui".




