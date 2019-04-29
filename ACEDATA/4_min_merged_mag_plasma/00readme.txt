ftp://nssdcftp.gsfc.nasa.gov/spacecraft_data/ace/4_min_merged_mag_plasma/00readme.txt


This data set was created at NSSDC after the downloading from the ACE 
Science Center of 4-min IMF data (from MAG instrument, PI: N.F. Ness) and 
64-s solar wind plasma parameters (from SWEPAM instrument, PI: D. McComas).
The plasma data were then averaged to 4-min resolution, taking any 64-s 
value whose start minute fell between the start minute of the IMF data and 
4 minutes later.  The data accessible herein were subsequently time shifted 
and used for making hourly averages "at Earth time" for adding to the OMNI
data set.  Time-shifted 4-min data are not accessible.  The Earth-time ACE 
hourly IMF data are accessible with concurrent hourly Earth-time Wind data
and IMP 8 data at http://nssdc.gsfc.nasa.gov/ftphelper/imp_wind_ace.html. 

The ACE MAG and SWEPAM Level 2 used herein are fully described at the ACE 
Science Center (http://www.srl.caltech.edu/ACE/ASC/level2/index.html) and 
at PI sites available therefrom.  Much ACE data and services are accessible 
from the ASC pages.

The contents of each of the 4-min merged records -

    WORD  FORMAT    MEANING                              UNITS/COMMENTS
                 
     1     I4       Year                                1998,1999, etc. 
     2     I4       Decimal Day                         January 1 = Day 1 
     3     I3       Hour                                0,1,...,23
     4     I3       Minute                              0,4,...,56

     5     I9       X GSE                                 km                        
     6     I9       Y GSE                                 km                        
     7     I9       Z GSE                                 km                        

     8    F9.3      Field Magnitude Average |B|         1/N SUM |B|, nT
     9    F9.3      BX GSE-Coordinate System             nanoteslas
    10    F9.3      BY GSE-Coordinate System             nanoteslas
    11    F9.3      BZ GSE-Coordinate System             nanoteslas
    12     I4       Number of vectors                    number   
    13     I3       Quality Flag                         See footnote

    14    F7.5      Plasma Proton density                [n/cc]
    15    F9.0      Plasma Temperature                   degrees, K   
    16    F7.4      He/H ratio
    17    F9.2      Plasma flow speed                    km/s
    18    F9.2      Vx Velocity                          km/s
    19    F9.2      Vy Velocity                          km/s
    20    F9.2      Vz Velocity                          km/s


Footnote:
Quality = 0 Normal data
        = 1 Spacecraft Maneuver & subsequent high-nutation period (~4 hr)
        = 2 Bad data/missing data

------------------------------------------------------------------
Related data and directories:
SPDF Data and Orbits Services <http://spdf.gsfc.nasa.gov/>
NSSDC Master Catalog <http://nssdc.gsfc.nasa.gov/nmc/sc-query.html>

NSSDC/SPDF Contact: Natalia Papitashvili <Natalia.E.Papitashvili.1@nasa.gov>,
 <request@mail630.gsfc.nasa.gov>

Please acknowledge NASA's National Space Science Data Center and
 Space Physics Data Facility, for data usage.

Authorizing NASA Official: Dr. R.E. McGuire, Head, SPDF, NASA Goddard Space Flight Center
                        301-286-7794, e-mail:  Robert.E.McGuire@nasa.gov
 
------------------------------------------
