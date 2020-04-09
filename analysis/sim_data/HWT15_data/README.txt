1. Redshifts from MRPlancksnaplist.txt
    1. snap  a      z    t/t0    t/yr
    2. 11 0.09  9.72 0.04 4.91e+08
    3. 12 0.10  8.93 0.04 5.50e+08
    4. 13 0.11  8.22 0.04 6.16e+08
    5. 15 0.13  6.97 0.06 7.66e+08
    6. 17 0.14  5.92 0.07 9.47e+08
    7. 19 0.17  5.03 0.08 1.16e+09
2. Database in http://gavo.mpa-garching.mpg.de/MyMillennium/MyDB
    1. Example:
        1. select Type, stellarMass
        2.     from Henriques2015a..MRscPlanck1    
        3. where snapnum=11 and StellarMass > 0.01
    2. Example for 2020:
        1. select Type, stellarMass
        2.     from Henriques2020a..MR   
        3. where snapnum=19 and StellarMass > 0.01
3. Units in http://gavo.mpa-garching.mpg.de/MyMillennium/Help?page=databases/henriques2015a/mrscplanck1 
    1. Convert to [Msun] by multiplying with stellarMass with MassUnits = 1e10 / hubble
4. Read with H = np.loadtxt('snap_11.txt',skiprows=11,delimiter=',')