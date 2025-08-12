This evaporation_experiments_identification_README.txt file was generated on 2021-04-30 by SASCHA IDEN

GENERAL INFORMATION

1. Title of Dataset
Capillary, film, and vapor flow in bare soil evaporation: Identification using experimental data

2. Author Information 
Name: Sascha C. Iden (1), Efstathios Diamantopoulos (2), Wolfgang Durner (1)
(1) Technische Universität Braunschweig, Institute of Geoecology, Division of Soil Science and Soil Physics, Langer Kamp 19c, 38106 Braunschweig, Germany
(2) Department of Plant and Environmental Science, University of Copenhagen, Denmark

3. Contact Information
Sascha C. Iden
Email: s.iden@tu-braunschweig.de

4. Information about funding sources that supported the collection of the data: 
German Research Foundation (DFG), grant DU283/10-1

5. Date Created
06.02.2015 - 06.03.2015 (experiments), 01.07.2020 - 25.07.2020 (simulations)

6. Geographic location of data generation
Technische Universität Braunschweig, Institute for Geoecology, Langer Kamp 19c, 38106 Braunschweig, Germany

7. Related publication
This dataset belongs to the OA publication: 
Authors: Iden, S.C., E. Diamantopoulos, and W. Durner 
Title: Capillary, film, and vapor flow in transient bare soil evaporation (2): Experimental identification of hydraulic conductivity in the medium to dry moisture range
Journal: Water Resources Research
Year: 2021
doi: 10.1029/2020WR028514


ABSTRACT
Evaporation experiments are frequently used to identify the hydraulic properties of soils by inverse simulations with the Richards equation or simplified calculation methods. Evaporation experiments with an extended instrumentation were conducted in a temperature-controlled lab. Evaporation rates were measured gravimetrically, soil water pressure head was measured using mini-tensiometers and relative humidity sensors, and soil temperature was measured using thermocouples. The measurements were evaluated by inverse modeling with the Richards equation and the soil hydraulic properties, i.e. the water retention and hydraulic conductivity curves, were identified. The dataset contains all experimental data and the results of the inverse simulation, i.e. fitted time series of pressure head and the identified soil hydraulic properties of two soils (sand and silt loam).


SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: 
The data are published under the creative commons license CC-BY 4.0

2. Links to publications that cite or use the data: 
https://doi.org/10.1029/2020WR028514


METHODOLOGICAL Information
The evaporation experiments were conducted in a temperature-controlled lab at the Institute of Geoecology, TU Braunschweig. The soil columns were packed and instrumented with (i) minitensiometers, (ii) relative humidity sensors, (iii) temperature sensors. Evaporation rates were measured gravimetrially, the entire system was put on a scale and the mass was recorded continuously. The inverse simulations were performed with the Hydrus-1D software code which is available at https://www.pc-progress.com/en/Default.aspx?hydrus-1d


DATA & FILE OVERVIEW

1. Folder sand_measurements

FILE: sand_ambient.txt
DESCRIPTION: ambient conditions (temperature, humidity) in the laborotory as shown in figure 1b,d 

FILE: sand_evap_cum.txt
DESCRIPTION: measured cumulative evaporation (cm) as shown in figure 1a 

FILE: sand_evap_rate.txt
DESCRIPTION: measured evaporation rate (cm/d) as shown in figure 1a 

FILE: sand_head_tensiometer.txt
DESCRIPTION: measured pressure head in the soil (tensiometer data) as shown in figure 1c 

FILE: sand_relative_humidity.txt
DESCRIPTION: measured relative humidity in the soil (rotronic hygroclip) as shown in figure 1d 

FILE: sand_temp.txt
DESCRIPTION: measured soil temperature as shown in figure 1b 

2. Folder sand_inverse

FILE: sand_head_fit_m1.txt
DESCRIPTION: fitted pressure head data using model M1 as shown in figure 3a 

FILE: sand_head_fit_m2.txt
DESCRIPTION: fitted pressure head data using model M2 as shown in figure 3b 

FILE: sand_head_fit_m3.txt
DESCRIPTION: fitted pressure head data using model M3 as shown in figure 3c 

FILE: sand_head_fit_m4.txt
DESCRIPTION: fitted pressure head data using model M4 as shown in figure 3d 

FILE: sand_head_measured.txt
DESCRIPTION: measured pressure head data (tensiometer and hygrometer) as shown in figure 3a-d 

FILE: sand_shp_fit_m1.txt
DESCRIPTION: identified soil hydraulic properties using model M1 as shown in figure 3e,f 

FILE: sand_shp_fit_m2.txt
DESCRIPTION: identified soil hydraulic properties using model M2 as shown in figure 3e,f 

FILE: sand_shp_fit_m3.txt
DESCRIPTION: identified soil hydraulic properties using model M3 as shown in figure 3e,f 

FILE: sand_shp_fit_m4.txt
DESCRIPTION: identified soil hydraulic properties using model M4 as shown in figure 3e,f 


3. Folder siltloam_measurements

FILE: siltloam_ambient.txt
DESCRIPTION: ambient conditions (temperature, humidity) in the laborotory as shown in figure 2b,d 

FILE: siltloam_evap_cum.txt
DESCRIPTION: measured cumulative evaporation (cm) as shown in figure 2a 

FILE: siltloam_evap_rate.txt
DESCRIPTION: measured evaporation rate (cm/d) as shown in figure 2a 

FILE: siltloam_head_tensiometer.txt
DESCRIPTION: measured pressure head in the soil (tensiometer data) as shown in figure 2c 

FILE: siltloam_relative_humidity.txt
DESCRIPTION: measured relative humidity in the soil (rotronic hygroclip) as shown in figure 2d 

FILE: siltloam_temp.txt
DESCRIPTION: measured soil temperature as shown in figure 2b 

4. Folder siltloam_inverse

FILE: siltloam_head_fit_m1.txt
DESCRIPTION: fitted pressure head data using model M1 as shown in figure 5a 

FILE: siltloam_head_fit_m2.txt
DESCRIPTION: fitted pressure head data using model M2 as shown in figure 5b 

FILE: siltloam_head_fit_m3.txt
DESCRIPTION: fitted pressure head data using model M3 as shown in figure 5c 

FILE: siltloam_head_fit_m4.txt
DESCRIPTION: fitted pressure head data using model M4 as shown in figure 5d 

FILE: siltloam_head_measured.txt
DESCRIPTION: measured pressure head data (tensiometer and hygrometer) as shown in figure 5a-d 

FILE: siltloam_shp_fit_m1.txt
DESCRIPTION: identified soil hydraulic properties using model M1 as shown in figure 5e,f 

FILE: siltloam_shp_fit_m2.txt
DESCRIPTION: identified soil hydraulic properties using model M2 as shown in figure 5e,f 

FILE: siltloam_shp_fit_m3.txt
DESCRIPTION: identified soil hydraulic properties using model M3 as shown in figure 5e,f 

FILE: siltloam_shp_fit_m4.txt
DESCRIPTION: identified soil hydraulic properties using model M4 as shown in figure 5e,f 

Glossary / Abbreviations 
The following terms occur in the data files:

ecum      cumulative evaporation from soil (cm)
erate     evaporation rate (cm/d)
headobs   measured soil water pressure head (cm)
head_[x]  fitted soil water pressure head (cm)
log10(K)  common log of total isothermal soil hydraulic conductivity (cm/d)
relH[x]   relative humidity in soil (%)
relHair   relative humidity of lab air (%) 
T[x]      soil temperature (°C)
Tair      temperature of lab air (°C)
theta     volumetric soil water content (cm/cm) 
time      time (d) 



