```python
'The following script using monthly averages of temperature data to simulate activity and energetics of a green salamander (Aneides aeneus)'
#-------- LIBRARIES---------#

from numpy import *
from math import *
from random import *
from pandas import *
import glob
import re
import itertools


#---------CONSTANTS----------#

Rv = 461.5 #J*K^-1*kg^-1
L = 2.5*10**6 #J per kg
gravity = 9.8 #meters per second
header = 'ncols        240\nnrows        231\nxllcorner    235489.9231464\nyllcorner    3803709.315658\ncellsize     890.16924311368\nNODATA_value -9999\n '

#---------DATAFRAMES----------#
results = pandas.DataFrame(columns = ['depth','threshold','body_size','skin_resistance','acclimation_status','epoch','humidity','month','mean_activity','sd_activity','IUCN_activity','sd_IUCN_activity','captures_activity','sd_captures_activity','mean_energy','sd_energy','IUCN_energy','sd_IUCN_energy','captures_energy','sd_captures_energy','%_positive_energy','%_positive_IUCN','%_positive_captures','extinct_in_range','Teh_active','Teh_inactive','Teh_total','dewpoint_temp'])
annual_results = pandas.DataFrame(columns = ['depth','threshold','body_size','skin_resistance','acclimation_status','epoch','humidity','mean_activity','sd_activity','IUCN_activity','sd_IUCN_activity','captures_activity','sd_captures_activity','mean_energy','sd_energy','IUCN_energy','sd_IUCN_energy','captures_energy','sd_captures_energy','%_positive_energy','%_positive_IUCN','%_positive_captures','extinct_in_range','Teh_active','Teh_inactive','Teh_total'])

#---------FUNCTIONS---------#
def numericalSort(value):
    'a function for sorting'
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

def create_env_lists(specific_file):
    'a function for creating lists of temperatures from asc files'
    data_list = []
    list_of_files = glob.glob(specific_file)
    sorted_list_of_files = sorted(list_of_files, key=numericalSort)
    indexer = 0
    for i in range(12):
        temporary_hourly_list = []
        for j in range(24):
            temporary_list = []
            mod1 = open(sorted_list_of_files[indexer])
            mod2 = mod1.readlines()
            for i in range(6):
                mod2.pop(0)
            for line in mod2:
                line = line.strip()
                rows = line.split()
                temporary_list.append(rows)
            for i in range(len(temporary_list)):
                for j in range(len(temporary_list[i])):
                    temporary_list[i][j] = float(temporary_list[i][j])
                    temporary_list[i][j] = round(temporary_list[i][j],4)
            temporary_hourly_list.append(temporary_list)
            indexer += 1
        data_list.append(temporary_hourly_list)
    return data_list

def read_coordinates(coordinate_file):
    'a function for reading coordinates'
    list_of_coords = []
    mod1 = open(coordinate_file)
    mod2 = mod1.readlines()
    for k in range(6):
        mod2.pop(0)
    for line in mod2:
        line = line.strip()
        rows = line.split()
        list_of_coords.append(rows)
    for i in range(len(list_of_coords)):
        for j in range(len(list_of_coords[i])):
            list_of_coords[i][j] = float(list_of_coords[i][j])
    return list_of_coords

def return_stats_from_coords(list_of_values,list_of_coords):
    'a function for returns statistics from specific coordinates'
    list_of_values_at_coords = []
    stats = []
    for m in range(len(list_of_coords)):
        if list_of_coords[m] > 0:
            if list_of_values[m] > -9999.0:
                list_of_values_at_coords.append(list_of_values[m])
            else:
                pass
        else:
            pass
    stats.append(mean(list_of_values_at_coords))
    stats.append(std(list_of_values_at_coords))
    return stats

def return_values_from_coords(list_of_values,list_of_coords):
    'a function for returns values from specific coordinates'
    list_of_values_at_coords = []
    for m in range(len(list_of_coords)):
            if list_of_coords[m] > 0:
                if list_of_values[m] > -9999.0:
                    list_of_values_at_coords.append(list_of_values[m])
                else:
                    pass
            else:
                pass
    return list_of_values_at_coords

def return_positive_energy(energy_list):
    'a function for returns regions with positive energy balance from a list'
    total_length = 0.0
    for j in range(len(energy_list)):
        if energy_list[j] > -9999:
            total_length += 1.0
    count = 0.0
    if total_length == 0:
        return 0.0
    else:
        for i in range(len(energy_list)):
            if energy_list[i] > 0.0:
                count += 1.0
            else:
                pass
        return count/total_length

def convert_to_annual(monthly_list,annual_list,month):
    'a function that converts monthly to annual values'
    if month == 0 or month == 2 or month == 4 or month == 6 or month == 7 or month == 9 or month == 11:
        new_monthly_list = [x * 31.0 for x in monthly_list]
        new_annual_list = [a + b for a, b in zip(annual_list, new_monthly_list)]
        return new_annual_list
    elif month == 1:
        new_monthly_list = [x * 28.0 for x in monthly_list]
        new_annual_list = [a + b for a, b in zip(annual_list, new_monthly_list)]
        return new_annual_list
    else:
        new_monthly_list = [x * 30.0 for x in monthly_list]
        new_annual_list = [a + b for a, b in zip(annual_list, new_monthly_list)]
        return new_annual_list

def convert_to_annual2(monthly_list,annual_list,month):
    'a function that converts monthly to annual values but ignores null values'
    if month == 0 or month == 2 or month == 4 or month == 6 or month == 7 or month == 9 or month == 11:
        for i in range(len(monthly_list)):
            if monthly_list[i] == -9999:
                pass
            else:
                monthly_list[i] = monthly_list[i] * 31.0
        new_annual_list = []
        for j in range(len(monthly_list)):
            if monthly_list[j] == -9999:
                new_annual_list.append(-9999)
            else:
                new_annual_list.append(monthly_list[j]+annual_list[j])
        return new_annual_list
    elif month == 1:
        for i in range(len(monthly_list)):
            if monthly_list[i] == -9999:
                pass
            else:
                monthly_list[i] = monthly_list[i] * 28.0
        new_annual_list = []
        for j in range(len(monthly_list)):
            if monthly_list[j] == -9999:
                new_annual_list.append(-9999)
            else:
                new_annual_list.append(monthly_list[j]+annual_list[j])
        return new_annual_list
    else:
        for i in range(len(monthly_list)):
            if monthly_list[i] == -9999:
                pass
            else:
                monthly_list[i] = monthly_list[i] * 30
        new_annual_list = []
        for j in range(len(monthly_list)):
            if monthly_list[j] == -9999:
                new_annual_list.append(-9999)
            else:
                new_annual_list.append(monthly_list[j]+annual_list[j])
        return new_annual_list

def local_extinction(annual_list,lipid_reserve):
    'a function that determines whether a location has become locally extirpated'
    for i in range(len(annual_list)):
        if annual_list[i] <= lipid_reserve:
            annual_list[i] = -9999
        else:
            pass
    return annual_list

def activity_extinction(activity_list,energy_list):
    'a function that converts monthly activity to extinction estimates'
    for i in range(len(energy_list)):
        if energy_list[i] == -9999:
            activity_list[i] = -9999
        else:
            pass
    return activity_list

def calculate_annual_activity(activity_list):
    'a function that calculates annual activity'
    for i in range(len(activity_list)):
        if activity_list[i] == -9999:
            pass
        else:
            activity_list[i] = activity_list[i]/365.0
    return activity_list

def extinct_in_range(energy_list,range_shape):
    'a function that estimates the proportion of areas extinct within an area defined by a raster'
    total_length = 0.0
    count = 0.0
    for j in range(len(range_shape)):
        if range_shape[j] >= 1.0:
            total_length += 1.0
            if energy_list[j] == -9999:
                count +=  1.0
            else:
                pass
        else:
            pass
    if total_length == 0:
        return 1.0
    else:
        return count/total_length

def make_and_export_grid(data_list,file_name,header):
    'a function that converts lists into asc files'
    _list = [data_list[i:i+231] for i in range(0, len(data_list), 231)]
    lines = []
    for row in _list:
            lines.append(' '.join(map(str, row)))
    result = '\n '.join(lines)
    f = open(file_name, 'w')
    f.writelines(str(header))
    f.writelines(str(result))
    f.close()

def dewpoint(temperature_list):
    'a function that calculates dewpoints from minimum temperatures'
    null_value  = -9999
    total_temp = []
    for i in range(len(temperature_list[0])):
        for j in range(len(temperature_list[0][0])):
            if temperature_list[6][i][j] > null_value:
                total_temp.append(temperature_list[6][i][j])
            else:
                pass
    return(min(total_temp))

def ignore_value_for_mean(list_with_values,value_to_ignore):
    'a function that ignores specified value when calculated the mean'
    list_without_values = []
    if mean(list_with_values) == value_to_ignore:
            list_without_values.append(value_to_ignore)
    else:
            for i in range(len(list_with_values)):
                if list_with_values[i] == value_to_ignore:
                    pass
                else:
                    list_without_values.append(list_with_values[i])
    return list_without_values

def convert_to_daily_average(annual_list):
    'a function that converts to daily average while ignoring null values'
    for i in range(len(annual_list)):
        if annual_list[i] == -9999:
            pass
        else:
            annual_list[i] = annual_list[i]/365
    return annual_list

def convert_to_integer(total_list):
    'a function that converts to integers'
    final_list = []
    for i in range(len(total_list)):
        if total_list[i] == -9999:
            final_list.append(-9999)
        else:
            total_list[i] = round(total_list[i],2)
            total_list[i] = int(total_list[i] * 100)
            final_list.append(total_list[i])
    return final_list

def run_sdm(temperatures_list,epoch_index,elevations,body_size_list,body_size_index,ri_list,ri_index,threshold_list,threshold_index,humidity_list,humidity_index,depth_list,depth_index):
    'a global function that runs the class for each month of the year and combines them into annual estimates'
    _annual_activity = [0] * len(temperatures_list[0][0][0]) * len(temperatures_list[0][0][0][0])
    _annual_energy = [0] * len(temperatures_list[0][0][0]) * len(temperatures_list[0][0][0][0])
    annual_Teh_active = []
    annual_Teh_inactive = []
    annual_Teh_total = []
    for month in range(len(temperatures_list[0])):
        dewpoint_temp = dewpoint(temperatures_list[0][month])
        sal = Individual(body_size_list[body_size_index],ri_list[ri_index],threshold_list[threshold_index])
        sal.calculate_nightly_activity(temperatures_list[epoch][month],temperatures_list[0][month],dewpoint_temp,elevations,humidity_list[humidity_index],depth_list[depth_index])
        if not sal.Teh_active:
            Teh_active = [0.0]
        else:
            Teh_active = [mean(sal.Teh_active)]
        Teh_inactive = [mean(sal.Teh_inactive)]
        Teh_total = [mean(sal.Teh_total)]
        annual_Teh_active.extend(Teh_active)
        annual_Teh_inactive.extend(Teh_inactive)
        annual_Teh_total.extend(Teh_total)
        IUCN_activity = return_stats_from_coords(sal.activity_grid,IUCN)
        captures_activity = return_stats_from_coords(sal.activity_grid,captures)
        IUCN_energy = return_stats_from_coords(sal.energy_grid,IUCN)
        captures_energy = return_stats_from_coords(sal.energy_grid,captures)
        temporary_dataframe= pandas.DataFrame([[depth_list[depth_index],threshold_list[threshold_index],body_size_list[body_size_index],ri_list[ri_index][0],list_of_acclimation_status[r_s],list_of_epochs[epoch],humidity_list[humidity_index],list_of_months[month],mean(ignore_value_for_mean(sal.activity_grid,-9999)),std(ignore_value_for_mean(sal.activity_grid,-9999)),IUCN_activity[0],IUCN_activity[1],captures_activity[0],captures_activity[1],mean(ignore_value_for_mean(sal.energy_grid,-9999)),std(ignore_value_for_mean(sal.energy_grid,-9999)),IUCN_energy[0],IUCN_energy[1],captures_energy[0],captures_energy[1],return_positive_energy(sal.energy_grid),return_positive_energy(return_values_from_coords(sal.energy_grid,IUCN)),return_positive_energy(return_values_from_coords(sal.energy_grid,captures)),extinct_in_range(_annual_energy,IUCN),Teh_active[0],Teh_inactive[0],Teh_total[0],dewpoint_temp]],columns = ['depth','threshold','body_size','skin_resistance','acclimation_status','epoch','humidity','month','mean_activity','sd_activity','IUCN_activity','sd_IUCN_activity','captures_activity','sd_captures_activity','mean_energy','sd_energy','IUCN_energy','sd_IUCN_energy','captures_energy','sd_captures_energy','%_positive_energy','%_positive_IUCN','%_positive_captures','extinct_in_range','Teh_active','Teh_inactive','Teh_total','dewpoint_temp'])
        global results
        results = results.append(temporary_dataframe)
        _annual_activity = convert_to_annual2(sal.activity_grid,_annual_activity,month)
        _annual_energy = convert_to_annual2(sal.energy_grid,_annual_energy,month)
    annual_activity = convert_to_daily_average(_annual_activity)
    annual_energy = convert_to_daily_average(_annual_energy)
    annual_energy_file = 'green_salamander_project/energy/annual_energy'+'_depth'+str(depth_list[depth_index])+'_threshold'+str(threshold_list[threshold_index])+'_mass'+str(body_size_list[body_size_index])+'_rs'+str(ri_list[ri_index][0])+'_acc'+str(list_of_acclimation_status[r_s])+'_epoch'+str(list_of_epochs[epoch])+'_humidity'+str(humidity_list[humidity_index])+'error.asc'
    annual_activity_file = 'green_salamander_project/activity/annual_activity'+'_depth'+str(depth_list[depth_index])+'_threshold'+str(threshold_list[threshold_index])+'_mass'+str(body_size_list[body_size_index])+'_rs'+str(ri_list[ri_index][0])+'_acc'+str(list_of_acclimation_status[r_s])+'_epoch'+str(list_of_epochs[epoch])+'_humidity'+str(humidity_list[humidity_index])+'error.asc'
    annual_energy_int = convert_to_integer(annual_energy)
    annual_activity_int = convert_to_integer(annual_activity)
    make_and_export_grid(annual_energy_int,annual_energy_file,header)
    make_and_export_grid(annual_activity_int,annual_activity_file,header)
    IUCN_annual_activity = return_stats_from_coords(annual_activity,IUCN)
    captures_annual_activity = return_stats_from_coords(annual_activity,captures)
    IUCN_annual_energy = return_stats_from_coords(annual_energy,IUCN)
    captures_annual_energy = return_stats_from_coords(annual_energy,captures)
    annual_Teh_active_no_zeros = ignore_value_for_mean(annual_Teh_active,0.0)
    temporary_dataframe_annual = pandas.DataFrame([[depth_list[depth_index],threshold_list[threshold_index],body_size_list[body_size_index],ri_list[ri_index][0],list_of_acclimation_status[r_s],list_of_epochs[epoch],humidity_list[humidity_index],mean(ignore_value_for_mean(annual_activity,-9999)),std(ignore_value_for_mean(annual_activity,-9999)),IUCN_annual_activity[0],IUCN_annual_activity[1],captures_annual_activity[0],captures_annual_activity[1],mean(ignore_value_for_mean(annual_energy,-9999)),std(ignore_value_for_mean(annual_energy,-9999)),IUCN_annual_energy[0],IUCN_annual_energy[1],captures_annual_energy[0],captures_annual_energy[1],return_positive_energy(annual_energy),return_positive_energy(return_values_from_coords(annual_energy,IUCN)),return_positive_energy(return_values_from_coords(annual_energy,captures)),extinct_in_range(annual_energy,IUCN),mean(annual_Teh_active_no_zeros),mean(annual_Teh_inactive),mean(annual_Teh_total)]],columns = ['depth','threshold','body_size','skin_resistance','acclimation_status','epoch','humidity','mean_activity','sd_activity','IUCN_activity','sd_IUCN_activity','captures_activity','sd_captures_activity','mean_energy','sd_energy','IUCN_energy','sd_IUCN_energy','captures_energy','sd_captures_energy','%_positive_energy','%_positive_IUCN','%_positive_captures','extinct_in_range','Teh_active','Teh_inactive','Teh_total'])
    global annual_results
    annual_results = annual_results.append(temporary_dataframe_annual)
    results.to_csv('green_salamander_project/results/climate_data_results_error.csv')
    annual_results.to_csv('green_salamander_project/results/annual_climate_data_results_error.csv')



#---------CLASSES-----------#

class Individual():
    def __init__(self,MASS,r_i,THRESHOLD):
        self.mass = MASS #grams
        self.diameter = 0.0016*log(self.mass) + 0.0061 #empirical formula for diameter from Riddell et al. 2017
        self.activity_threshold = self.mass - (THRESHOLD*self.mass) #grams
        self.mass_to_lose = self.mass - self.activity_threshold #grams
        self.Rs = r_i[0] #sec cm^-1
        self.rb = 0.0 #sec cm^-1
        self.activity = 0.0 #hours
        self.surface_area = 8.42*(self.mass**0.694) #cm^2
        self.EWL = 0.0 #g sec^-1
        self.T_eh = 0.0
        self.activity_grid = []
        self.activity_status = 0.0 #zero is active, one is inactive
        self.elevm = 760.0
        self.energy_status = 0.0
        self.energy_grid = []
        self.acclimation_status = r_i[1]
        self.Teh_active = []
        self.Teh_inactive = []
        self.Teh_total = []
        self.E_G = 1.0 #emmisitivity for soil, blackbody conditions
        self.E_S = 0.96 #emissivity of salamander
        self.A_L = 0.965 #absorptance of organism to longwave radiation

    def longwave_sky(self,temperature):
        'this function calculates longwave radiation from the skin based on Campbell and Norman 2010'
        return 1.0*(5.670373*10**-8 * (temperature + 273.15)**4) #88

    def longwave_ground(self,temperature):
        'this function calculates longwave radiation from the ground based on Campbell and Norman 2010'
        b = 5.670373*10**-8
        return self.E_G*b*(temperature+273.15)**4.

    def Rabs(self,temperature):
        'this function calculates the radiation absorbed from the sky and ground based on Campbell and Norman 2010'
        return (0.5*(self.A_L*(self.longwave_sky(temperature)+self.longwave_ground(temperature))))

    def radiative_conductance(self,temperature):
        'this function calculates radiative conductance from Campbell and Norman 2010'
        return (4.*(5.670373*10**-8)*(temperature+273.15)**3.)/29.3

    def calculate_Teh(self,r_i,r_b,diameter,temp,elev,vpd,Rabs,radiative_conductance):
        'this function calculates the humid operative temperature from Campbell and Norman 2010'
        gamma = 0.000666
        gvs = 1/((r_i*100.0)/41.4) #to convert from s/m to m^2 s/mol, divide s/m resistance by the molar volume of air, 41.4 mol m^-3, Campbell and Norman 2010, pg 80
        gva = 1/(((r_b*0.93)*100)/41.4) #rb for water is 0.93 of boundary layer
        gHa = 1/((r_b*100)/41.4)
        gamma_naut = gamma * (1 / gvs + 1 / gva) / (1 / gHa + 1 / radiative_conductance)
        s = ((((17.502*240.97))*0.611*exp((17.502*temp)/(temp+240.97)))/(240.97+temp)**2)/((101.3*exp(-elev/8200)))
        T_eh = temp+(gamma_naut/(gamma_naut+s))*(((Rabs - (self.E_S*(5.670373*10**-8)*((temp+273.15)**4)))/(29.3*(radiative_conductance+gHa)))-(vpd/(gamma_naut*(101.3*exp(-elev/8200)))))
        return T_eh

    def calculate_ea(self,dewpoint):
        'this function calculates the actual vapor pressure'
        ea = ((2.71828182845904**(((1.0/273.0)-(1.0/(dewpoint + 273.15)))*5422.9939))*0.611)
        return ea

    def calculate_es(self,temp_K):
        'this function calculates saturation vapor pressure'
        es = 0.611*exp((L/Rv)*((1./273.15)-(1./temp_K)))
        return es

    def calculate_soil(self,tmax,tmin,hour,depth):
        'this function calculates soil temperatures at the provided depth based on temperatures at the surface, from Campbell and Norman 2010'
        return ((tmax+tmin)/2.0)+((tmax-tmin)/2.0)*(2.71**(-depth/17.0))*sin((3.14/12.)*(hour-8)-depth/17.0)

    def update_rb(self,temp_K,e_s,e_a):
        'this function calculats the boundary layer resistance based upon calculations from Gates 1980 (Biophysical Ecology) and Monteith 2013 (Principles of Environmental Physics)'
        air_pressure = (101325.*(1.-(2.2569*10**-5)*self.elevm)**5.2553)
        air_density = air_pressure/(287.04*temp_K)
        dynamic_viscosity = (1.8325*10**-5)*((296.16+120.)/(temp_K+120.))*((temp_K/296.16)**1.5)
        kinematic_viscosity = dynamic_viscosity/air_density
        T_surface = (temp_K)*(1.+0.38*((e_s*1000.)/air_pressure))
        T_air = (temp_K)*(1.+0.38*((e_a*1000.)/air_pressure))
        coef_thermal_expansion = 1.0/temp_K
        Grashof = (coef_thermal_expansion*gravity*(self.diameter**3)*(abs(T_surface-T_air)))/(kinematic_viscosity**2)
        Nusselt_free = 0.48*((Grashof)**0.25)
        Reynolds_forced = (0.1*self.diameter)/kinematic_viscosity #dimensionless, eq. 9.31 Gates (1980)
        Nusselt_forced = 0.615*Reynolds_forced**0.466 #dimensionless, coefficients from Table 9.3 in Gates (1980)
        thermal_conductivity = (2.4525*10**-2)+((7.038*10**-5)*(temp_K-273.15))
        Nusselt_combined = (Nusselt_free**3 + Nusselt_forced**3)**(1. / 3.) #dimensionless, Bird, Stewart, & Lightfoot's mixed convection formula (p. 445, Transport Phenomena, 2002)
        hc_combined = (Nusselt_combined*thermal_conductivity)/self.diameter #W m^-2 C^-1
        mixing_ratio = (0.6257*(e_a*1000))/(air_pressure-(1.006*(e_a*1000)))
        specific_heat = (1004.84+(1846.4*mixing_ratio))/(1+mixing_ratio)
        self.rb = ((specific_heat * air_density)/hc_combined)/100  # sec m^-1


    def calculate_nightly_activity(self,temp,temp_current,dewpoint,elev,humidity_scenario,depth): #figure out what to do about year and climate scenario
        'calculates total activity time (seconds) and energy balance (in joules) for the average day of the month based on hourly temps, vpds, and soil temps'
        nocturnal_hours = [20,21,22,23,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
        for row in range(len(temp[0])):
            for column in range(len(temp[0][0])):
                for hour in nocturnal_hours:
                    if temp[hour][row][column] == -9999:
                        self.activity_grid.append(-9999)
                        self.energy_grid.append(-9999)
                        break
                    else:
                        pass
                    if hour >= 20 or hour <= 5:
                        if temp[hour][row][column] > 5.0 and temp[hour][row][column] < 28.0:
                            if self.activity_status == 0.0: #0.0 means the salamander is active
                                original_temp = temp[hour][row][column] #units: C
                                temp_K_original = original_temp + 273.15 #units: K
                                e_a_current = float(self.calculate_ea(dewpoint)) #units: kPa
                                e_s_current = float(self.calculate_es(temp_current[hour][row][column] + 273.15)) #units: kPa
                                rh_current = e_a_current/e_s_current #proportion
                                e_s_original = float(self.calculate_es(temp_K_original)) #units: kPa
                                e_a_original = rh_current * e_s_original #units: kPa
                                e_a_new = e_a_original + (e_a_original*humidity_scenario) #units: kPa
                                vpd = (e_s_original-e_a_new) #units: kPa
                                self.elevm = elev[row][column] #units: m
                                self.update_rb(temp_K_original,e_s_original,e_a_new) #units: s/cm
                                self.T_eh = self.calculate_Teh(self.Rs,self.rb,self.diameter,original_temp,self.elevm,vpd,self.Rabs(original_temp),self.radiative_conductance(original_temp)) #units: C
                                if self.acclimation_status == 0.0: #no thermal sensitivity of resistance to water loss
                                    Rs = self.Rs #units: s/cm
                                else: #thermal sensitivity of resistance to water loss
                                    Rs = self.ri_acclimation(self.T_eh) #units: s/cm
                                rho = (e_s_original/(temp_K_original*Rv))-(e_a_new/(temp_K_original*Rv)) #g/cm^3
                                EWL = self.surface_area*(rho/(Rs+(self.rb*0.93))) #adjust boundary layer here for water vapor #g/s
                                self.EWL = EWL
                                hourly_loss = self.EWL*3600.0 #units: g/h
                                if self.mass_to_lose > hourly_loss:
                                    self.mass_to_lose -= hourly_loss
                                    self.activity += 1.0
                                    energy_intake = ((((0.00383582 + (-0.00252254)*self.T_eh + 0.00090893*self.T_eh**2 + (-2.52723e-5)*self.T_eh**3)*1000)*self.mass)/24.) #kJ/g/h
                                    self.energy_status += energy_intake
                                    self.mass_to_lose += (energy_intake/22000.0)*2.33 #g/h of water from food
                                    volume_oxygen = ((((10.0**((0.04618974*self.T_eh)+0.59925591*log10(self.mass)+0.86009768))/1000.0)*1.5)*20.1) #kJ/g/h
                                    self.energy_status -= volume_oxygen
                                    self.Teh_active.append(self.T_eh)
                                    self.Teh_total.append(self.T_eh)
                                else:
                                    self.activity += (self.mass_to_lose/hourly_loss)
                                    energy_intake = ((((0.00383582 + (-0.00252254)*self.T_eh + 0.00090893*self.T_eh**2 + (-2.52723e-5)*self.T_eh**3)*1000)*self.mass)/24.)*(self.mass_to_lose/hourly_loss)
                                    self.energy_status += energy_intake
                                    volume_oxygen = ((((10.0**((0.04618974*self.T_eh)+0.59925591*log10(self.mass)+0.86009768))/1000.0)*1.5)*20.1)*(self.mass_to_lose/hourly_loss)
                                    self.energy_status -= volume_oxygen
                                    self.mass_to_lose = 0.0
                                    self.activity_status += 1.0
                                    self.Teh_active.append(self.T_eh)
                                    self.Teh_total.append(self.T_eh)
                                    tmax = temp[14][row][column]
                                    tmin = temp[6][row][column]
                                    soil_T = self.calculate_soil(tmax,tmin,hour,depth)
                                    volume_oxygen = (((10.0**((0.04618974*soil_T)+0.59925591*log10(self.mass)+0.86009768))/1000.0)*20.1)*(1-(self.mass_to_lose/hourly_loss))
                                    self.energy_status -= volume_oxygen
                                    self.Teh_inactive.append(soil_T)
                                    self.Teh_total.append(soil_T)
                            else:#if individual is inactive
                                tmax = temp[14][row][column]
                                tmin = temp[6][row][column]
                                soil_T = self.calculate_soil(tmax,tmin,hour,depth)
                                volume_oxygen = (((10.0**((0.04618974*soil_T)+0.59925591*log10(self.mass)+0.86009768))/1000.0)*20.1)
                                self.energy_status -= volume_oxygen
                                self.Teh_inactive.append(soil_T)
                                self.Teh_total.append(soil_T)
                        elif temp[hour][row][column] > 28.0 and temp[hour][row][column] <= 33.3: #if temperatures are too warm for activity
                            tmax = temp[14][row][column]
                            tmin = temp[6][row][column]
                            soil_T = self.calculate_soil(tmax,tmin,hour,depth)
                            volume_oxygen = (((10.0**((0.04618974*soil_T)+0.59925591*log10(self.mass)+0.86009768))/1000.0)*20.1)
                            self.energy_status -= volume_oxygen
                            self.Teh_inactive.append(soil_T)
                            self.Teh_total.append(soil_T)
                            pass
                        elif temp[hour][row][column] > 33.3: #CTmax
                            self.energy_status = -9999
                            break
                        else:
                            pass
                    else:#for daytime calculations for energetic costs
                        tmax = temp[14][row][column]
                        tmin = temp[6][row][column]
                        soil_T = self.calculate_soil(tmax,tmin,hour,depth)
                        if self.energy_status == -9999:
                            break
                        elif soil_T > 0.0:
                            volume_oxygen = (((10.0**((0.04618974*soil_T)+0.59925591*log10(self.mass)+0.86009768))/1000.0)*20.1)
                            self.energy_status -= volume_oxygen
                            self.Teh_inactive.append(soil_T)
                            self.Teh_total.append(soil_T)
                        elif soil_T > 33.3: #CTmax
                            self.energy_status = -9999
                            break
                        else:
                            pass
                if temp[hour][row][column] == -9999:
                    pass
                else:
                    self.activity_grid.append(round(self.activity,4)) #after you complete all hours for the month
                    self.energy_grid.append(round(self.energy_status,4))
                    self.activity = 0.0 #reset activity counter between coordinates
                    self.activity_status = 0.0 #reset activity status
                    self.energy_status = 0.0 #reset between coordinates
                    self.mass_to_lose = self.mass - self.activity_threshold #reset mass to loss

    def ri_acclimation(self,temp):
        'this function adjusts the resistance to water loss based upon the air temperature'
        if temp > 24:
            ri = 15.815 + (-1.4928)*24.0 + 0.054181*24.0**2
            return ri
        elif temp >= 12 and temp <= 24:
            ri = 15.815 + (-1.4928)*temp + 0.054181*temp**2
            return ri
        else:
            ri = 15.815 + (-1.4928)*12.0 + 0.054181*12.0**2
            return ri

#-------SCRIPT-----------#
#environmental data#
temperatures = []

#compile temperature data
#working directory: manuscripts
temperatures.append(create_env_lists('green_salamander_project/climate_data/current/*.asc'))
temperatures.append(create_env_lists('green_salamander_project/climate_data/future/cc/*.asc'))
temperatures.append(create_env_lists('green_salamander_project/climate_data/future/hg/*.asc'))

#create elevation values from digital elevation map
dem = open('green_salamander_project/dem_asc.asc')
dem_0 = dem.readlines()
for i in range(6):
    dem_0.pop(0)

elev = []
for line in dem_0:
    line = line.strip()
    rows = line.split()
    elev.append(rows)

for i in range(len(elev)):
    for j in range(len(elev[i])):
        elev[i][j] = float(elev[i][j])

#read coords from captures
IUCN = list(itertools.chain.from_iterable(read_coordinates('green_salamander_project/dem_asc.asc')))
captures = list(itertools.chain.from_iterable(read_coordinates('green_salamander_project/dem_asc.asc')))

#run model for body sizes, skin resistances, epochs, and humidities
list_of_acclimation_status = ['nac','yac']
list_of_resistances = [[6.0,0],[6.0,1]] #[[3.0,0.0],[5.0,0.0],[7.0,0.0]]
list_of_sizes = [2.0,3.0,4.0] #[2.0,3.0,4.0]
list_of_epochs = ["CU","CC","HG"] #["CU","CC","HG"]
list_of_months = [1,2,3,4,5,6,7,8,9,10,11,12]
list_of_thresholds = [0.05,0.075,0.1] # [0.035,0.07,0.1]
list_of_humidities = [0.0,0.25,-0.25] #[0.0,0.25,-0.25]
list_of_depths = [30] #[5.0,20.0]
counter = 0.0
total = len(list_of_resistances) * len(list_of_sizes) * len(list_of_epochs) * len(list_of_thresholds) * len(list_of_humidities) * len(list_of_depths)
for epoch in range(len(temperatures)):
    for threshold in range(len(list_of_thresholds)):
        for body_size in range(len(list_of_sizes)):
            for r_s in range(len(list_of_resistances)):
                print('you are ' + str((counter/total)*100) + '% complete')
                for depth in range(len(list_of_depths)):
                    for humidity in range(len(list_of_humidities)):
                        run_sdm(temperatures,epoch,elev,list_of_sizes,body_size,list_of_resistances,r_s,list_of_thresholds,threshold,list_of_humidities,humidity,list_of_depths,depth)
                        counter += 1.0
print('Your simulation is complete')
```
