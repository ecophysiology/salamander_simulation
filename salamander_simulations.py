#---------------------------#
#--------DESCRIPTION--------#
#---------------------------#
'''The following script generates estimates of activity, energy balance, and extinction
throughout the southern Appalachian Mountains for high elevation lungless salamanders 
within the Genus Plethodon. Each section contains a heading that describes the various
components of the script, from libraries to simulations. The script was built using Python
(version 2.7).

The script generates ASCII files and a dataframe of summary statistics. Be sure to enter
the correct path for your temperature data as well as the output paths for the results.
The paths that require the users attention have been labelled as your_path_here.'''

#---------------------------#
#-------- LIBRARIES---------#
#---------------------------#

from numpy import *
from math import *
from random import *
from pandas import *
import glob
import re
import itertools

#---------------------------#
#---------CONSTANTS---------#
#---------------------------#

Rv = 461.5 #J*K^-1*kg^-1
L = 2.5*10**6 #J per kg
gravity = 9.8 #meters per second
header = 'ncols        600\nnrows        480\nxllcorner    -85.000560000000\nyllcorner    34.000560000000\ncellsize     0.008333333333\nNODATA_value -9999\n '

#-----------------------------#
#---------DATAFRAMES----------#
#-----------------------------#

results = pandas.DataFrame(columns = ['year','behavior','scenario','threshold','depth','body_size','skin_resistance','acclimation_status','epoch','humidity','month','mean_activity','sd_activity','IUCN_activity','sd_IUCN_activity','captures_activity','sd_captures_activity','mean_energy','sd_energy','IUCN_energy','sd_IUCN_energy','captures_energy','sd_captures_energy','%_positive_energy','%_positive_IUCN','%_positive_captures','extinct_in_range','Teh_active','Teh_inactive','Teh_total'])
annual_results = pandas.DataFrame(columns = ['year','behavior','scenario','threshold','depth','body_size','skin_resistance','acclimation_status','epoch','humidity','mean_activity','sd_activity','IUCN_activity','sd_IUCN_activity','captures_activity','sd_captures_activity','mean_energy','sd_energy','IUCN_energy','sd_IUCN_energy','captures_energy','sd_captures_energy','%_positive_energy','%_positive_IUCN','%_positive_captures','extinct_in_range','Teh_active','Teh_inactive','Teh_total'])        

#---------------------------#  
#---------FUNCTIONS---------#
#---------------------------#
   
def numericalSort(value):
    'function for reading asc files'
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts
        
def create_env_lists(specific_file):
    'function for reading asc files'
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
    'function for reading asc files'
    list_of_coords = []
    mod1 = open(coordinate_file)
    mod2 = mod1.readlines()
    for k in range(5):
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
    'function for statistics'   
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
    'function for statistics'  
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
    'function for statistics'
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
    'function converting average daily values to annual'
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

def local_extinction(annual_list,lipid_reserve):
    'function for determing energy depletion'
    for i in range(len(annual_list)):
        if annual_list[i] <= lipid_reserve:
            annual_list[i] = -9999
        else:
            pass
    return annual_list

def activity_extinction(activity_list,energy_list):
    'function for statistics'
    for i in range(len(energy_list)):
        if energy_list[i] == -9999:
            activity_list[i] = -9999
        else:
            pass
    return activity_list
    
def calculate_annual_activity(activity_list):
    'function for statistics'
    for i in range(len(activity_list)):
        if activity_list[i] == -9999:
            pass
        else:
            activity_list[i] = activity_list[i]/365.0
    return activity_list

def extinct_in_range(energy_list,range_shape):
    'function for statistics'
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
    'turns list into a map'
    _list = [data_list[i:i+480] for i in range(0, len(data_list), 480)]
    lines = []
    for row in _list:
            lines.append(' '.join(map(str, row)))
    result = '\n '.join(lines)
    f = open(file_name, 'wb')
    f.writelines(str(header))
    f.writelines(str(result))
    f.close()

def calculate_warming(scenario,year):
    'a function to estimate IPCC warming scenarios, scenario 0 is RCP2.6, 1 is RCP4.5, and 2 is RCP8.5'
    if scenario == 0:
        new_temp = (-729.703646061232) + 0.7060042314301908*year + ((-0.0001705273427662923)*year**2)
    elif scenario == 1:
        new_temp = (238.06436957591 + ((-0.2541705332812829)*year) + (6.761168634581547e-5*year**2))
    else:
        new_temp = (528.3628274053061 + ((-0.5551376971298737)*year) + (0.0001454898537077816*year**2))
    return new_temp

def dewpoint(temperature_list):
    'a function to estimate the dewpoint in a list of minimum temperatures'
    minimums = []
    for i in range(len(temperature_list[0])):
        minimums.append(min(temperature_list[10][i]))
    return(min(minimums))
    
def run_sdm(temperatures_list,epoch,elevations,body_size_list,body_size_index,ri_list,ri_index,threshold_list,threshold_index,humidity_list,humidity_index,depth_list,depth_index,behavior_list,behavior_index):    
        'function to run the mechanistic SDM'
        temp_rise = calculate_warming(climate_scenarios[scenario],list_of_years[year])
        _annual_activity = [0] * len(temperatures_list[0][0][0]) * len(temperatures_list[0][0][0][0])
        annual_energy = [0] * len(temperatures_list[0][0][0]) * len(temperatures_list[0][0][0][0])
        lipid_reserves = (1.9855*body_size_list[body_size_index] + 0.1663)*-1000.0 #determine extinction area
        annual_Teh_active = []
        annual_Teh_inactive = []
        annual_Teh_total = []
        for month in range(len(temperatures_list[0])):
            dewpoint_temp = dewpoint(temperatures_list[epoch][month])
            sal = Individual(body_size_list[body_size_index],ri_list[ri_index],threshold_list[threshold_index])
            sal.calculate_nightly_activity(temperatures_list[epoch][month],elevations,temp_rise,dewpoint_temp,depth_list[depth_index],behavior_list[behavior_index],humidity_list[humidity_index])
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
            temporary_dataframe= pandas.DataFrame([[list_of_years[year],behaviors[behavior_status],names_of_scenarios[scenario],list_of_thresholds[threshold],list_of_depths[depth],list_of_sizes[body_size],list_of_resistances[r_s][0],list_of_acclimation_status[r_s],list_of_epochs[epoch],humidity_list[humidity_index],list_of_months[month],mean(sal.activity_grid),std(sal.activity_grid),IUCN_activity[0],IUCN_activity[1],captures_activity[0],captures_activity[1],mean(sal.energy_grid),std(sal.energy_grid),IUCN_energy[0],IUCN_energy[1],captures_energy[0],captures_energy[1],return_positive_energy(sal.energy_grid),return_positive_energy(return_values_from_coords(sal.energy_grid,IUCN)),return_positive_energy(return_values_from_coords(sal.energy_grid,captures)),extinct_in_range(annual_energy,IUCN),Teh_active[0],Teh_inactive[0],Teh_total[0]]],columns = ['year','behavior','scenario','threshold','depth','body_size','skin_resistance','acclimation_status','epoch','humidity','month','mean_activity','sd_activity','IUCN_activity','sd_IUCN_activity','captures_activity','sd_captures_activity','mean_energy','sd_energy','IUCN_energy','sd_IUCN_energy','captures_energy','sd_captures_energy','%_positive_energy','%_positive_IUCN','%_positive_captures','extinct_in_range','Teh_active','Teh_inactive','Teh_total'])
            global results
            results = results.append(temporary_dataframe)
            _annual_activity = convert_to_annual(sal.activity_grid,_annual_activity,month)
            annual_energy = convert_to_annual(sal.energy_grid,annual_energy,month)
            annual_energy = local_extinction(annual_energy,lipid_reserves)
        annual_activity = [_annual_activity[i]/365.0 for i in range(len(_annual_activity))]
        annual_energy_file = 'your_path_here/annual_energy_'+str(list_of_years[year])+'_beh'+str(behaviors[behavior_status])+'_scenario_'+str(names_of_scenarios[scenario])+'_threshold'+str(list_of_thresholds[threshold])+'_depth'+str(list_of_depths[depth])+'_mass'+str(list_of_sizes[body_size])+'_rs'+str(list_of_resistances[r_s][0])+'_acc'+str(list_of_acclimation_status[r_s])+'_epoch'+str(list_of_epochs[epoch])+'_vpd'+str(humidity_list[humidity_index])+'.asc'
        annual_activity_file = 'your_path_here/annual_activity_'+str(list_of_years[year])+'_beh'+str(behaviors[behavior_status])+'_scenario_'+str(names_of_scenarios[scenario])+'_threshold'+str(list_of_thresholds[threshold])+'_depth'+str(list_of_depths[depth])+'_mass'+str(list_of_sizes[body_size])+'_rs'+str(list_of_resistances[r_s][0])+'_acc'+str(list_of_acclimation_status[r_s])+'_epoch'+str(list_of_epochs[epoch])+'_vpd'+str(humidity_list[humidity_index])+'.asc'
        make_and_export_grid(annual_energy,annual_energy_file,header)
        make_and_export_grid(annual_activity,annual_activity_file,header)
        IUCN_annual_activity = return_stats_from_coords(annual_activity,IUCN)
        captures_annual_activity = return_stats_from_coords(annual_activity,captures)
        IUCN_annual_energy = return_stats_from_coords(annual_energy,IUCN)
        captures_annual_energy = return_stats_from_coords(annual_energy,captures)
        temporary_dataframe_annual = pandas.DataFrame([[list_of_years[year],behaviors[behavior_status],names_of_scenarios[scenario],list_of_thresholds[threshold],list_of_depths[depth],list_of_sizes[body_size],list_of_resistances[r_s][0],list_of_acclimation_status[r_s],list_of_epochs[epoch],humidity_list[humidity_index],mean(annual_activity),std(annual_activity),IUCN_annual_activity[0],IUCN_annual_activity[1],captures_annual_activity[0],captures_annual_activity[1],mean(annual_energy),std(annual_energy),IUCN_annual_energy[0],IUCN_annual_energy[1],captures_annual_energy[0],captures_annual_energy[1],return_positive_energy(annual_energy),return_positive_energy(return_values_from_coords(annual_energy,IUCN)),return_positive_energy(return_values_from_coords(annual_energy,captures)),extinct_in_range(annual_energy,IUCN),mean(annual_Teh_active),mean(annual_Teh_inactive),mean(annual_Teh_total)]],columns = ['year','behavior','scenario','threshold','depth','body_size','skin_resistance','acclimation_status','epoch','humidity','mean_activity','sd_activity','IUCN_activity','sd_IUCN_activity','captures_activity','sd_captures_activity','mean_energy','sd_energy','IUCN_energy','sd_IUCN_energy','captures_energy','sd_captures_energy','%_positive_energy','%_positive_IUCN','%_positive_captures','extinct_in_range','Teh_active','Teh_inactive','Teh_total'])
        global annual_results
        annual_results = annual_results.append(temporary_dataframe_annual)
        results.to_csv('your_path_here/climate_data_results.csv')
        annual_results.to_csv('your_path_here/annual_climate_data_results.csv')

        
#---------------------------#
#---------CLASSES-----------#
#---------------------------#
            
class Individual():
    def __init__(self,MASS,r_i,THRESHOLD):#requires mass, skin resistance, and dehydration threshold
        self.mass = MASS
        self.diameter = 0.0016*log(self.mass) + 0.0061 #empirical formula for diameter, Riddell et al. 2017
        self.activity_threshold = self.mass - (THRESHOLD*self.mass)
        self.mass_to_lose = self.mass - self.activity_threshold
        self.Rs = r_i[0] #sec/cm^1
        self.rb = 0.0 #sec/cm^1
        self.activity = 0.0 #hours, keeps track of total activity time
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
        self.E_G = 0.97 # emmisitivity for more soil types (Campbell and Norman 1998)
        self.E_S = 0.96 #emissivity of salamander (Campbell and Norman 1998)
        self.A_L = 0.965 #absorptance of organism to longwave radiation (Bartlett and Gates 1967, Buckley 2008)
    
    def longwave_sky(self,temperature):
        'longwave radiation from sky function, Campbell and Norman 1998'
        return 53.1*10**-14*(temperature+273.15)**6.

    def longwave_ground(self,temperature):
        'longwave radiation from ground function, Campbell and Norman 1998'
        b = 5.670373*10**-8
        return self.E_G*b*(temperature+273.15)**4.
        
    def Rabs(self,temperature):
        'radiation absorbed function, adapted from Campbell and Norman 1998'
        return (0.5*(self.A_L*(self.longwave_sky(temperature)+self.longwave_ground(temperature))))
        
    def radiative_conductance(self,temperature):
        'radiative conductance function, Campbell and Norman 1998'
        return (4.*(5.670373*10**-8)*(temperature+273.15)**3.)/29.3
    
    def calculate_Teh(self,r_i,r_b,diameter,temp,elev,vpd,Rabs,radiative_conductance):
        'calculate humid operative temperature, Campbell and Norman 1998'
        gamma_naut = 0.000666
        a = (r_i*100.0)/41.4
        gva = (r_b*100)/41.4
        rad = (4*5670373*10**(-8)*(temp+273.15)*3.)/29.3
        gamma = gamma_naut*((a+(1./gva))/((1./rad)+(1./gva)))
        s = ((((17.502*240.97))*0.611*exp((17.502*temp)/(temp+240.97)))/(240.97+temp)**2)/(101.3*exp(-elev/8200))
        T_eh = temp+(gamma/(gamma+s))*(((Rabs - (self.E_S*(5.670373*10**-8)*((temp+273.15)**4)))/(29.3*(radiative_conductance+gva)))-(vpd/(gamma*(101.3*exp(-elev/8200)))))
        return T_eh
        
    def calculate_ea(self,dewpoint):
        'calculate actual vapor pressure from dewpoint temperature, Stull 2000'
        ea = ((2.71828182845904**(((1.0/273.0)-(1.0/(dewpoint + 273.15)))*5422.9939))*0.611)
        return ea
        
    def calculate_es(self,temp_K):
        'calculate saturation vapor pressure in kPa, adapted from Stull 2000'
        es = 0.611*exp((L/Rv)*((1./273.15)-(1./temp_K))) 
        return es
        
    def calculate_soil(self,tmax,tmin,hour,depth):
        'a function to estimate soil temperatures; adapted from Campbell and Norman 1998'
        if hour == 0 or hour == 1 or hour == 2 or hour == 3:  
            return ((tmax+tmin)/2.0)+((tmax-tmin)/2.0)*(2.71**(-depth/5.238))*sin((3.14/12.)*(hour-(-13.))-depth/5.238)
        else:
            return ((tmax+tmin)/2.0)+((tmax-tmin)/2.0)*(2.71**(-depth/5.238))*sin((3.14/12.)*(hour-11.)-depth/5.238)
    
    def update_rb(self,temp_K,e_s,e_a):
        'a function to estimate the boundary layer resistance, Riddell et al. 2017'
        air_pressure = (101325.*(1.-(2.2569*10**-5)*self.elevm)**5.2553)
        air_density = air_pressure/(287.04*temp_K)
        dynamic_viscosity = (1.8325*10**-5)*((296.16+120.)/(temp_K+120.))*((temp_K/296.16)**1.5)
        kinematic_viscosity = dynamic_viscosity/air_density
        T_surface = (temp_K)*(1.+0.38*((e_s*1000.)/air_pressure))
        T_air = (temp_K)*(1.+0.38*((e_a*1000.)/air_pressure))
        coef_thermal_expansion = 1.0/temp_K
        Grashof = (coef_thermal_expansion*gravity*(self.diameter**3)*(abs(T_surface-T_air)))/(kinematic_viscosity**2)
        Nusselt = 0.48*((Grashof)**0.25)
        thermal_conductivity = (2.4525*10**-2)+((7.038*10**-5)*(temp_K-273.15))
        hc = (Nusselt*thermal_conductivity)/self.diameter
        mixing_ratio = (0.6257*(e_a*1000))/(air_pressure-(1.006*(e_a*1000)))
        specific_heat = (1004.84+(1846.4*mixing_ratio))/(1+mixing_ratio)
        self.rb = 0.93*((specific_heat*air_density)/hc)/100
         
    def calculate_nightly_activity(self,temp,elev,temperature_rise,dewpoint,depth,behavior_toggle,humidity_scenario): #figure out what to do about year and climate scenario (0.00383582 + (-0.00252254)*self.T_eh + 0.00090893*self.T_eh^2 + (-2.52723e-5)*self.T_eh^3
        'calculates total activity time (seconds) for the week, and a list nightly activity times, based on hourly temps, vpds, and soil temps'
        for row in range(len(temp[0])): #for each row
            for column in range(len(temp[0][0])): #for each column value
                for hour in range(len(temp)):
                    if hour <= 9:
                        new_temp = temp[hour][row][column] + temperature_rise
                        if new_temp > 5.0 and new_temp < 25.0: #range of temperatures for activity
                            if self.activity_status == 0.0: #0.0 means the salamander is active
                                original_temp = temp[hour][row][column]
                                temp_K_new = new_temp + 273.15 #converts Celsius to Kelvin
                                temp_K_original = original_temp + 273.15
                                self.elevm = elev[row][column] #update elevation for rb calculations
                                e_a_original = float(self.calculate_ea(dewpoint))
                                e_s_original = float(self.calculate_es(temp_K_original))
                                e_s_new = float(self.calculate_es(temp_K_new))
                                rh_original = float((e_a_original/e_s_original)*100.0)
                                e_a_new = ((rh_original/100.0) * e_s_new) + (((rh_original/100.0) * e_s_new)*humidity_scenario)
                                vpd = (e_s_original-e_a_original)
                                self.update_rb(temp_K_new,e_s_new,e_a_new) #update boundary layer
                                self.T_eh = self.calculate_Teh(self.Rs,self.rb,self.diameter,new_temp,self.elevm,vpd,self.Rabs(new_temp),self.radiative_conductance(new_temp)) #calculate humid operative body temperature
                                self.Teh_active.append(self.T_eh)
                                self.Teh_total.append(self.T_eh)
                                if self.acclimation_status == 0.0:#a function that helps you to adjust skin resistance based upon VPD
                                    Rs = self.Rs
                                else:
                                    Rs = self.ri_acclimation(vpd)    
                                rho = (e_s_new/(temp_K_new*Rv))-(e_a_new/(temp_K_new*Rv)) #calculate water vapor density gradient
                                EWL = self.surface_area*(rho/(Rs+self.rb))
                                self.EWL = EWL
                                hourly_loss = self.EWL*3600.0
                                energy_intake = ((((0.00383582 + (-0.00252254)*self.T_eh + 0.00090893*self.T_eh**2 + (-2.52723e-5)*self.T_eh**3)*1000)*self.mass)/24.) #joules per hr
                                volume_oxygen = ((((10.0**((0.04094914*self.T_eh)+0.5867291*log10(self.mass)+1.04958053))/1000.0)*1.5)*21.1) #joules per hour
                                behavior = behavior_toggle
                                if behavior > 0: #1 means the salamander exhibits avoidance behavior, 0 means that they will not avoid harsh climates
                                    if energy_intake > volume_oxygen:
                                        if self.mass_to_lose > hourly_loss:
                                            self.mass_to_lose -= hourly_loss #mass to lose is deducted
                                            self.activity += 1.0 #you gain an hour of activity
                                            energy_intake = ((((0.00383582 + (-0.00252254)*self.T_eh + 0.00090893*self.T_eh**2 + (-2.52723e-5)*self.T_eh**3)*1000)*self.mass)/24.) #joules per hr
                                            self.energy_status += energy_intake #gain your joules per hour
                                            self.mass_to_lose += (energy_intake/22000.0)*2.33 #add water from food per joule
                                            volume_oxygen = ((((10.0**((0.04094914*self.T_eh)+0.5867291*log10(self.mass)+1.04958053))/1000.0)*1.5)*21.1) #joules per hour
                                            self.energy_status -= volume_oxygen
                                            self.Teh_active.append(self.T_eh)
                                            self.Teh_total.append(self.T_eh)
                                        else:
                                            self.activity += (self.mass_to_lose/hourly_loss)
                                            energy_intake = ((((0.00383582 + (-0.00252254)*self.T_eh + 0.00090893*self.T_eh**2 + (-2.52723e-5)*self.T_eh**3)*1000)*self.mass)/24.)*(self.mass_to_lose/hourly_loss)
                                            self.energy_status += energy_intake
                                            volume_oxygen = ((((10.0**((0.04094914*self.T_eh)+0.5867291*log10(self.mass)+1.04958053))/1000.0)*1.5)*21.1)*(self.mass_to_lose/hourly_loss) #joules per hour
                                            self.energy_status -= volume_oxygen
                                   	    self.mass_to_lose = 0.0
                                            self.activity_status += 1.0
                                            self.Teh_active.append(self.T_eh)
                                            self.Teh_total.append(self.T_eh)
                                            tmax = temp[19][row][column] + temperature_rise
                                            tmin = temp[10][row][column] + temperature_rise
                                            soil_T = self.calculate_soil(tmax,tmin,hour,depth)
                                            volume_oxygen = (((10.0**((0.04094914*soil_T)+0.5867291*log10(self.mass)+1.04958053))/1000.0)*21.1)*(1-(self.mass_to_lose/hourly_loss)) #joules per hour
                                            self.energy_status -= volume_oxygen
                                            self.Teh_inactive.append(soil_T)
                                            self.Teh_total.append(soil_T)
                                    else:
                                        tmax = temp[19][row][column] + temperature_rise
                                        tmin = temp[10][row][column] + temperature_rise
                                        soil_T = self.calculate_soil(tmax,tmin,hour,depth)
                                        volume_oxygen = (((10.0**((0.04094914*soil_T)+0.5867291*log10(self.mass)+1.04958053))/1000.0)*21.1) #joules per hour
                                        self.energy_status -= volume_oxygen
                                        self.Teh_inactive.append(soil_T)
                                        self.Teh_total.append(soil_T)
                                else:
                                    if self.mass_to_lose > hourly_loss:
                                            self.mass_to_lose -= hourly_loss #mass to lose is deducted
                                            self.activity += 1.0 #you gain an hour of activity
                                            energy_intake = ((((0.00383582 + (-0.00252254)*self.T_eh + 0.00090893*self.T_eh**2 + (-2.52723e-5)*self.T_eh**3)*1000)*self.mass)/24.) #joules per hr
                                            self.energy_status += energy_intake #gain your joules per hour
                                            self.mass_to_lose += (energy_intake/22000.0)*2.33 #add water from food per joule
                                            volume_oxygen = ((((10.0**((0.04094914*self.T_eh)+0.5867291*log10(self.mass)+1.04958053))/1000.0)*1.5)*21.1) #joules per hour
                                            self.energy_status -= volume_oxygen
                                            self.Teh_active.append(self.T_eh)
                                            self.Teh_total.append(self.T_eh)
                                    else:
                                            self.activity += (self.mass_to_lose/hourly_loss)
                                            energy_intake = ((((0.00383582 + (-0.00252254)*self.T_eh + 0.00090893*self.T_eh**2 + (-2.52723e-5)*self.T_eh**3)*1000)*self.mass)/24.)*(self.mass_to_lose/hourly_loss)
                                            self.energy_status += energy_intake
                                            volume_oxygen = ((((10.0**((0.04094914*self.T_eh)+0.5867291*log10(self.mass)+1.04958053))/1000.0)*1.5)*21.1)*(self.mass_to_lose/hourly_loss) #joules per hour
                                            self.energy_status -= volume_oxygen
                                   	    self.mass_to_lose = 0.0
                                            self.activity_status += 1.0
                                            self.Teh_active.append(self.T_eh)
                                            self.Teh_total.append(self.T_eh) 
                                            tmax = temp[19][row][column] + temperature_rise
                                            tmin = temp[10][row][column] + temperature_rise
                                            soil_T = self.calculate_soil(tmax,tmin,hour,depth)
                                            volume_oxygen = (((10.0**((0.04094914*soil_T)+0.5867291*log10(self.mass)+1.04958053))/1000.0)*21.1)*(1-(self.mass_to_lose/hourly_loss)) #joules per hour
                                            self.energy_status -= volume_oxygen
                                            self.Teh_inactive.append(soil_T)
                                            self.Teh_total.append(soil_T)
                            else:#if individual is inactive
                                tmax = temp[19][row][column] + temperature_rise
                                tmin = temp[10][row][column] + temperature_rise
                                soil_T = self.calculate_soil(tmax,tmin,hour,depth)
                                volume_oxygen = (((10.0**((0.04094914*soil_T)+0.5867291*log10(self.mass)+1.04958053))/1000.0)*21.1) #joules per hour
                                self.energy_status -= volume_oxygen
                                self.Teh_inactive.append(soil_T)
                                self.Teh_total.append(soil_T)
                        elif new_temp > 25.0 and new_temp <= 34.9: #if temperatures are too warm for activity
                            tmax = temp[19][row][column] + temperature_rise
                            tmin = temp[10][row][column] + temperature_rise
                            soil_T = self.calculate_soil(tmax,tmin,hour,depth)
                            volume_oxygen = (((10.0**((0.04094914*soil_T)+0.5867291*log10(self.mass)+1.04958053))/1000.0)*21.1) #joules per hour
                            self.energy_status -= volume_oxygen
                            self.Teh_inactive.append(soil_T)
                            self.Teh_total.append(soil_T)
                        elif new_temp > 34.9: #CTmax
                            self.energy_status = -9999
                            break
                        else:
                            pass
                    else:#this is for daytime calculations for energetic costs
                        tmax = temp[19][row][column] + temperature_rise
                        tmin = temp[10][row][column] + temperature_rise
                        soil_T = self.calculate_soil(tmax,tmin,hour,depth)
                        if self.energy_status == -9999:
                            break
                        elif soil_T > 0.0:
                            volume_oxygen = (((10.0**((0.04094914*soil_T)+0.5867291*log10(self.mass)+1.04958053))/1000.0)*21.1) #joules per hour
                            self.energy_status -= volume_oxygen
                            self.Teh_inactive.append(soil_T)
                            self.Teh_total.append(soil_T)
                        elif soil_T > 34.9: #CTmax
                            self.energy_status = -9999
                            break
                        else:
                            pass
                self.activity_grid.append(round(self.activity,4)) #after you complete all hours for the month
                self.energy_grid.append(round(self.energy_status,4))
                self.activity = 0.0 #reset activity counter between coordinates
                self.activity_status = 0.0 #reset activity status
                self.energy_status = 0.0 #reset between coordinates
                self.mass_to_lose = self.mass - self.activity_threshold #reset mass to lose
    
    def ri_acclimation(self,vpd):
        if vpd > 0.12:
            ri = 7.0
            return ri
        elif vpd < 0.12 and vpd > 0.02:
            ri = (vpd*11.168) + 5.302
            return ri
        else:
            ri = 5.0
            return ri 
            
#-------------------------------#                     
#---------ENVIRONMENT-----------#
#-------------------------------# 

'''The following script requires a list of average hourly temperature data
for each month of the year in ASCII (.asc) format as well as a digital
elevation map of the same resolution and extent. We also analyzed the proportion
of sites with positive energy balance based upon the geographic distribution of
the Plethodon jordani species complex from the IUCN and the all historic captures
of these species using rasters. We will provide these maps, but you must point to 
locations of the maps on your system.'''

temperatures = []
temperatures.append(create_env_lists('your_path_here/*.asc'))


#create elevation values from DEM
dem = open('your_path_here/dem_for_rb.asc')
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
        
#read coordinates of IUCN geographic range and historic capture coordinates for statistics
IUCN = list(itertools.chain.from_iterable(read_coordinates('your_path_here/raster.asc')))
captures = list(itertools.chain.from_iterable(read_coordinates('your_path_here/points.asc')))

#-------------------------------#                     
#---------SIMULATION------------#
#-------------------------------# 

list_of_acclimation_status = ['nac','nac']#nac indicates constant value of skin resistance, yac will adjust skin resitance value based upon VPD based upon empirical evidence of acclimation
list_of_resistances = [[4.0,0.0],[7.0,0.0]]#must be entered as a list of lists, first value indicates skin resistances, second value indicates whether skin resistance is constant (0) or varies with VPD (1)
list_of_sizes = [2.0,3.0,4.0,5.0]#run simulations for all of these masses (g)
list_of_epochs = ["CC"]
list_of_months = [1,2,3,4,5,6,7,8,9,10,11,12] #run simulations for all of these months of the year
list_of_thresholds = [0.035,0.07,0.1] #run simulations for all of these dehydration thresholds
list_of_years = range(2000,2101,10) #run simulations for all of these years
climate_scenarios = [1,2]#run simulations for all of these climate change scenarios
list_of_depths = [2.5,30.0] #run simulations for all of these soil depths of inactive salamanders
behaviors = [0,1]#run simulations for all of these avoidance behaviors (0 = no avoidance; 1 = avoidance)
list_of_humidities = [0.0]#run simulations for all of these humidities (varied by percent from normal)
names_of_scenarios = ['rcp85'] #['rcp26','rcp45','rcp85']
results = pandas.DataFrame(columns = ['year','behavior','scenario','threshold','depth','body_size','skin_resistance','acclimation_status','epoch','humidity','month','mean_activity','sd_activity','IUCN_activity','sd_IUCN_activity','captures_activity','sd_captures_activity','mean_energy','sd_energy','IUCN_energy','sd_IUCN_energy','captures_energy','sd_captures_energy','%_positive_energy','%_positive_IUCN','%_positive_captures','extinct_in_range','Teh_active','Teh_inactive','Teh_total'])
annual_results = pandas.DataFrame(columns = ['year','behavior','scenario','threshold','depth','body_size','skin_resistance','acclimation_status','epoch','humidity','mean_activity','sd_activity','IUCN_activity','sd_IUCN_activity','captures_activity','sd_captures_activity','mean_energy','sd_energy','IUCN_energy','sd_IUCN_energy','captures_energy','sd_captures_energy','%_positive_energy','%_positive_IUCN','%_positive_captures','extinct_in_range','Teh_active','Teh_inactive','Teh_total'])
for epoch in range(len(temperatures)):
    for behavior_status in range(len(behaviors)):
        for threshold in range(len(list_of_thresholds)):
            for body_size in range(len(list_of_sizes)):
                for r_s in range(len(list_of_resistances)):
                    for scenario in range(len(climate_scenarios)):
                        for depth in range(len(list_of_depths)):
                            for year in range(len(list_of_years)):
                                for humidity in range(len(list_of_humidities)):
                                    run_sdm(temperatures,epoch,elev,list_of_sizes,body_size,list_of_resistances,r_s,list_of_thresholds,threshold,list_of_humidities,humidity,list_of_depths,depth,behaviors,behavior_status)
print 'Your simulation is complete'                                  