# Made by Ivo Cools and Attila Geresdi, Geresdi group, Chalmers. Good of them!
# Contact: cools@chalmers.se and Geresdi@chalmers.se (2022)
# Updated regularly (semi)


# Important with the P5004A: turn on "hislip" from the Network Analyser program. (system - system setup - remote interface)
# Of course also remote drive access. The VISA adress will be your "qcodes"-adress you give in 

import logging
log = logging.getLogger(__name__)

import numpy as np

import pyvisa as visa
import os
from time import sleep

from qcodes import VisaInstrument
from qcodes.instrument.base import Instrument
from qcodes.instrument.channel import InstrumentChannel

from qcodes.utils.validators import Numbers, Enum, Ints
from qcodes.utils.helpers import create_on_off_val_mapping
from qcodes.instrument.parameter import (
    ParamRawDataType
)

from Geresdi_lab_code.lakeshore.Model_372 import Model_372



# Main application code
# Heavily based on
# https://github.com/QCoDeS/Qcodes_contrib_drivers/blob/master/qcodes_contrib_drivers/drivers/CopperMountain/M5180.py


class P5004A(VisaInstrument):
#    P5004_lake = Model_372('lakeshore_372_P5004', 'TCPIP::192.168.0.115::7777::SOCKET')
#    P5004_heater = P5004_lake.sample_heater
    
    def __init__(self, name: str, address: str, timeout: int = 1000, **kwargs):            
        super().__init__(name= name,
                         address= address,
                         timeout= timeout,
                         **kwargs)

        # Args:
        #     name (str): Name you give to the instrument.
        #     address (str): Address of the instrument.
        #     timeout (int): VISA timeout 
        #         for long measurements, potentiall we need to make this longer
        #         default 100000

        self.add_function('reset', call_cmd='*RST')
        
        self.add_parameter(name='output',
                           label='Output',
                           get_parser=str,
                           get_cmd='OUTP:STAT?',
                           set_cmd='OUTP:STAT {}',
                           val_mapping=create_on_off_val_mapping(on_val='1\n',
                                                                 off_val='0\n'))

        self.add_parameter(name='power',
                           label='Power',
                           get_parser=float,
                           get_cmd='SOUR:POW?',
                           set_cmd='SOUR:POW {}',
                           unit='dBm',
                           vals=Numbers(min_value=-80,
                                        max_value=10))

        self.add_parameter(name='if_bandwidth',
                           label='IF Bandwidth',
                           get_parser=float,
                           get_cmd='SENSe:BANDwidth?',
                           set_cmd='SENSe:BANDwidth {}',
                           unit='Hz')

        self.add_parameter('average_amount',
                           label='Averages amount',
                           get_cmd='SENS1:AVER:COUN?',
                           set_cmd='SENS1:AVER ON;:SENS1:AVER:COUN {}',
                           get_parser=int,
                           set_parser=int,
                           unit='',
                           vals=Numbers(min_value=1, max_value=65e3))
    
        #Comm average: set also turns on averaging
        self.add_parameter('averages_enabled',
                           label='Average',
                           get_cmd='SENS1:AVER:STAT?',
                           set_cmd='SENS1:AVER:STAT {}',
                           val_mapping=create_on_off_val_mapping(on_val='1\n',
                                                                 off_val='0\n'))
        
        self.add_parameter('average_on',
                           label='Averages on or off',
                           get_cmd='SENS:AVER?',
                           set_cmd='SENS1:AVER {}',
                           get_parser=str,
                           set_parser=str,
                           unit='',
                           val_mapping=create_on_off_val_mapping(on_val='1\n',
                                                                 off_val='0\n'))
        
        self.add_parameter('average_mode',
                           label='Average mode to sweep or point',
                           get_cmd='SENS:AVER:MODE?',
                           set_cmd='SENS:AVER:MODE {}',
                           get_parser=str,
                           set_parser=str,
                           vals= Enum("POINt", "SWEEP")
                          )
    
            
        self.add_parameter(name='set_Sxx',
                           label='set_Sxx',
                           get_cmd='CALCulate:MEASure:PARameter?',
                           set_cmd='CALCulate:MEASure:PARameter {}',
                           vals=Enum('S11', 'S21', 'S12', 'S22')
                           )

        
        self.add_parameter(name='sweep_type',
                           label='Determine the sweep type/unit',
                           get_cmd='SENS:SWE:TYPE?',
                           set_cmd='SENS:SWE:TYPE {}',
                           vals=Enum('LIN', 'LOG', 'POW','CW', 'SEGM')
                          )

        #####################################################
        self.add_parameter('electrical_delay',
                           label='Electrical delay',
                           get_cmd='CALC:CORRection:EDELay:TIME?',
                           set_cmd='CALC:CORR:EDEL:TIME {}',
                           get_parser=float,
                           set_parser=float,
                           unit='s',
                           vals=Numbers(-10, 10))

        self.add_parameter('electrical_distance',
                           label='Electrical distance',
                           get_cmd='CALC1:CORR:EDEL:DIST?',
                           set_cmd='CALC1:CORR:EDEL:DIST {}',
                           get_parser=float,
                           set_parser=float,
                           unit='m',
                           vals=Numbers())
        self.electrical_distance(53e-9)
        
        self.add_parameter('clock_source',
                           label='Clock source',
                           get_cmd='SENSe:ROSCillator:SOURce?',
                           set_cmd='SENSe:ROSCillator:SOURce {}',
                           get_parser=str,
                           set_parser=str,
                           vals=Enum('int', 'Int', 'INT',
                                     'internal', 'Internal', 'INTERNAL',
                                     'ext', 'Ext', 'EXT',
                                     'external', 'External', 'EXTERNAL'))

        self.add_parameter(name='start_freq',
                           label='Start Frequency',
                           get_parser=float,
                           get_cmd='SENS:FREQ:STAR?',
                           set_cmd=self._set_start,
                           unit='Hz',
                           vals=Numbers(min_value=10e4,
                                        max_value=19.99999e9))

        self.add_parameter(name='stop_freq',
                           label='Stop Frequency',
                           get_parser=float,
                           get_cmd='SENS:FREQ:STOP?',
                           set_cmd=self._set_stop,
                           unit='Hz',
                           vals=Numbers(min_value=10e4 + 1,
                                        max_value=20e9))

        self.add_parameter(name='center_freq',
                           label='Center Frequency',
                           get_parser=float,
                           get_cmd='SENS:FREQ:CENT?',
                           set_cmd=self._set_center,
                           unit='Hz',
                           vals=Numbers(min_value=10e4 + 1,
                                        max_value=19.99999e9))

        self.add_parameter(name='span',
                           label='Frequency Span',
                           get_parser=float,
                           get_cmd='SENS:FREQ:SPAN?',
                           set_cmd=self._set_span,
                           unit='Hz',
                           vals=Numbers(min_value=1,
                                        max_value=19.99999e9))

        self.add_parameter('npts',
                           label='Number of points',
                           get_parser=int,
                           set_parser=int,
                           get_cmd='SENS:SWE:POIN?',
                           set_cmd=self._set_npts,
                           unit='',
                           vals=Ints(min_value=2,
                                     max_value=200001))

        self.add_parameter('traces_amount',
                           label='Number of traces',
                           get_parser=int,
                           set_parser=int,
                           get_cmd='CALC:PAR:COUN?',
                           set_cmd='CALC:PAR:COUN {}',
                           unit='',
                           vals=Ints(min_value=1,
                                     max_value=16))
        #Note: trace amount changing will cancel all ongoing measurements
        
        self.add_parameter(name='trigger_source',
                           label='Trigger source',
                           get_parser=str,
                           get_cmd=self._get_trigger,
                           set_cmd=self._set_trigger,
                           vals=Enum('IMMediate', 'EXTernal', 'MANual'))

        self.add_parameter(name='data_transfer_format',
                           # this can add a speedup to readout
                           label='Data format during transfer',
                           get_parser=str,
                           get_cmd='FORM:DATA?',
                           set_cmd='FORM:DATA {}',
                           vals=Enum('ascii', 'real', 'real32'))
    
        self.add_parameter(name='receiver_gain',
                           label='Sets the gain settings to all ports. Use SENS:SOUR:REC:GAIN:CAT? to return a list of available gain states for the VNA',
                           get_parser=str,
                           get_cmd='SENSe:SOURce:RECeiver:GAIN:ALL?',
                           set_cmd='SENSe:SOURce:RECeiver:GAIN:ALL \'{}\'',
                           vals=Enum('Auto', 'High', 'Low'))
        
        self.connect_message()
        self.startup(address)
        
# --------------------------------------------- #
# ----------------- functions ----------------- #
# --------------------------------------------- #

    def startup(self, address) -> None:
        print("Startup commencing")
        VISApath = r'C:\Program Files\IVI Foundation\VISA\Win64\agvisa\agbin\visa32.dll'
        rm = visa.ResourceManager(VISApath)
        self.device_vna = rm.open_resource(address)
        try:
            print('test')
            # self.device_vna = rm.open_resource('TCPIP0::DHCP2-100217::hislip_PXI10_CHASSIS1_SLOT1_INDEX0::INSTR')
            # Query identification string *IDN?
            self.device_vna.query("*IDN?")
            self.device_vna.query("SYST:ERR?")
            self.device_vna.write("*CLS")

            # start fresh
            self.device_vna.write('CALCulate:PARameter:DELete:ALL')
            sleep(0.25)
            self.device_vna.query("*OPC?")
            self.device_vna.write("CALC:PAR:EXT "'ch1_S21'", 'S21'")
            self.device_vna.query("*OPC?")
            self.device_vna.write("DISP:MEAS:FEED 1")

            print("Startup finished")
        except visa.VisaIOError as e:
            print(f"An error occurred: {e}")
    def connect_lakeshore(self):
        P5004_lake = Model_372('lakeshore_372_P5004', 'TCPIP::192.168.0.115::7777::SOCKET')
        P5004_heater = P5004_lake.sample_heater
        self.P5004_lake = P5004_lake
        self.P5004_heater =  P5004_heater
        
    def disconnect_lakeshore(self):
        P5004_lake = self.P5004_lake
        P5004_lake.close()
        
#      Add perhaps something as a safety under user interaction when cancelling measurement

    def measure_S21(self, 
                    start_freq_or_center: float,
                    stop_freq_or_span: float,
                    points: int, 
                    power: float, 
                    average: int,
                    if_bandwidth: float,
                    type_of_sweep: str,
                    el_delay = 52.90 # in ns
                    ):
           
        ''' 
        The main function. Depending on type_of_sweep, the first and second values are start/stop or center/span resp.
        This value can be either "ss" (start stop) or "cs" (center span)
        we only do S21's, and we only use a single channel at once
        '''

        ##### make sure the window is on
        self.device_vna.write("DISPlay:WINDow1:STATE ON")
        sleep(0.5)
         ##### we name the measurement and make sure it's a S21 measurement
                    # for some reason here we need to use directly talk to the VNA in stead of below, where
                    # it is not necessary. Super weird.. Maybe has to do with the query?
                    # wasn't worth it to check out
              # open a trace if it's not there or not named 
              # hence, only do single trace measurements with this code

        if self.device_vna.query('CALC:PAR:CAT:EXT?') != '"meas,S21"\n':
                # sometimes the name is like "Ch1_S11_1, meas S21".. avoid this by deleting
            self.device_vna.write('CALCulate:PARameter:DELete:ALL')
            self.device_vna.query('*OPC?')
            self.device_vna.write('CALCulate1:PARameter:DEFine:EXT \'meas\',S21')
            self.device_vna.write('*WAI')
            self.device_vna.query('*OPC?')
            self.device_vna.write('CALCulate1:PARameter:SELect \'meas\'')
            self.device_vna.write("*CLS")
            self.device_vna.write('DISPlay:WINDow1:TRACe1:FEED \'meas\'')
    
            # Usually look at data with MLOG unit, most natural
        self.write("CALC:MEAS:FORM MLOG")
            # clear the previous measurement by turning off the averaging - measurement sets averaging to a
            # minimal value
        self.write('SENSe1:AVERage:STATe OFF')
    
        ##### set parameters for sweep
        if type_of_sweep =='ss':
            self.write("SENS1:FREQ:STAR {}".format(start_freq_or_center))
            self._set_stop(stop_freq_or_span)
            # self.write("SENS1:FREQ:STOP {}".format(stop_freq_or_span))
        elif type_of_sweep == 'cs':
            self.write("SENSe1:FREQuency:CENTer {}".format(start_freq_or_center))
            self.write("SENS1:FREQ:SPAN {}".format(stop_freq_or_span))
        else:
            raise Exception('Only allowed sweep types are \'ss\' (start-stop frequency set) and \'cs\' (center-span frequency set)')
            
        self.write('SENSe1:SWEep:POINts {}'.format(points))
        self.write('SENSe1:BANDwidth {}'.format(if_bandwidth))
        self.write('CALCulate1:CORRection:EDELay:TIME {}NS'.format(el_delay))
        self.write('SOUR:POW1 {}'.format(power))
        self.write('SENSe1:SWEep:TIME:AUTO ON') # usually best to have this on
        self.output('on')
        sleep(0.5) # give some time for the command queue to clear
            #If you are looking at the VNAnalyser, then autoscale is always nice
        self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
        
            #averaging should always be > 0, and a natural number at that
        if(average < 1):
            average = 1
        average = round(average//1)
        self.write('SENSe1:AVERage:COUnt {}'.format(average))
        
                ##### start the measurement by starting the averaging
        self.write('SENSe1:AVERage:STATe ON')
        
        
        ##### Wait until the measurement is done. 
        sleep(2.0) # cannot queue measurements without this badboy here, somehow
        self.write("*WAI")
        
            # we check here the status of the averaging every 5 seconds
        sleepcounter = 0
        while self.device_vna.query("STAT:OPER:AVER1:COND?") == '+0\n':
            sleep(2.5)
            sleepcounter += 2.5
            # Let us not average/ rescale too often
            if not sleepcounter % 25:
                self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
            
        cnt = 0
        sleep(1.0)

   
    def save_data(self, 
                  save_path: str, 
                  filename: str, 
                  temperature: float,
                  added_attenuation = 0,
                  vflag = False
                 ):    
        ''' Path does not end in '\' !!
        made for Windows OS
        
        read in the data, and send it to a file
        Usually comes right after the "measure" function
        Filename is just the initial prefix - added text is freq, P, T in final filename
        Added attenuation: att. of 20 dB extra is "-20 dB"
        '''
        
        ###### if temperature is 0, then it means we are underrange of calibration and hence we set the T to 7ish mK
        if temperature == 0.0:
            temperature = 7.0
            
        ##### read in frequency
        start_f = float(self.start_freq())
        stop_f = float(self.stop_freq())
        pts = int(self.npts())
        freq = np.linspace(start_f, stop_f, pts)

        ##### read in REAL
            # may look funny because you literally read what's on the screen
            # I did not manage for now to do it otherwise, but this works... without bugs for now!
        self.write("CALCulate1:MEASure:FORM REAL")
        self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
        self.write("*WAI")
        sleep(1)
        real = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")

        # read in IMAG
        self.write("CALCulate1:MEASure:FORM IMAG")
        self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
        self.write("*WAI")
        sleep(1)
        imag = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")
        
        fullpath = os.path.join(save_path, filename).replace(os.sep, '/')
        # for some nice naming, considering usually we measure in GHz - we round on MHz
        start_f_str = str(round((start_f/1e6)))
        stop_f_str =  str(round((stop_f/1e6)))
        
        # open output file and put data points into the file
        header = 'P = {}dBm \n T = {}mK\n IF_BW = {}Hz, # averages = {}, elec. delay = {} ns \n frequency [Hz], S21 (real), S21 (imaginary), S21 (logmag), S21 (phase)'.format( self.power() + added_attenuation, temperature, round(self.if_bandwidth()), self.average_amount(), round(self.electrical_delay()*1e9, 2) )
        
        if( vflag == True ):
            file = open(fullpath + '.csv',"w")
        if( vflag == False ):
            file = open(fullpath + '_' + str(start_f_str) + 'MHz - ' + str(stop_f_str) + 'MHz_P' + str(self.power() + added_attenuation) + 'dBm' + '_T' + str(temperature) + 'mK' + '.csv',"w")       
        file.write(header + '\n')
        
        S21_mag   = 10*np.log10( np.square(real) + np.square(imag))
        S21_phase = np.arctan(np.divide(real,imag)* 180/np.pi)
        count = 0
        for i in freq:
            file.write( str(i) + ',' + str(real[count]) + ',' + str(imag[count]) + ',' + str(S21_mag[count]) + ',' + str(S21_phase[count]) + '\n')
            count = count + 1
        file.close()
        
        # when we save data, we change the plot appearance, 
        # so we change it back after data saving!
        # it also is more intuitive to look at data this way
        self.write("CALC:MEAS:FORM MLOG")
        self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
        
    def get_data_Vitto(self,
                        variables_ext = []
                    ):    

        ##### read in frequency
        start_f = float(self.start_freq())
        stop_f = float(self.stop_freq())
        pts = int(self.npts())
        freq = np.linspace(start_f, stop_f, pts)

        ##### read in REAL
            # may look funny because you literally read what's on the screen
            # I did not manage for now to do it otherwise, but this works... without bugs for now!
        self.write("CALCulate1:MEASure:FORM REAL")
        self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
        self.write("*WAI")
        sleep(1)
        real = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")

        # read in IMAG
        self.write("CALCulate1:MEASure:FORM IMAG")
        self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
        self.write("*WAI")
        sleep(1)
        imag = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")
                
        S21_mag   = 10*np.log10( np.square(real) + np.square(imag))
        S21_phase = np.arctan(np.divide(real,imag)* 180/np.pi)
        count = 0
        
        data = []
        
        if( len( variables_ext ) != 0 ):
            
            for i in freq:
                data.append( str( variables_ext )[ 1: -1] + ' ,' + str(i) + ' ,' + str(real[count]) + ' ,' + str(imag[count]) + ' ,' + str(S21_mag[count]) + ' ,' + str(S21_phase[count]) )
                count = count + 1
        else:
        
            for i in freq:
                data.append( str(i) + ' ,' + str(real[count]) + ' ,' + str(imag[count]) + ' ,' + str(S21_mag[count]) + ' ,' + str(S21_phase[count]) )
                count = count + 1
        
        # when we get data, we change the plot appearance, 
        # so we change it back after data saving!
        # it also is more intuitive to look at data this way
        self.write("CALC:MEAS:FORM MLOG")
        self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
                
        return data
        # comment    
    def powersweep(self, 
                   start_pwr: float, 
                   end_pwr: float, 
                   amount_of_points: int, 
                   start_freq_or_center: float, 
                   stop_freq_or_span: float,
                   type_of_sweep: str,
                   points: int, 
                   if_bandwidth: float,
                   general_file_path: str, 
                   filename: str, 
                   temperature: float,
                   live_temperature = False,
                   el_delay = 52.90, # in ns
                   added_attenuation = 0,
                   extra_average_factor = 1,
                   user_results_folder = r'\power_sweep_results'
                  ):
        '''
        If you want to do a powersweep, this is your function
        added_attenuation: if you added 20 db Att, then this is "-20"
        Example of user_results_folder: r'\power_sweep_results'
        Live temperature lets you update the temperature for every point; uses channel 6
        '''
        # Power sweep, taking the data at every power point.
        # Adjusted from https://github.com/Boulder-Cryogenic-Quantum-Testbed/measurement/.../self_control/self_control.py
        # First, do temperature stuff:
        if live_temperature:
            P5004_lake = self.P5004_lake
            P5004_lake.ch06.units('kelvin');
            
        ##### create an array with the values of power for each sweep
        sweeps = np.linspace(start_pwr, end_pwr, amount_of_points)
        stepsize = sweeps[0] - sweeps[1]

        ##### create a new directory for the output to be put into
            # make this whatever you want it to be
        if user_results_folder:
            results_folder = user_results_folder
            results_pth = general_file_path + results_folder
        else:
            results_pth = general_file_path
        ##### It should not be a problem that there is already a folder with this name
            # but the user should be told
        try: 
            os.mkdir(results_pth) 
        except OSError as error: 
            print(error)
            print('\n ... but continuing... \n')
        sleep(1)
        
        ##### do the stuff
        itera = 0
        for power_value in sweeps:
            # To be changed how this scales! Probably should depend on the IF_bandwidth
            print("Right now, we are measuring the spectrum at {}dBm applied/inferred power. Value {}/{}.".format(power_value + added_attenuation, itera+1, len(sweeps)) )
          
            average = 5.57846902e-05 * power_value ** 4 + 2.22722094e-03 * power_value ** 3 - 4.88449821e-02 * power_value ** 2 - 2.89754544e+00 * power_value - 1.80165131e+01
            # set to: 
            # 1200 avg at -70 dBm
            # 200      at -40
            # 30       at -20
            # 15       at -10
            
            # with IFBW 1k:
            # 750 avg  at -75 dBm
            # 500      at -70 dBm
            # 50       at -50 dBm
            # 20       at -40 dBm
            # 10       at -30 dBm
            # Use a higher ~6 or so averaging factor if using low-post amplfication measurements
            
            average = extra_average_factor * average
            # Note: due to rounding this can be a bit different than expected
            if average < 5:
                average = 5
            self.measure_S21(start_freq_or_center, stop_freq_or_span, 
                             points, power_value, average, if_bandwidth, type_of_sweep, el_delay)
                             
            if live_temperature:
                temperature = P5004_lake.ch06.temperature()
            self.save_data(results_pth, filename, temperature, added_attenuation)
            sleep(0.5)
            itera += 1
        print("Finished power sweep from {}dBm to {}dBm.".format(start_pwr, end_pwr))
        
        
        
    def powersweep_Vitto(self, 
                   start_pwr: float, 
                   end_pwr: float, 
                   amount_of_points: int, 
                   start_freq_or_center: float, 
                   stop_freq_or_span: float,
                   type_of_sweep: str,
                   points: int, 
                   if_bandwidth: float, 
                   el_delay = 52,
                   added_attenuation = 0,
                   variables_ext = [],
                   silent = False
                  ):
        '''
        If you want to do a powersweep, this is your function
        added_attenuation: if you added 20 db Att, then this is "-20"
        Example of user_results_folder: r'\power_sweep_results'
        Live temperature lets you update the temperature for every point; uses channel 6
        '''
        # Power sweep, taking the data at every power point.
        # Adjusted from https://github.com/Boulder-Cryogenic-Quantum-Testbed/measurement/.../self_control/self_control.py
            
        ##### create an array with the values of power for each sweep
        sweeps = np.linspace(start_pwr, end_pwr, amount_of_points)
        stepsize = sweeps[0] - sweeps[1]
               
        ##### do the stuff
        itera = 0
        
        data = []
        
        for power_value in sweeps:
            # To be changed how this scales! Probably should depend on the IF_bandwidth
          
            #average = 5.57846902e-05 * (power_value + added_attenuation ) ** 4 + 2.22722094e-03 * (power_value + added_attenuation ) ** 3 - 4.88449821e-02 * (power_value + added_attenuation ) ** 2 - 2.89754544e+00 * (power_value + added_attenuation ) - 1.80165131e+01
            average = 5.57846902e-05 * (power_value ) ** 4 + 2.22722094e-03 * (power_value ) ** 3 - 4.88449821e-02 * (power_value ) ** 2 - 2.89754544e+00 * (power_value ) - 1.80165131e+01
            # set to: 
            # 1200 avg at -70 dBm
            # 200      at -40
            # 30       at -20
            # 15       at -10
            
            # with IFBW 1k:
            # 750 avg  at -75 dBm
            # 500      at -70 dBm
            # 50       at -50 dBm
            # 20       at -40 dBm
            # 10       at -30 dBm
            # Use a higher ~6 or so averaging factor if using low-post amplfication measurements
            

            # Note: due to rounding this can be a bit different than expected
            if average < 5:
                average = 5
            self.measure_S21(start_freq_or_center, stop_freq_or_span, 
                             points, power_value, np.round(average), if_bandwidth, type_of_sweep, el_delay)
                             
            data.append( self.get_data_Vitto( variables_ext + [ power_value + added_attenuation, np.round(average) ] ) )
            sleep(0.5)
            itera += 1
        
        if silent == False:
            print("VNA power output sweeped from {}dBm to {}dBm".format(start_pwr, end_pwr ))        
        
        return data
        
    def CW_measurement_UWphase(self, 
                       points: int,
                       center_frequency: float,
                       power_list: list,
                       times_list: list,
                       if_bandwidth: float,
                       save_path: str, 
                       filename: str,
                       temperature: float,
                       el_delay = 60.974,
                       added_attenuation = 0,
                       vflag = False
                      ):
        '''
        Eldelay in ns, standard value is on inclusive room T amplifiers
        added_attenuation: if you added 20 db Att, then this is "-20"
        '''
        fullpath = os.path.join(save_path, filename).replace(os.sep, '/')
        
        self.write('CALCulate1:CORRection:EDELay:TIME {}NS'.format(el_delay))
        BW = round(self.if_bandwidth())
        ED = round(self.electrical_delay()*1e9, 2)
        ## starting stuff 
        self.write("CALC:MEAS:FORM UPHase")
        self.write('SENSe1:AVERage:STATe OFF')
        self.write("SENS:SWE:TYPE CW")
                
        self.write(f"SENSe1:FREQuency:CENTer {center_frequency}")
        self.write(f"SENS:SWE:POIN {points}")
        self.write(f'SENSe1:BANDwidth {if_bandwidth}')
        self.write('DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO')
       
        # reset at end
        self.write("TRIGger:SOURce MANual")
        self.write('SENSe1:SWEep:TIME:AUTO OFF')
        for time_v in times_list:
            self.write(f'SENS:SWE:TIME {time_v}')
            times = np.linspace(0, time_v, points)
            
            for power_v in power_list:
                
                self.write("CALC:MEAS:FORM UPHase")
                self.write('DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO')
                self.write(f'SOUR:POW1 {power_v}')
                sleep(0.5)
                
                self.write("INITiate:IMM")
                #sleep(1) # cannot queue measurements without this badboy here, somehow
                self.device_vna.query("*OPC?")
                sleep(time_v)
                phaseU = self.device_vna.query_ascii_values("CALC:MEAS:DATA:FDATA?")

                ##### read in REAL
                    # may look funny because you literally read what's on the screen
                    # I did not manage for now to do it otherwise, but this works... without bugs for now!
                self.write("CALCulate1:MEASure:FORM REAL")
                self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
                self.write("*WAI")
                sleep(1)
                real = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")

                # read in IMAG
                self.write("CALCulate1:MEASure:FORM IMAG")
                self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
                self.write("*WAI")
                sleep(1)
                imag = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")

                
                header = 'P = {} dBm \n T = {} mK \n IF_BW = {} Hz \n elec. delay = {} ns \n freq_0 = {} Hz \n time [s], S21 (PhaseU), S21 (real), S21 (imaginary)'.format( power_v + added_attenuation, temperature, BW, ED, center_frequency )
                
                                  
                if( vflag == True ):
                    file = open( fullpath + '.csv', "w" )
                if( vflag == False ):
                    file = open( fullpath  + '_P' + str(power_v + added_attenuation) + 'dBm' +'_t' + str(time_v) + 's' + '_T' + str(temperature) + 'mK' + '_f' + str(int(center_frequency)) + 'Hz'+ '.csv',"w")
                file.write(header + '\n')
                
                count = 0
                for i in times:
                    file.write( str(i) + ',' + str(phaseU[count]) + ',' + str(real[count]) + ',' + str(imag[count]) + '\n')
                    count = count + 1
                file.close()
                
                
                
    def CW_measurement_UWphase_Vitto(self, 
                   points: int,
                   center_frequency: float,
                   power_list: list,
                   times_list: list,
                   if_bandwidth: float,
                   el_delay = 60.974,
                  ):
        '''
        Eldelay in ns, standard value is on inclusive room T amplifiers

        '''
        self.write('CALCulate1:CORRection:EDELay:TIME {}NS'.format(el_delay))
        ## starting stuff 
        self.write("CALC:MEAS:FORM UPHase")
        self.write('SENSe1:AVERage:STATe OFF')
        self.write("SENS:SWE:TYPE CW")
                
        self.write(f"SENSe1:FREQuency:CENTer {center_frequency}")
        self.write(f"SENS:SWE:POIN {points}")
        self.write(f'SENSe1:BANDwidth {if_bandwidth}')
        self.write('DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO')

        BW = round(self.if_bandwidth())
        ED = round(self.electrical_delay()*1e9, 2)

        # reset at end
        self.write("TRIGger:SOURce MANual")
        self.write('SENSe1:SWEep:TIME:AUTO OFF')

        for time_v in times_list:
            self.write(f'SENS:SWE:TIME {time_v}')
            times = np.linspace(0, time_v, points)
            up, re, im, uperr, reerr, imerr = [], [], [], [], [], []
            
            for power_v in power_list:
                
                self.write("CALC:MEAS:FORM UPHase")
                self.write('DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO')
                self.write(f'SOUR:POW1 {power_v}')
                sleep(0.5)
                
                self.write("INITiate:IMM")
                #sleep(1) # cannot queue measurements without this badboy here, somehow
                self.device_vna.query("*OPC?")
                sleep(time_v)
                phaseU = self.device_vna.query_ascii_values("CALC:MEAS:DATA:FDATA?")

                ##### read in REAL
                    # may look funny because you literally read what's on the screen
                    # I did not manage for now to do it otherwise, but this works... without bugs for now!
                self.write("CALCulate1:MEASure:FORM REAL")
                self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
                self.write("*WAI")
                sleep(1)
                real = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")

                # read in IMAG
                self.write("CALCulate1:MEASure:FORM IMAG")
                self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
                self.write("*WAI")
                sleep(1)
                imag = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")

                up.append( np.mean( phaseU ) )
                re.append( np.mean( real ) )
                im.append( np.mean( imag ) )
                uperr.append( np.std( phaseU ) )
                reerr.append( np.std( real ) )
                imerr.append( np.std( imag ) )                
            return up, re, im, uperr, reerr, imerr



    def CW_measurement_UWphase_Vitto_2(self, 
                   points: int,
                   center_frequency: float,
                   power_list: list,
                   times_list: list,
                   if_bandwidth: float,
                   el_delay = 60.974
                  ):
        '''
        Eldelay in ns, standard value is on inclusive room T amplifiers
        added_attenuation: if you added 20 db Att, then this is "-20"
        '''
        
        self.write('CALCulate1:CORRection:EDELay:TIME {}NS'.format(el_delay))
        BW = round(self.if_bandwidth())
        ED = round(self.electrical_delay()*1e9, 2)
        ## starting stuff 
        self.write("CALC:MEAS:FORM UPHase")
        self.write('SENSe1:AVERage:STATe OFF')
        self.write("SENS:SWE:TYPE CW")
                
        self.write(f"SENSe1:FREQuency:CENTer {center_frequency}")
        self.write(f"SENS:SWE:POIN {points}")
        self.write(f'SENSe1:BANDwidth {if_bandwidth}')
        self.write('DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO')
       
        # reset at end
        self.write("TRIGger:SOURce MANual")
        self.write('SENSe1:SWEep:TIME:AUTO OFF')
        for time_v in times_list:
            self.write(f'SENS:SWE:TIME {time_v}')
            times = np.linspace(0, time_v, points)
            
            up, re, im, uperr, reerr, imerr = [], [], [], [], [], []
            
            for power_v in power_list:
                
                self.write("CALC:MEAS:FORM UPHase")
                self.write('DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO')
                self.write(f'SOUR:POW1 {power_v}')
                sleep(0.5)
                
                self.write("INITiate:IMM")
                #sleep(1) # cannot queue measurements without this badboy here, somehow
                self.device_vna.query("*OPC?")
                sleep(time_v)
                phaseU = self.device_vna.query_ascii_values("CALC:MEAS:DATA:FDATA?")

                ##### read in REAL
                    # may look funny because you literally read what's on the screen
                    # I did not manage for now to do it otherwise, but this works... without bugs for now!
                self.write("CALCulate1:MEASure:FORM REAL")
                self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
                self.write("*WAI")
                sleep(0.5)
                real = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")

                # read in IMAG
                self.write("CALCulate1:MEASure:FORM IMAG")
                self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
                self.write("*WAI")
                sleep(0.5)
                imag = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")

                up.append( np.mean( phaseU ) )
                re.append( np.mean( real ) )
                im.append( np.mean( imag ) )
                uperr.append( np.std( phaseU ) )
                reerr.append( np.std( real ) )
                imerr.append( np.std( imag ) )                
            return up, re, im, uperr, reerr, imerr

    def CW_measurement_UWphase_Vitto_3(self, 
                   points: int,
                   center_frequency: float,
                   power_list: list,
                   times_list: list,
                   if_bandwidth: float,
                   el_delay = 60.974
                  ):
        '''
        Eldelay in ns, standard value is on inclusive room T amplifiers
        added_attenuation: if you added 20 db Att, then this is "-20"
        '''
        
        self.write('CALCulate1:CORRection:EDELay:TIME {}NS'.format(el_delay))
        BW = round(self.if_bandwidth())
        ED = round(self.electrical_delay()*1e9, 2)
        ## starting stuff 
        self.write("CALC:MEAS:FORM UPHase")
        self.write('SENSe1:AVERage:STATe OFF')
        self.write("SENS:SWE:TYPE CW")
                
        self.write(f"SENSe1:FREQuency:CENTer {center_frequency}")
        self.write(f"SENS:SWE:POIN {points}")
        self.write(f'SENSe1:BANDwidth {if_bandwidth}')
        self.write('DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO')
       
        # reset at end
        self.write("TRIGger:SOURce MANual")
        self.write('SENSe1:SWEep:TIME:AUTO OFF')
        for time_v in times_list:
            self.write(f'SENS:SWE:TIME {time_v}')
            times = np.linspace(0, time_v, points)
            
            up, re, im, uperr, reerr, imerr = [], [], [], [], [], []
            
            for power_v in power_list:
                
                self.write("CALC:MEAS:FORM UPHase")
                self.write('DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO')
                self.write(f'SOUR:POW1 {power_v}')
                #sleep(0.5)
                
                self.write("INITiate:IMM")
                #sleep(1) # cannot queue measurements without this badboy here, somehow
                self.device_vna.query("*OPC?")
                sleep(time_v)
                phaseU = self.device_vna.query_ascii_values("CALC:MEAS:DATA:FDATA?")

                ##### read in REAL
                    # may look funny because you literally read what's on the screen
                    # I did not manage for now to do it otherwise, but this works... without bugs for now!
                self.write("CALCulate1:MEASure:FORM REAL")
                self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
                self.write("*WAI")
                #sleep(0.5)
                real = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")

                # read in IMAG
                self.write("CALCulate1:MEASure:FORM IMAG")
                self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
                self.write("*WAI")
                #sleep(0.5)
                imag = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")

                up.append( np.mean( phaseU ) )
                re.append( np.mean( real ) )
                im.append( np.mean( imag ) )
                uperr.append( np.std( phaseU ) )
                reerr.append( np.std( real ) )
                imerr.append( np.std( imag ) )                
            return up, re, im, uperr, reerr, imerr


    def CW_measurement_UWphase_Vitto_avg(self, 
                   points: int,
                   center_frequency: float,
                   power_list: list,
                   times_list: list,
                   if_bandwidth: float,
                   average: int,
                   el_delay = 60.974
                  ):
        '''
        Eldelay in ns, standard value is on inclusive room T amplifiers
        points: amount on data points in one measurement
        Average: amount of measurements to average
        '''
        self.write('CALCulate1:CORRection:EDELay:TIME {}NS'.format(el_delay))
        ## starting stuff 
        self.write("CALC:MEAS:FORM UPHase")
        self.write('SENSe1:AVERage:STATe OFF')
        self.write("SENS:SWE:TYPE CW")
                
        self.write(f"SENSe1:FREQuency:CENTer {center_frequency}")
        self.write(f"SENS:SWE:POIN {points}")
        self.write(f'SENSe1:BANDwidth {if_bandwidth}')
        self.write('SENSe1:SWEep:TIME:AUTO ON') # usually best to have this on
        self.output('on')
        self.write('DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO')

                
        # reset at end
        self.write("TRIGger:SOURce MANual")
        self.write('SENSe1:SWEep:TIME:AUTO OFF')

        BW = round(self.if_bandwidth())
        ED = round(self.electrical_delay()*1e9, 2)

        if(average < 1):
            average = 1
        average = round(average//1)
        self.write('SENSe1:AVERage:COUnt {}'.format(average))
        
                ##### start the measurement by starting the averaging
        self.write('SENSe1:AVERage:STATe ON')


        ##### Wait until the measurement is done.
        sleep(2.0) # cannot queue measurements without this badboy here, somehow
        self.write("*WAI") 

            # we check here the status of the averaging every 5 seconds
        sleepcounter = 0
        while self.device_vna.query("STAT:OPER:AVER1:COND?") == '+0\n':
            sleep(2.5)
            sleepcounter += 2.5
            print( self.device_vna.query("STAT:OPER:AVER1:COND?") )
            # Let us not average/ rescale too often
            if not sleepcounter % 25:
                self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")

        # reset at end
        self.write("TRIGger:SOURce MANual")
        self.write('SENSe1:SWEep:TIME:AUTO OFF')

        for time_v in times_list:
            self.write(f'SENS:SWE:TIME {time_v}')
            times = np.linspace(0, time_v, points)
            
            up, re, im, uperr, reerr, imerr = [], [], [], [], [], []
            
            for power_v in power_list:               
                
                self.write("CALC:MEAS:FORM UPHase")
                self.write('DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO')
                self.write(f'SOUR:POW1 {power_v}')
                sleep(0.5)       
                                
                self.write("INITiate:IMM")
                

                
                
                        ##### Wait until the measurement is done. 
                sleep(2.0) # cannot queue measurements without this badboy here, somehow
                self.write("*WAI") 

                    # we check here the status of the averaging every 5 seconds
                #sleepcounter = 0
                #while self.device_vna.query("STAT:OPER:AVER1:COND?") == '+0\n':
                 #   sleep(2.5)
                  #  sleepcounter += 2.5
                   # print( self.device_vna.query("STAT:OPER:AVER1:COND?") )
                    # Let us not average/ rescale too often
                    #if not sleepcounter % 25:
                     #   self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
                
                
                
                #sleep(1) # cannot queue measurements without this badboy here, somehow
                self.device_vna.query("*OPC?")
                sleep(time_v)
                phaseU = self.device_vna.query_ascii_values("CALC:MEAS:DATA:FDATA?")

                ##### read in REAL
                    # may look funny because you literally read what's on the screen
                    # I did not manage for now to do it otherwise, but this works... without bugs for now!
                self.write("CALCulate1:MEASure:FORM REAL")
                self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
                self.write("*WAI")
                sleep(1)
                real = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")

                # read in IMAG
                self.write("CALCulate1:MEASure:FORM IMAG")
                self.write("DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO")
                self.write("*WAI")
                sleep(1)
                imag = self.device_vna.query_ascii_values("CALC:MEAS1:DATA:FDATA?")               
                
                up.append( np.mean( phaseU ) )
                re.append( np.mean( real ) )
                im.append( np.mean( imag ) )
                uperr.append( np.std( phaseU ) )
                reerr.append( np.std( real ) )
                imerr.append( np.std( imag ) )                
            return up, re, im, uperr, reerr, imerr         
        
        
    def reset_averages(self) -> None:
            """
            Resets average count to 0
            """
            self.write("SENS1.AVER.CLE")
            
    def _set_start(self, val: float) -> None:
        """Sets the start frequency and updates linear trace parameters.
        Args:
            val (float): start frequency to be set
        Raises:
            ValueError: If start > stop
        """
        stop_freq = self.stop_freq()
        if val >= stop_freq:
            raise ValueError("Stop frequency must be larger than start "
                             "frequency.")
        self.write("SENS1:FREQ:STAR {}".format(val))
        # we get start as the vna may not be able to set it to the
        # exact value provided.
        start_freq = self.start_freq()  #ask for what the VNA has set the starting frequency
        if abs(val - start_freq) >= 1:
            log.warning(
                "Could not set start frequency to {} setting it to "
                "{}".format(val, start_freq)
            )
        

    def _set_stop(self, val: float) -> None:
        """Sets the start frequency and updates linear trace parameters.
        Args:
            val (float): start frequency to be set
        Raises:
            ValueError: If stop < start
        """
        start_freq = self.start_freq()
        if val <= start_freq:
            raise ValueError("Stop frequency must be larger than start "
                             "frequency.")
        self.write("SENS1:FREQ:STOP {}".format(val))
        # We get stop as the vna may not be able to set it to the
        # exact value provided.
        stop_freq = self.stop_freq() #ask to what the VNA has set the end frequency
        if abs(val - stop_freq) >= 1:
            log.warning(
                "Could not set stop frequency to {} setting it to "
                "{}".format(val, stop_freq)
            )
    
    def _set_span(self, val: float) -> None:
        # rather useless
        """Sets frequency span and updates linear trace parameters.
        Args:
            val (float): frequency span to be set
        """
        self.write("SENS1:FREQ:SPAN {}".format(val))

    def _set_center(self, val: float) -> None:
        # rather useless
        """Sets center frequency and updates linear trace parameters.
        Args:
            val (float): center frequency to be set
        """
        self.write("SENS1:FREQ:CENT {}".format(val))
        

    def _set_npts(self, val: int) -> None:
        # rather useless
        """Sets number of points and updates linear trace parameters.
        Args:
            val (int): number of points to be set.
        """
        self.write("SENS1:SWE:POIN {}".format(val))

    def _get_trigger(self) -> str:
        """Gets trigger source.
        Returns:
            str: Trigger source.
        """
        r = self.ask('TRIG:SOUR?')

        if r.lower() == 'man/n': # sends one trigger signal when manually triggered from the front panel or INIT:IMM is sent
            return 'MANual'
        elif r.lower() == 'ext/n': # external (rear panel) source.
            return 'EXTernal'
        elif r.lower() == 'imm/n': # internal source sends continuous trigger signals
            return 'IMMediate'
        else: 
            # we should not get this
            return 'bus'

    def _set_trigger(self, trigger: str) -> None:
        """Sets trigger source.
        Args:
            trigger (str): Trigger source
        """
        self.write('TRIG:SOUR ' + trigger.upper())


    @staticmethod
    def _db(data: np.ndarray) -> np.ndarray:
        """
        Return dB from magnitude
        Args:
            data (np.ndarray): data to be transformed into dB.
        Returns:
            data (np.ndarray): data transformed in dB.
        """

        return 20. * np.log10(np.abs(data))

   
        