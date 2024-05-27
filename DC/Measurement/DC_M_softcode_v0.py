''' DC_M_softcode_v0.py contains functions and classes that can be used to improve DC measurement on 'Ask' 
It is called 'softcode' because it is not meant to control the actual measurement, but to collect all the info we need for that and pass them to the 'hardcode'.
GENERAL RULES:

1) New functions/classes have to be written following a readable spacing/naming strategy...
2) If New funcitons/classis are added or major changes are implemented a new version '..._vx.py' has to be created because we should somehow keep track of the scripts
we used
...

'''

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
''' Here we just import all we basics libraries '''

import os
import csv
import glob
import re
from time import sleep
import pathlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
import math as math
import lmfit
from lmfit import Model, Parameter, report_fit
import pyvisa
from pandas import Series
import scipy
from scipy.interpolate import interp1d
from scipy.fft import fft, fftfreq
from scipy.signal import blackman
from scipy import interpolate

import copy 
from scipy.signal import find_peaks
import statistics
from tabulate import tabulate
from scipy.signal import lfilter

from matplotlib import cm
from matplotlib.ticker import LinearLocator
import Geresdi_lab_code.DC.Measurement.qtplot_data_vVitto as qt
from Geresdi_lab_code.DC.Measurement.qtplot_data_vVitto  import *

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
''' Here we collect in a dictionary the rack configuration as we want them to be saved in the data file. Some of them can be further customize, for example 
with 'dac' names or amplification factor. 

GENERAL RULES:

1)
2)
...

'''
#----------------------------------------------------------------------------------------------------------------------------------------

V_ivvi_config = r'V{}: {} _ S3b {}mV/V.'            
I_ivvi_config = r'I{}: {} _ IVd {}nA/V on.'
Gx1_ivvi_config = r'G{}: {} _ 1V/V.'
Gx5_ivvi_config = r'G{}: {} _ S1f 5V/V.'
Gx30_ivvi_config = r'G{}: {} _ GATE DRIVE 30V/V.'

Vm_ivvi_config_dc = r'Vm: M2b 100V/V dc _ M0a 1kHz Out1 _ Keithley 101.'
Vm_ivvi_config_ac = r'Vm: M2b 100V/V dc + ac {}kV/V_ M0a 1kHz Out1 _ Keithley 101.'

Im_ivvi_config_ln = r'Im: M1h Low Noise x1 {}MV/A on _ M0a 1kHz Out2 _ Keithley 102.'  
Im_ivvi_config_ln_x100dc = r'Im: M1h Low Noise x100dc {}MV/A on _ M0a 1kHz Out2 _ Keithley 102.' 
Im_ivvi_config_ln_x100ac = r'Im: M1h Low Noise x100ac {}MV/A on _ M0a 1kHz Out2 _ Keithley 102.' 


Im_ivvi_config_lrin = r'Im: M1h Low Rin x1 {}MV/A on _ M0a 1kHz Out2 _ Keithley 102.'  
Im_ivvi_config_lrin_x100dc = r'Im: M1h Low Rin x100dc {}MV/A on _ M0a 1kHz Out2 _ Keithley 102.' 
Im_ivvi_config_lrin_x100ac = r'Im: M1h Low Rin x100ac {}MV/A on _ M0a 1kHz Out2 _ Keithley 102.' 


lockin_bias_dV_to_dV = r'dVb: lockin signal output (+V) On Range 1V Offset 0V frequency {}Hz _ S0a In1 1kHz 10mV/V _ S3b (see Vs or Vb).'
lockin_bias_dV_to_dI = r'dIb: lockin signal output (+V) On Range 1V Offset 0V frequency {}Hz _ S0a In1 1kHz 10mV/V _ IVd (see Is or Ib).'

lockin_measure_dV = r'dVm: M2b (see Vm) _ M0a 1kHz Out1 _ lockin signal input +V Range 3 Scaling 1V/V AC _ Low-Pass F. order 8 TC {}s.'
lockin_measure_dI = r'dIm: M1h (see Im) _  M0a 1kHz Out2 _ lockin signal input +V (Range 3 Scaling 1V/V AC) _ Low-Pass F. order 8 TC {}s.'

rack = {
    
    'V': V_ivvi_config,
    'I': I_ivvi_config,
    'Gx1': Gx1_ivvi_config,
    'Gx5': Gx5_ivvi_config,
    'Gx30': Gx30_ivvi_config,
    
    'Vm_dc': Vm_ivvi_config_dc,
    'Vm_ac': Vm_ivvi_config_ac,    
    
    'Im_ln': Im_ivvi_config_ln,
    'Im_ln_100dc': Im_ivvi_config_ln_x100dc,
    'Im_ln_100ac': Im_ivvi_config_ln_x100ac,
    
    'Im_lrin': Im_ivvi_config_lrin,
    'Im_lrin_100dc': Im_ivvi_config_lrin_x100dc,
    'Im_lrin_100ac': Im_ivvi_config_lrin_x100ac,
    
    'dVb': lockin_bias_dV_to_dV,
    'dIb': lockin_bias_dV_to_dI,
    
    'dVm': lockin_measure_dV,
    'dIm': lockin_measure_dI
    
}

G_source_line = lambda  rack_config, s_or_b, dac, c: rack_config.format( s_or_b, dac ) + '\n\t\tp = {}'.format( c )

source_line = lambda  rack_config, s_or_b, dac, mV_or_nA_per_V, c4, c2 : rack_config.format( s_or_b, dac, mV_or_nA_per_V ) + '\n\t\tp4 = {}\n\t\tp2 = {}'.format( c4, c2 )
measure_line = lambda rack_config, kV_or_MV_per_V_or_A, c4, c2 : rack_config.format( kV_or_MV_per_V_or_A ) + '\n\t\tp4 = {}\n\t\tp2 = {}'.format( c4, c2 )  

Vm_line = lambda rack_config, kV_over_V_ac, c4, c2 : rack_config.format( kV_over_V_ac ) + '\n\t\tp4 = {}\n\t\tp2 = {}'.format( c4, c2 )  
Im_line = lambda rack_config, MV_over_A, c4, c2 :  rack_config.format( MV_over_A ) + '\n\t\tp4 = {}\n\t\tp2 = {}'.format( c4, c2 ) 

lockin_measure_line = lambda rack_config, tc: rack_config.format( tc )
lockin_source_line = lambda rack_config, frequency: rack_config.format( frequency )

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------


''' Below we difine some classes that help us collecting the info of our setup and what kind of measurement we want.
We have Device, Gate, Source, Measure, Lock_IN and Lock_OUT. At some point we need to add Magnet and Temperature.'''

        
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------        
        
class Gate:
    
    def __init__( self, rack, dac, sweep_or_bias, sweep_back = False, mV_per_s = 50, unit = 'mV', c1 = '' ):
        
        self.rack0 = rack
        self.dac = dac
        self.c1 = c1        
        self.V_per_V = int( self.line.split('V/V')[ 0 ].split( ' ' )[-1] )
        self.sweep_or_bias = sweep_or_bias
        self.sweep_back = sweep_back
        self.mV_per_s = mV_per_s
        self.unit = unit

        
    
    @property
    def name( self ):
            
            try: 
            
                len( self.sweep_or_bias )
                return 'Gs'
            
            except: return 'Gb'       
    
#    @property
#    def sweep( self ):
#        
#        if self.sweep_back and self.name[ -1 ] == 's':
#            
#            return np.concatenate( [ self.sweep_or_bias, self.sweep_or_bias[ ::-1 ] ] )
#        
#        return self.sweep_or_bias
        
    @property
    def column( self ):
        
        try: size = len( self.sweep_or_bias )
        except: size = 1
        
        do_times = self.V_per_V*1e-3
        
        return '\n\tname: {}\n\tdo_times: {}\n\tsize: {}\n\ttype: coordinate\n'.format( self.name, do_times, size )
    
    @property
    def line( self ):
        
        return G_source_line( self.rack0, self.name[ 1 ], self.dac, self.c1 )    
    
    @property
    def dacn( self ):
        
        return ReadDacNumber( self.dac )
    
    def update_c( self, c1 ):
        
        self.c1 = c1  
        
    def update_sweep_or_bias( self, sweep_or_bias, sweep_back = False ):
        
        self.sweep_or_bias = sweep_or_bias
        self.sweep_back = sweep_back

    def check( self ):
        
        print( self.line )
        print('--------------------')
        print( '\tCOLUMN:' + self.column )
            
#----------------------------------------------------------------------------------------------------------------------------------------    
        
class Source:
    
    def __init__( self, rack, dac, mV_or_nA_per_V, sweep_or_bias, sweep_back = False, mV_per_s = False, unit = 'mV', c4 = '', c2 = '' ):
        
        self.rack0 = rack
        self.dac = dac
        self.sweep_or_bias = sweep_or_bias
        self.mV_or_nA_per_V = mV_or_nA_per_V
        self.sweep_back = sweep_back
        self.mV_per_s = mV_per_s
        self.unit = unit
        self.c4 = c4
        self.c2 = c2

    
    @property
    def name( self ):
            
            try: 
            
                len( self.sweep_or_bias )
                return self.rack0.split( ':' )[ 0 ].format( 's')
            
            except: return self.rack0.split( ':' )[ 0 ].format( 'b')       
    
#    @property
#    def sweep( self ):
#        
#        if self.sweep_back and self.name[ -1 ] == 's':
#            
#            return np.concatenate( [ self.sweep_or_bias, self.sweep_or_bias[ ::-1 ] ] )
#        return self.sweep_or_bias
        
    @property
    def column( self ):
        
        try: size = len( self.sweep_or_bias )
        except: size = 1
        
        if self.name[ 0 ] == 'I':
                
            do_times = 1e-3/self.mV_or_nA_per_V
        
        if self.name[ 0 ] == 'V':
            
            do_times = 1e-3*self.mV_or_nA_per_V
        
        return '\n\tname: {}\n\tdo_times: {}\n\tsize: {}\n\ttype: coordinate\n'.format( self.name, do_times, size )
    
    @property
    def line( self ):
        
        return source_line( self.rack0, self.name[ 1 ], self.dac, self.mV_or_nA_per_V, self.c4, self.c2 )    
    
    @property
    def dacn( self ):
        
        return ReadDacNumber( self.dac )

    def update_c( self, c4, c2 ):
        
        self.c4 = c4
        self.c2 = c2
        
    def update_sweep_or_bias( self, sweep_or_bias, sweep_back = False ):
        
        self.sweep_or_bias = sweep_or_bias
        self.sweep_back = sweep_back
                
    def check( self ):
        
        print( self.line )
        print('--------------------')
        print( '\tCOLUMN:' + self.column )
        
#----------------------------------------------------------------------------------------------------------------------------------------   
        
class Measure:
    
    def __init__( self, rack, kV_or_MV_per_V_or_A, unit = 'V', c4 = '', c2 = '' ):
        
        self.rack0 = rack
        self.kV_or_MV_per_V_or_A = kV_or_MV_per_V_or_A 
        self.unit = unit
        self.c4 = c4
        self.c2 = c2
        self.name = self.rack0.split( ':' )[ 0 ]
    
    @property
    def column( self ):
        
        
        if self.name[ 0 ] == 'I':
            
            additional_amplification = 1
            
            if len( self.rack0.split( 'dc' ) ) == 2:
                
                additional_amplification = 100
                
            do_times = 1/( self.kV_or_MV_per_V_or_A*additional_amplification*1e6 )
        
        if self.name[ 0 ] == 'V':
            
            do_times = 1/100
        
        return '\n\tname: {}\n\tdo_times: {}\n\ttype: value\n'.format( self.name, do_times )
    
    @property
    def line( self ):
        
        return measure_line( self.rack0, self.kV_or_MV_per_V_or_A, self.c4, self.c2 )  
    
    @property
    def keithley( self ):
        
        if '101' in self.line:
            
            return 'k101'
        
        if '102' in self.line:
            
            return 'k102'
        
        if '103' in self.line:
            
            return 'k103'
        

    def update_c( self, c4, c2 ):
        
        self.c4 = c4
        self.c2 = c2
        
    def check( self ):
        
        print( self.line )
        print('--------------------')
        print( '\tCOLUMN:' + self.column )
        
#----------------------------------------------------------------------------------------------------------------------------------------     
        
class Lock_OUT:
    
    def __init__( self, rack, VI_rack_conversion, signal, frequency, unit = 'Vpk' ):
        
        self.rack0 = rack
        self.VI_rack_conversion = VI_rack_conversion
        self.name = self.rack0.split( ':' )[ 0 ]
        self.signal = signal 
        self.frequency = frequency
        self.unit = unit
        self.lockin = 'lockin'
    
    @property
    def column( self ):
            
        do_times = ( 1/100 )*( 1/np.sqrt( 2 ) )*self.VI_rack_conversion
               
        return '\n\tname: {}\n\tdo_times: {}\n\tsize: {}\n\ttype: coordinate\n'.format( self.name, do_times, 1 )
    
    @property
    def line( self ):
        
        return lockin_source_line( self.rack0, self.frequency )          
     
    def update_amplitude( self, amplidute ):
        
        self.signal = amplidute 
        
    def update_frequency( self, frequency ):
        
        self.frequency = frequency 
    
    def check( self ):
        
        print( self.line )
        print('--------------------')
        print( '\tCOLUMN:' + self.column )
        
#----------------------------------------------------------------------------------------------------------------------------------------
        
class Lock_IN:
    
    def __init__( self, rack, ac_amplification, tc, unit = 'V' ):
        
        self.rack0 = rack
        self.ac_amplification = ac_amplification 
        self.unit = unit
        self.tc = tc
        self.name = self.rack0.split( ':' )[ 0 ]
        self.lockin = 'lockin'
    
    @property
    def column( self ):
                     
        do_times = 1/self.ac_amplification
        
        column_re = '\n\tname: {}\n\tdo_times: {}\n\ttype: value\n'.format( self.name + '_re', do_times )
        column_im = '\n\tname: {}\n\tdo_times: {}\n\ttype: value\n'.format( self.name + '_im', do_times )
        
        return column_re, column_im
    
    @property
    def line( self ):
        
        return lockin_measure_line( self.rack0, self.tc )   
    
    def update_tc( self, tv ):
        
        self.tc = tc 
        
    def check( self ):
        
        print( self.line )
        print('--------------------')
        print( '\tCOLUMNs:' + self.column[ 0 ] + self.column[ 1 ] )

#----------------------------------------------------------------------------------------------------------------------------------------

class Device:
    
    def __init__( self, name, config, comments = '' ):
        
        self.name = name
        self.config = config
        self.comments = comments
        
    def update( self, config ):
        
        self.config = config 

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
''' Below we create a SetUp class wchich contains all the info abount the device we are measuring, the tools we are using and the 
rack configuration. 

NOTE: It is important to respect the position of the setup.device attribute when creating the SetUp. the rest of the tools (Sources, Gates, MEasures...)
can be put in any order. However they will be automatically orderes in a way that is consistent with the rest of the code.

'''
    

class SetUp:
    
    def __init__( self, device, *args ):
        
        self.device = device
        
        self.sources = []
        self.gates = []
        self.measures = []
        self.lock_in = []
        self.lock_out = []
        
        for i in [*args]:
            
            if type( i ) == Source:
                
                self.sources.append( i )
                
            if type( i ) == Gate:
                
                self.gates.append( i )
                
            if type( i ) == Measure:
                
                self.measures.append( i )
                
            if type( i ) == Lock_IN: 
                
                self.lock_in.append( i )
                
            if type( i ) == Lock_OUT:
                
                self.lock_out.append( i )
        

        self.tools = self.gates + self.sources + self.measures + self.lock_out + self.lock_in

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
''' Below we collect all the functions we need. Mainly, they concern data file managment and how to write setup info in the the file header. 

GENERAL RULES:

1)
2)
...

'''
#----------------------------------------------------------------------------------------------------------------------------------------


def InsertText( text: str,
                filepath: str ):

#At the given file path in inserts the text at the beginning of the file.
    
    with open( filepath, 'r+') as file:
        content = file.read()
        file.seek(0, 0)
        file.write( text + content )  

#
#----------------------------------------------------------------------------------------------------------------------------------------
#

def ReadDacNumber( dac: str ):
    
#Returns the int number x of the dac form the scring 'dacx'
    
    n = int( dac[ 3:: ] )
    
    return n  

#
#----------------------------------------------------------------------------------------------------------------------------------------
#

def NumberExists( path: str,
                  number: int
                ):

#It checks if the Meas_ number already exists.

    datafilelist = glob.glob( path + '/*' + '.txt' )
    
    for j in datafilelist:

        if ( float( j.split('_')[ -1 ][ : - 4 ] ) == number ):
              
            return True
    
    return False

#
#----------------------------------------------------------------------------------------------------------------------------------------
#

def SetMeasNumber( 
                    path: str
                ):
#It looks for the max Meas_ number in the path and return tha numebr +1 
        
    datafilelist = glob.glob( path + '/*' + '.txt' )
    
    if len( datafilelist ) == 0:
        
        return 0
    
    numbers = [ float( j.split('_')[ -1 ][ : - 4 ] ) for j in datafilelist ]

    return int ( max( numbers ) + 1 )

#
#----------------------------------------------------------------------------------------------------------------------------------------
#

def CreateFilePath( filename: str,
                    filefolder: str,
                    filepath: str,
                    number = False
                  ):

# It looks on the given subfolder path and looks for the max Meas_ numer in the folder and automatically gives a number = number + 1 file
# if number = False. if the file with max Meas_ number has no data inside, it asks you if you want to overwrite that file without creating a new one.
# You can decide to input the Meas_ number after the filepath so that the funciton will not give it automatically. However there if an interlock if that
# number already exixts 

   
    path = os.path.join( filepath, filefolder )
    pathlib.Path( path ).mkdir( parents = True, exist_ok = True )
    
    if not number: 
        
        meas_number = SetMeasNumber( path )
        
    else:
        
        if NumberExists( path, number ):
        
            print( "WARNING: File number already exist!\n" )
        
            YorN = False
        
            while  not YorN:        
            
                flag = input( 'Keep writing on the same file? (y/n)\n')
            
                if flag == 'y':
                    
                    meas_number = number
                    break
            
                if flag == 'n':
                
                    YorN == True
                
                    return -1992
                
        else: meas_number = number
        
    newpath = path + '\\' + filename + '_Meas_{}.txt'.format( meas_number )
    
    
    if NumberExists( path, meas_number -1 ):
        
        try: 
            check = qt.DatFile( path + '\\' + filename + '_Meas_{}.txt'.format( meas_number - 1 ) )
            if len( check.df[ check.df.columns[ 0 ] ] ) == 0:

                YorN = False

                while  not YorN:        

                    flag = input( r'_Meas_{} has no data. Do you want to overwrite? (y/n)'.format( meas_number - 1 ) + '\n' )

                    if flag == 'y':

                        return path + '\\' + filename + '_Meas_{}.txt'.format( meas_number - 1 )
                        break

                    if flag == 'n':

                        return newpath
        except: 
            
            return newpath
    
    return newpath
                
#
#----------------------------------------------------------------------------------------------------------------------------------------
#

def PutLine( 
                filepath: str
            ):
    
#Just write a sepaaretion line of -- in the comment part on the file to help reading. 

    with open( filepath, 'a+') as f:

        header = '\t-------------------------------------------------------------------------------'
        np.savetxt( f, [], header = header, comments = '#' )
 

#
#----------------------------------------------------------------------------------------------------------------------------------------
# 

def WriteInfo( 
                  filepath: str,
                  setup: SetUp
              ):

#This function write the info concerning the SetUp file.  
    
    if filepath == -1:
        
        return 'Error in the path!\n'

    error = WriteDeviceInfo( filepath, setup )
    
    if error != 0:

        print( 'ERROR: WriteDeviceInfo gave error {}'.format( error ) )

    PutLine( filepath )
    
    error = WriteSetUpInfo( filepath, setup )
    
    if error != 0:

        print( 'ERROR: WriteDeviceInfo gave error {}'.format( error ) )    
            
    PutLine( filepath )
    
    error = WriteColumnInfo( filepath, setup )

    if error != 0:

        print( 'ERROR: WriteDeviceInfo gave error {}'.format( error ) )
        
    PutLine( filepath )
    
    error = WriteHeaderInfo( filepath, setup )
    
    if error != 0:

        print( 'ERROR: WriteDeviceInfo gave error {}'.format( error ) )
    
    InsertText( '#\t-------------------------------------------------------------------------------\n', filepath )
#
#----------------------------------------------------------------------------------------------------------------------------------------
#

def WriteDeviceInfo(
                 filepath: str,
                 setup: SetUp
                ):
    
#This function write the info concerning the Device in the SetUp. 
    
    try:
        with open( filepath, 'w+') as f:

            header = '\tDEVICE:'
            np.savetxt( f, [], header = header, comments = '#' )

        with open( filepath, 'a+') as f:

            header = '\t' + setup.device.name
            np.savetxt( f, [], header = header, comments = '#' )
            header = '\t' + setup.device.config
            np.savetxt( f, [], header = header, comments = '#' )
            header = '\t' + setup.device.comments
            np.savetxt( f, [], header = header, comments = '#' )
            
        return 0
    
    except: return -1992
#  
#----------------------------------------------------------------------------------------------------------------------------------------
#

def WriteSetUpInfo(
                 filepath:str,
                 setup: SetUp
                ):
    
#This function write the info concerning rack configuration of the tools in the SetUp. 
    
    try:
        with open( filepath, 'a+') as f:

            header = '\tSETUP:'
            np.savetxt( f, [], header = header, comments = '#' )

        for i in setup.tools :

            with open( filepath, 'a+') as f:

                header = '\t' + i.line
                np.savetxt( f, [], header = header, comments = '#' )
        
        return 0
    
    except: return -1992
            
#
#----------------------------------------------------------------------------------------------------------------------------------------
#

def WriteColumnInfo(
                     filepath:str,
                     setup: SetUp
                    ):

#This function write the column info of the file in a way that the qt analysis code can be used.

    try:
        
        for i in range( 0, len( setup.tools ) ):

            with open( filepath, 'a+') as f:

                if type( setup.tools[ i ] ) == Lock_IN:

                    header = '\tCOLUMN {}:'.format( i ) + setup.tools[ i ].column[ 0 ]
                    np.savetxt( f, [], header = header, comments = '#' )
                    header = '\tCOLUMN {}:'.format( i + 1 ) + setup.tools[ i ].column[ 1 ]
                    np.savetxt( f, [], header = header, comments = '#' )
                    continue

                header = '\tCOLUMN {}:'.format( i ) + setup.tools[ i ].column
                np.savetxt( f, [], header = header, comments = '#' )
        
        return 0
    
    except: return -1992
            
#
#----------------------------------------------------------------------------------------------------------------------------------------
#

def WriteHeaderInfo(
                     filepath:str,
                     setup: SetUp
                    ):
    
# This function write the header of the file, i.e. the name of the columns and their units as they are saved 
    
    column_names, units = [], []
    
    try:
        
        for i in range( 0, len( setup.tools ) ):

            if type( setup.tools[ i ] ) == Lock_IN:

                column_names.append( setup.tools[ i ].name + '_re'  )
                units.append( setup.tools[ i ].unit )
                column_names.append( setup.tools[ i ].name + '_im'  )
                units.append( setup.tools[ i ].unit )
                continue

            column_names.append( setup.tools[ i ].name  )
            units.append( setup.tools[ i ].unit )

        with open( filepath, 'a+') as f:

            header = "\t" + ','.join( [ str( elem ) + "\t" for elem in column_names ] )
            np.savetxt(f, [], header = header, comments = '#' )
            header = "\t" + ','.join( [ str( elem ) + '\t' for elem in units ] )+'\n'
            np.savetxt(f, [], header = header, comments = '#' )

        return 0
    
    except: return -1992

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
'''  -------------------------------- THIS IS THE END ------------------------------------ '''
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
print( 'DC_M_softcode_v0.py lib imported!' )