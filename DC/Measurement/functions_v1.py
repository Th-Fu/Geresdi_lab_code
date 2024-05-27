import os
import pathlib
import glob
from time import sleep
import time
import matplotlib.pyplot as plt
import numpy as np
import qcodes as qc
import math as math
from qcodes import (
    Measurement,
    experiments,
    initialise_database,
    initialise_or_create_database_at,
    load_by_guid,
    load_by_run_spec,
    load_experiment,
    load_last_experiment,
    load_or_create_experiment,
    new_experiment,
)
from qcodes.dataset.plotting import plot_dataset
from qcodes.logger.logger import start_all_logging
from qcodes.tests.instrument_mocks import DummyInstrument
import qcodes.instrument_drivers.tektronix.Keithley_6500
import qcodes_contrib_drivers.drivers.QuTech.IVVI
import qcodes.instrument_drivers.american_magnetics.AMI430
import zhinst.qcodes
import datetime

import lmfit
from pandas import Series

from lmfit import Model, Parameter, report_fit

from qcodes.instrument_drivers.Lakeshore.Model_372 import Model_372


def InsertText( text: str,
                filepath: str ):
    
    with open( filepath, 'r+') as file:
        content = file.read()
        file.seek(0, 0)
        file.write( text + content )  
        
def ReadDacNumber( dac: str ):
    
    n = int( dac[ 3: ] )
    
    return n  
        
def NumberExists( path: str,
                  number: int
                ):
        
    datafilelist = glob.glob( path + '/*' + '.txt' )
    
    for j in datafilelist:

        if ( float( j.split('_')[ -1 ][ : - 4 ] ) == number ):
              
            return True
    
    return False
    

def CreateDataFile( filename: str,
                    filefolder: str,
                    filepath: str,
                    number: int
                  ):
   
    path = os.path.join( filepath, filefolder )
    pathlib.Path( path ).mkdir( parents = True, exist_ok = True )
    
    if NumberExists( path, number ):
        
        print( "WARNING: File number already exist!\n" )
        
        YorN = False
        
        while  not YorN:        
            
            flag = input( 'Keep writing on the same file? (y/n)\n')
            
            if flag == 'y':
                
                break
            
            if flag == 'n':
                
                YorN == True
                
                return -1
        
    newpath = path + '\\' + filename + '_Meas_{}.txt'.format( number )
    
    return newpath

def WriteSetupInfo( 
              *args: str
            ):
    
   
    setup = ['SETUP']
        
    for i in args:
        
        setup.append( '\t' + i )
    
    return setup


def PutLine( filepath:str ):
    
    with open( filepath, 'a+') as f:

        header = '\t-------------------------------------------------------------------------------'
        np.savetxt( f, [], header = header, comments = '#' )

def CreateHeader( 
                  columns_name: list,
                  amplifications: list,
                  units: list,
                  filepath: str,
                  setupinfo: list,
                 ):

    if filepath == -1:
        
        return 'Overwriting aborted\n'
    
    PutLine( filepath )
    
    for i in setupinfo:

        with open( filepath, 'a+') as f:

            header = '\t' + i
            np.savetxt( f, [], header = header, comments = '#' )
            
    PutLine( filepath )
        
    for i in range( 0, len( columns_name ) ):

        with open( filepath, 'a+') as f:

            header = '\tCOLUMN {}:\n\tname: {}\n\tdo_times: {}\n'.format( str( i ), columns_name[ i ], amplifications[ i ] )
            np.savetxt( f, [], header = header, comments = '#' )
    
    PutLine( filepath )
    
    with open( filepath, 'a+') as f:

        header = "\t" + ','.join( [ str( elem ) + "\t" for elem in columns_name ] ) + '\n' + "\t" + ','.join( [ str( elem ) + '\t' for elem in units ] )+'\n'
        np.savetxt(f, [], header = header, comments = '#' ) 


class SM_Tool:

    def __init__( self, name, do_times, tool_info, type_sweep ):
        
        self.name = name
        self.do_times = do_times
        self.tool_info = tool_info
        self.type_sweep = type_sweep

    def BackSweep( self ):
        
        return SM_tool( name, do_times, tool_info, type_sweep[ ::-1 ] )

class SetUp:

    def __init__( self, lake, ivvi , kV, kI, kX, lockin, mag_x, mag_y, mag_z ):
        
        self.lake = lake
        self.ivvi = ivvi 
        self.kV = kV
        self.kI = kI
        self.kX = kX 
        self.lockin = lockin
        self.mag_x = mag_x
        self_mag_y = mag_y
        self_mag_z = mag_z
        
        
    def RampGate(
          self,
          dac: SM_Tool,
          set_point: float,
          mV_per_second = 50,
          rounding = 0
       ):
       #It ramps the input dac from the current value to the set_point with a speed of (50/amp)-ish mV/s.
       #You can change the speed but you can't go faster than 83-ish mV/s.

        dac_number = ReadDacNumber( dac.name )

        starting_point = round( self.ivvi.dac_voltages()[ dac_number - 1 ], rounding  )

        if mV_per_second < 83:

            points = abs( set_point/dac.do_times - starting_point )/( mV_per_second/dac.do_times )

            if int( points ) == 0:

                self.ivvi.set( dac.name, set_point/dac.do_times )
                return '{} is already there'.format( dac.name )

            sweep = np.linspace( starting_point, set_point, int( points ) )

            for i in sweep:

                self.ivvi.set( dac.name, i/amp.do_times )

                sleep( 1 )

        if mV_per_second >= 83:

            self.ivvi.set( dac.name, set_point/dac.do_times )

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def Plot_Vs_Im( 
                    self,
                    filepath: str,
                    s_dac: SM_Tool,
                    m_k: SM_Tool,
                    sweepback: bool,
                    fig = False,
                    sleep_time = 0.1
                    pcolor = 'bo',
                    dac_factor = 1e-3
                    ):


        dacnumber = ReadDacNumber( s_dac.name )

        for i in s_dac.sweep:

            self.ivvi.set( Vs_dac.name, i/s_dac.do_times )

            time.sleep( sleep_time )

            V_s = self.ivvi.dac_voltages()[ dacnumber - 1 ]

            I_m = self.kI.amplitude() 
           
            plt.plot( V_s*dac_factor*s_dac.do_times, I_m*m_k.do_times, pcolor )
            plt.xlabel( 'V_s (V)' )
            plt.ylabel( 'I_m (A)' )
            plt.ticklabel_format( scilimits=(-3, 4) )
            fig.canvas.draw()

            if filepath == '':

                continue

            data = np.column_stack( ( V_s, I_m ) )

            with open( path, 'a+' ) as f:

                np.savetxt( f, data, delimiter = '\t' ) 
                f.flush()   

        if sweepback:
        
            self.Sweep_Vs_Im( filepath, s_dac.BackSweep(), m_k, False, sleep_time, fig, 'ro' )
            
    
    def Plot_Is_Vm( 
                    self,
                    filepath: str,
                    s_dac: SM_Tool,
                    m_k: SM_Tool,
                    sweepback: bool,
                    fig = False,
                    sleep_time = 0.1
                    pcolor = 'bo',
                    dac_factor = 1e-3
                    ):


        dacnumber = ReadDacNumber( s_dac.name )

        for i in s_dac.sweep:

            self.ivvi.set( s_dac.name, i/s_dac.do_times )

            time.sleep( sleep_time )

            I_s = self.ivvi.dac_voltages()[ dacnumber - 1 ]

            V_m = self.kV.amplitude() 
           
            plt.plot( I_s*dac_factor*s_dac.do_times, V_m*m_k.do_times, pcolor )
            plt.xlabel( 'I_s (A)' )
            plt.ylabel( 'V_m (V)' )
            plt.ticklabel_format( scilimits=(-3, 4) )
            fig.canvas.draw()

            if filepath == '':

                continue

            data = np.column_stack( ( I_s, V_m ) )

            with open( path, 'a+' ) as f:

                np.savetxt( f, data, delimiter = '\t' ) 
                f.flush()   

        if sweepback:
        
            self.Plot_Is_Vm( filepath, s_dac.BackSweep(), m_k, False, sleep_time, fig, 'ro' )
            
    
    def Sweep_Vs_Im( 
                    self,
                    filepath: str,
                    s_dac: SM_Tool,
                    m_k: SM_Tool,
                    sweepback: bool,
                    sweep_time = 0.1,
                    save_more = 0
                    ):


        dacnumber = ReadDacNumber( s_dac.name )

        for i in Vs_dac.sweep:

            self.ivvi.set( s_dac.name, i/s_dac.do_times )
            
            time.sleep( sleep_time )
            
            V_s = self.ivvi.dac_voltages()[ dacnumber - 1 ]

            I_m = self.kI.amplitude() 
           
            if save_more == 0:
                
                data = np.column_stack( ( V_s, I_m ) )
                
            if save_more != 0:
           
                data = np.column_stack( [ V_s, I_m ] + save_more )

            with open( path, 'a+' ) as f:

                np.savetxt( f, data, delimiter = '\t' ) 
                f.flush()   

        if sweepback:
        
            self.Sweep_Vs_Im( filepath, s_dac.BackSweep(), m_k, False, sleep_time, save_more )
            
            
    def Sweep_Is_Vm( 
                    self,
                    filepath: str,
                    s_dac: SM_Tool,
                    m_k: SM_Tool,
                    sweepback: bool,
                    sleep_time = 0.1
                    dac_factor = 1e-3
                    ):


        dacnumber = ReadDacNumber( s_dac.name )

        for i in s_dac.sweep:

            self.ivvi.set( s_dac.name, i/s_dac.do_times )
            
            time.sleep( sleep_time )

            I_s = self.ivvi.dac_voltages()[ dacnumber - 1 ]

            V_m = self.kV.amplitude() 

            if save_more == 0:
                
                data = np.column_stack( ( V_s, I_m ) )
                
            if save_more != 0:
           
                data = np.column_stack( [ V_s, I_m ] + save_more )

            with open( path, 'a+' ) as f:

                np.savetxt( f, data, delimiter = '\t' ) 
                f.flush()   

        if sweepback:
        
            self.Sweep_Is_Vm( filepath, s_dac.BackSweep(), m_k, False, sleep_time, save_more )


    