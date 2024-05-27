
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



class SetUp:

    def __init__( self, lake, ivvi , kI, kV, kIG, lockin, mag_x, mag_y, mag_z ):
        
        self.lake = lake
        self.ivvi = ivvi 
        self.kI = kI
        self.kV = kV
        self.kIG = kIG 
        self.lockin = lockin
        self.mag_x = mag_x
        self_mag_y = mag_y
        self_mag_z = mag_z
        
        
    def RampGate(
          self,
          dac: str,
          set_point: float,
          amp: float,
          mV_per_second = 50,
          rounding = 0
       ):
       #It ramps the input dac from the current value to the set_point with a speed of (50/amp)-ish mV/s.
       #You can change the speed but you can't go faster than 83-ish mV/s.

        dac_number = ReadDacNumber( dac )

        starting_point = round( self.ivvi.dac_voltages()[ dac_number - 1 ], rounding  )

        if mV_per_second < 83:

            points = abs( set_point/amp - starting_point )/( mV_per_second/amp )

            if int( points ) == 0:

                self.ivvi.set( dac, set_point/amp )
                return '{} is already there'.format( dac )

            sweep = np.linspace( starting_point, set_point, int( points ) )

            for i in sweep:

                self.ivvi.set( dac, i/amp )

                sleep( 1 )

        if mV_per_second >= 83:

            self.ivvi.set( dac, set_point/amp )

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def Sweep_Vs_Im( 
                    self,
                    filepath: str,
                    dac: str,
                    V_over_V_s: float,
                    V_over_A_m: float,
                    sweep: list,
                    sweepback: bool,
                    fig = False,
                    pcolor = 'bo',
                    dac_factor = 1e-3
                    ):


        dacnumber = ReadDacNumber( dac )

        for i in sweep:

            self.ivvi.set( dac, i )

            V_s = self.ivvi.dac_voltages()[ dacnumber - 1 ]

            I_m = self.kI.amplitude() 
           
            plt.plot( V_s*dac_factor*V_over_V_s, I_m*( 1/V_over_A_m ), pcolor )
            plt.xlabel('V_s (V)')
            plt.ylabel('I_m (A)')
            plt.ticklabel_format( scilimits=(-3, 4) )
            fig.canvas.draw()

            time.sleep( 0.1 )

            if filepath == '':

                continue

            data = np.column_stack( ( V_s, I_m ) )

            with open( path, 'a+' ) as f:

                np.savetxt( f, data, delimiter = '\t' ) 
                f.flush()   

        if sweepback:
        
            self.Measurement_Vs_Im( filepath, dac, V_over_V_s, V_over_A_m, sweep[ :: -1 ], False, fig, 'ro' )


    def Sweep_Vg_Vb_Im_Vm(
                        filepath: str,
                        V_g_dac: str,
                        V_b_dac:str,
                        V_b: float,
                        V_g_sweep: list,
                        V_over_V_g: float,
                        V_overo_V_m: float,
                        V_over_A_m: float,
                        fig = False,
                        back = True,
                        V_g_ramping_mV_per_second = 50, 
                        pcolor = 'bo',
                        dac_factor = 1e-3
                    ):

        V_b_dac_number = ReadDacNumber( V_b_dac )
        V_g_dac_number = ReadDacNumber( V_g_dac )
        
        print( 'ramping dacs...')
        
        self.ivvi.set( V_b_dac , V_b )
        RampGate( V_g_dac, V_g_sweep[ 0 ], V_over_V_g, V_g_ramping_mV_per_second )
        
        print( 'dacs at target' )
                
        for i in V_g_sweep:

            self.ivvi.set( V_g_dac, i/V_over_V_g )

            V_b = self.ivvi.dac_voltages()[ V_b_dac_number - 1 ]
            V_g = self.ivvi.dac_voltages()[ V_g_dac_number - 1 ]
            
            I_m = self.kI.amplitude()
            V_m = self.kV.amplitude()

            plt.plot( V_g*V_over_V_g, I_m*( 1/V_over_A_m ), pcolor )
            plt.xlabel('V_g (mV)')
            plt.ylabel('I_m (A)')
            plt.ticklabel_format( scilimits=(-3, 4) )
            fig.canvas.draw()

            time.sleep( 0.1 )
            
            if filepath == '':
                
                continue
            
            data = np.column_stack( ( V_g, V_b, I_m, V_m ) )  
            
            with open( path, 'a+' ) as f:

                np.savetxt( f, data, delimiter = '\t' ) 
                f.flush()  
        
        if sweepback:
            
            self.Sweep_Vg_Vb_Im_Vm( filepath, V_g_dac, V_b_dac, V_b, V_g_sweep[ :: -1 ], V_over_V_g, V_over_A_m, V_over_V_m, fig, False, V_g_ramping_mV_per_second, pcolor = 'ro' )
            
            
    def Sweep_Vg_Vb_Im_Vm_Igm(
                            filepath: str,
                            V_g_dac: str,
                            V_b_dac:str,
                            V_b: float,
                            V_g_sweep: list,
                            V_over_V_g: float,
                            V_over_A_m: float,
                            V_over_V_m: float,
                            V_over_A_m_Ig: float,
                            fig = False,
                            back = True,
                            V_g_ramping_mV_per_second = 50, 
                            pcolor = 'bo',
                            dac_factor = 1e-3
                    ):

        V_b_dac_number = ReadDacNumber( V_b_dac )
        V_g_dac_number = ReadDacNumber( V_g_dac )
        
        print( 'ramping dacs...')
        
        self.ivvi.set( V_b_dac , V_b )
        RampGate( V_g_dac, V_g_sweep[ 0 ], V_over_V_g, V_g_ramping_mV_per_second )
        
        print( 'dacs at target' )
                
        for i in V_g_sweep:

            self.ivvi.set( V_g_dac, i/V_over_V_g )

            V_b = self.ivvi.dac_voltages()[ V_b_dac_number - 1 ]
            V_g = self.ivvi.dac_voltages()[ V_g_dac_number - 1 ]
            
            I_m = self.kI.amplitude()
            V_m = self.kV.amplitude()
            Ig_m = self.kIG.amplitude()
            
            plt.subplot( 1, 2, 1 )
            plt.plot( V_g*V_over_V_g, I_m*( 1/V_over_A_m ), pcolor )
            plt.xlabel('V_g (mV)')
            plt.ylabel('I_m (A)')
            plt.ticklabel_format( scilimits=(-3, 4) )
            fig.canvas.draw()
            
            plt.subplot( 1, 2, 2 )
            plt.plot( V_g*V_over_V_g, Ig_m*( 1/V_over_A_m_Ig ), pcolor )
            plt.xlabel('V_g (mV)')
            plt.ylabel('I_m (A)')
            plt.ticklabel_format( scilimits=(-3, 4) )
            fig.canvas.draw()

            time.sleep( 0.1 )
            
            if filepath == '':
                
                continue
            
            data = np.column_stack( ( V_g, V_b, I_m, V_m, Ig_m ) )  
            
            with open( path, 'a+' ) as f:

                np.savetxt( f, data, delimiter = '\t' ) 
                f.flush()  
        
        if sweepback:
            
            self.Sweep_Vg_Vb_Im_Vm_Igm( filepath, V_g_dac, V_b_dac, V_b, V_g_sweep[ :: -1 ], V_over_V_g, V_over_A_m, V_over_V_m, V_overA_m_Ig, fig, False, V_g_ramping_mV_per_second, pcolor = 'ro' )
            
            
            
            
    def Sweep_Vs_Vg_Im_Vm(
                        filepath: str,
                        V_s_dac: str,
                        V_g_dac:str,
                        V_g: float,
                        V_s_sweep: list,
                        V_over_V_s: float,
                        V_over_V_g: float,
                        V_overo_V_m: float,
                        V_over_A_m: float,
                        fig = False,
                        back = True,
                        V_g_ramping_mV_per_second = 50, 
                        pcolor = 'bo',
                        dac_factor = 1e-3
                    ):

        V_s_dac_number = ReadDacNumber( V_s_dac )
        V_g_dac_number = ReadDacNumber( V_g_dac )
        
        print( 'ramping dacs...')
        
        self.ivvi.set( V_s_dac , V_s_sweep[ 0 ] )
        RampGate( V_g_dac, V_g, V_over_V_g, V_g_ramping_mV_per_second )
        
        print( 'dacs at target' )
                
        for i in V_s_sweep:

            self.ivvi.set( V_s_dac, i )

            V_s = self.ivvi.dac_voltages()[ V_s_dac_number - 1 ]
            V_g = self.ivvi.dac_voltages()[ V_g_dac_number - 1 ]
            
            I_m = self.kI.amplitude()
            V_m = self.kV.amplitude()

            plt.plot( V_s*V_over_V_s, I_m*( 1/V_over_A_m ), pcolor )
            plt.xlabel('V_g (mV)')
            plt.ylabel('I_m (A)')
            plt.ticklabel_format( scilimits=(-3, 4) )
            fig.canvas.draw()

            time.sleep( 0.1 )
            
            if filepath == '':
                
                continue
            
            data = np.column_stack( ( V_s, V_g, V_m, I_m ) )  
            
            with open( path, 'a+' ) as f:

                np.savetxt( f, data, delimiter = '\t' ) 
                f.flush()  
        
        if sweepback:
            
            self.Sweep_Vs_Vg_Im_Vm( filepath, V_s_dac, V_g_dac, V_g, V_s_sweep[ :: -1 ], V_over_V_s, V_over_V_g, V_over_A_m, V_over_V_m, fig, False, V_g_ramping_mV_per_second, pcolor = 'ro' )
            
            
    def Sweep_Vs_Vg_Im_Vm_Igm(
                            filepath: str,
                            V_s_dac: str,
                            V_g_dac:str,
                            V_g: float,
                            V_s_sweep: list,
                            V_over_V_s: float,
                            V_over_V_g: float,
                            V_over_A_m: float,
                            V_over_V_m: float,
                            V_over_A_m_Ig: float,
                            fig = False,
                            back = True,
                            V_g_ramping_mV_per_second = 50, 
                            pcolor = 'bo',
                            dac_factor = 1e-3
                    ):

        V_s_dac_number = ReadDacNumber( V_s_dac )
        V_g_dac_number = ReadDacNumber( V_g_dac )
        
        print( 'ramping dacs...')
        
        self.ivvi.set( V_s_dac , V_s_sweep[ 0 ] )
        RampGate( V_g_dac, V_g, V_over_V_g, V_g_ramping_mV_per_second )
        
        print( 'dacs at target' )
                
        for i in V_g_sweep:

            self.ivvi.set( V_g_dac, i/V_over_V_g )

            V_s = self.ivvi.dac_voltages()[ V_s_dac_number - 1 ]
            V_g = self.ivvi.dac_voltages()[ V_g_dac_number - 1 ]
            
            I_m = self.kI.amplitude()
            V_m = self.kV.amplitude()
            Ig_m = self.kIG.amplitude()
            
            plt.subplot( 1, 2, 1 )
            plt.plot( V_s*V_over_V_s, I_m*( 1/V_over_A_m ), pcolor )
            plt.xlabel('V_g (mV)')
            plt.ylabel('I_m (A)')
            plt.ticklabel_format( scilimits=(-3, 4) )
            fig.canvas.draw()
            
            plt.subplot( 1, 2, 2 )
            plt.plot( V_s*V_over_V_s, Ig_m*( 1/V_over_A_m_Ig ), pcolor )
            plt.xlabel('V_g (mV)')
            plt.ylabel('I_m (A)')
            plt.ticklabel_format( scilimits=(-3, 4) )
            fig.canvas.draw()

            time.sleep( 0.1 )
            
            if filepath == '':
                
                continue
            
            data = np.column_stack( ( V_s, V_g, I_m, V_m, Ig_m ) )  
            
            with open( path, 'a+' ) as f:

                np.savetxt( f, data, delimiter = '\t' ) 
                f.flush()  
        
        if sweepback:
            
            self.Sweep_Vs_Vg_Im_Vm_Igm( filepath, V_g_dac, V_b_dac, V_b, V_g_sweep[ :: -1 ], V_over_V_s, V_over_V_g, V_over_A_m, V_over_V_m, V_overA_m_Ig, fig, False, V_g_ramping_mV_per_second, pcolor = 'ro' )
            


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def Measurement_Vs_Im( 
                            self,
                            filepath: str,
                            dac: str,
                            V_over_V_s: float,
                            V_over_A_m: float,
                            sweep: list,
                            sweepback: bool,
                            fig = False,
                            pcolor = 'bo',
                            dac_factor = 1e-3
                            ):
                            
        if filepath != '':

            start = datetime.datetime.now()
            timestart = '#\tStart: {}\n'.format( start )
            InsertText( timestart, filepath )

        
        self.Measurement_Vs_Im( filepath, dac, V_over_V_s, V_over_A_m, sweep, False, fig )
        
        
        if filepath != '':

            stopt = datetime.datetime.now()
            timestop = '#\tEnd: {}\n'.format( stopt )
            InsertText( timestop, filepath )
            
            
    def Measurement_Vg_Vb_Im_Vm( 
                            filepath: str,
                            V_g_dac: str,
                            V_b_dac:str,
                            V_b: float,
                            V_g_sweep: list,
                            V_over_V_g: float,
                            V_over_A_m: float,
                            V_over_V_m: float,
                            fig = False,
                            back = True,
                            V_g_ramping_mV_per_second = 50, 
                            pcolor = 'bo',
                            dac_factor = 1e-3
                        ):
                            
        if filepath != '':

            start = datetime.datetime.now()
            timestart = '#\tStart: {}\n'.format( start )
            InsertText( timestart, filepath )

        
        self.Sweep_Vg_Vb_Im_Vm( filepath, V_g_dac, V_b_dac, V_b, V_g_sweep, V_over_V_g, V_over_A_m, V_over_V_m, fig, back, V_g_ramping_mV_per_second )
        
        
        if filepath != '':

            stopt = datetime.datetime.now()
            timestop = '#\tEnd: {}\n'.format( stopt )
            InsertText( timestop, filepath )


    def Measurement_Vg_Vb_Im_Vm_Igm( 
                            filepath: str,
                            V_g_dac: str,
                            V_b_dac:str,
                            V_b: float,
                            V_g_sweep: list,
                            V_over_V_g: float,
                            V_over_A_m: float,
                            V_overV_m: float,
                            V_over_A_m_Ig: float,
                            fig = False,
                            back = True,
                            V_g_ramping_mV_per_second = 50, 
                            pcolor = 'bo',
                            dac_factor = 1e-3
                        ):
                            
        if filepath != '':

            start = datetime.datetime.now()
            timestart = '#\tStart: {}\n'.format( start )
            InsertText( timestart, filepath )

        
        self.Sweep_Vg_Vb_Im_Igm( filepath, V_g_dac, V_b_dac, V_b, V_g_sweep, V_over_V_g, V_over_A_m, V_overV_m, V_over_A_m_Ig, fig, back, V_g_ramping_mV_per_second )
        
        
        if filepath != '':

            stopt = datetime.datetime.now()
            timestop = '#\tEnd: {}\n'.format( stopt )
            InsertText( timestop, filepath )
            
            
            
            
    def Measurement_Vs_Vg_Im_Vm( 
                                filepath: str,
                                V_s_dac: str,
                                V_g_dac:str,
                                V_g: float,
                                V_s_sweep: list,
                                V_over_V_s: float,
                                V_over_V_g: float,
                                V_overo_V_m: float,
                                V_over_A_m: float,
                                fig = False,
                                back = True,
                                V_g_ramping_mV_per_second = 50, 
                                pcolor = 'bo',
                                dac_factor = 1e-3
                            ):
                            
        if filepath != '':

            start = datetime.datetime.now()
            timestart = '#\tStart: {}\n'.format( start )
            InsertText( timestart, filepath )

        
        self.Sweep_Vs_Vg_Im_Igm( filepath, V_s_dac, V_g_dac, V_g, V_s_sweep, V_over_V_s, V_overoV_g, V_over_A_m, V_overV_m, fig, back, V_g_ramping_mV_per_second )
        
        
        if filepath != '':

            stopt = datetime.datetime.now()
            timestop = '#\tEnd: {}\n'.format( stopt )
            InsertText( timestop, filepath )


    def Measurement_Vs_Vg_Im_Vm_Igm( 
                                    filepath: str,
                                    V_s_dac: str,
                                    V_g_dac:str,
                                    V_g: float,
                                    V_s_sweep: list,
                                    V_over_V_s: float,
                                    V_over_V_g: float,
                                    V_over_A_m: float,
                                    V_over_V_m: float,
                                    V_over_A_m_Ig: float,
                                    fig = False,
                                    back = True,
                                    V_g_ramping_mV_per_second = 50, 
                                    pcolor = 'bo',
                                    dac_factor = 1e-3
                            ):
                            
        if filepath != '':

            start = datetime.datetime.now()
            timestart = '#\tStart: {}\n'.format( start )
            InsertText( timestart, filepath )

        
        self.Sweep_Vs_Vg_Im_Igm( filepath, V_s_dac, V_g_dac, V_g, V_s_sweep, V_over_V_s, V_overoV_g, V_over_A_m, V_overV_m, V_over_A_m_Ig, fig, back, V_g_ramping_mV_per_second )
        
        
        if filepath != '':

            stopt = datetime.datetime.now()
            timestop = '#\tEnd: {}\n'.format( stopt )
            InsertText( timestop, filepath )