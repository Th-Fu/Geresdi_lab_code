### I. Cools - R583423141-1_6way mw switch driver using SPI rack DACs
### 12/04/2022

import numpy as np
import threading
from time import sleep

from qcodes.instrument.channel import InstrumentChannel
from qcodes.instrument.parameter import ManualParameter
from qcodes import validators

from functools import partial
from typing import List, Optional

from spirack.D5a_module import D5a_module as D5a_api
from spi_rack.spi_module_base import spi_module_base
from spi_rack.spi_rack import d5a_module


####################################################        
            #### Added some stuff ####
####################################################

class mw_switch_class(d5a_module):
    def __init__(
        self,
        name: str,
        dac_list: list,
        turn_ramp_off: bool = True,
        reset_switch: bool = True
    ):
        
        self.dac_list = dac_list
        if len(set(dac_list)) != 6:
             raise ValueError(f"Length of dac_list must be 6, not {len(set(dac_list))}. Names have to be unique.")
        
        if turn_ramp_off:
            self.disable_ramp_for_all()
            
        #### Changing span   
        print("Hopefully everything is grounded and off. You've still got 10s..\nSetting the range on the DACs to 4V_uni and 0V output.\nVoltage spikes MAY and probably WILL occur!")
        sleep(10)
        for dac in self.dac_list:
            # This sets uni span "range_4V_uni"
            dac.span('range_4V_uni')
            dac.voltage(0)
            # corresponding: SPAN_MAPPING = 
            # {"range_4V_uni": 0, "range_8V_uni": 1, "range_4V_bi": 2, "range_8V_bi": 3, "range_2V_bi": 4}

        if reset_switch:
            print('By the way, resetting the switch\n Wait atleast 6s!')
            self.reset_all_ports()
            
    def open_port(self, port: int):
        '''
        Opens a single port. Uses dac list for that.
        '''
        # for numbering of ports, see Radiall R583423141 datasheet
        if port > 6 or port < 1:
            raise TypeError(f"Port should be inclusive between 1 and 6")
        dac = self.dac_list[port-1]    
        print(f"Opening port {port}, voltage provided by: {dac}")
        dac.voltage(4)

        
    def close_port(self, port: int):
        '''
        Closes a single port. Uses dac list for that.
        '''
        # for numbering of ports, see Radiall R583423141 datasheet
        # just shuts of the voltage really
        if port > 6 or port < 1:
            raise TypeError(f"Port should be inclusive between 1 and 6")   
        dac = self.dac_list[port-1]
        print(f"Closing port {port}, shutting down: {dac} to 0V")
        dac.voltage(0)
        
    def close_all_ports(self):
        '''
        Basically closes all ports. Nomen est omen.
        '''
        for port in range(1,7):
            self.close_port(port)
            
    def reset_all_ports(self):
        '''
        Resets all ports to closed. First it opens all of them though, in case one port was left open somehow.
        Good to do at the start of your measurements!
        '''
        for port in range(1,7):
            self.open_port(port)
            sleep(0.5)
            self.close_port(port)
            sleep(0.5)
            
    def enable_ramp_for_all(self):
        '''I don't know why you'd need this function but it's there for ya!'''
        for dac in self.dac_list:
            dac.ramping_enabled(False)
            
    def disable_ramp_for_all(self):
        '''I don't know why you'd need this function but it's there for ya!'''
        for dac in self.dac_list:
            dac.ramping_enabled(False)
            
    def open_port_list(self):
        '''This function is fallible: essentially it just measures the applied voltages on each dac and reacts accordingly.'''
        counter = 0
        open_ports = []
        
        for dacco in self.dac_list:
            check_val = eval(f"dacco.get_settings(counter)[0]")
            
            if float(check_val) != float(0.0):
                open_ports.append(f"port {counter+1}")
            if float(check_val) < 0.0:
                print("Excuse me? Please only apply positive voltages...")
                raise ValueError 
            
            counter +=1
        if len(open_ports) == 0:
            print("No open ports! (i.e. no applied voltages)")
        elif len(open_ports) == 1:
            print(f"The only open port is {open_ports[0]}.")
        else:
            print(f"The open ports are {', '.join(open_ports)}.")
            
        
            
####################################################        
            #### ----- end ----- ####
####################################################    