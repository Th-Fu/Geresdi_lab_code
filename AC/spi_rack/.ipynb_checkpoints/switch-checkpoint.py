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
        print("Hopefully everything is grounded and off. You've still got 8s.. setting the range to 4V")
        sleep(8)
        for dac in self.dac_list:
            # This sets uni span "range_4V_uni"
            dac.span('range_4V_uni')
            
        if reset_switch:
            print('By the way, resetting the switch')
            self.reset_all_ports()
            

    def open_port(self, port: int):
        # for numbering of ports, see Radiall R583423141 datasheet
        if port > 6 or port < 1:
            raise TypeError(f"Port should be inclusive between 1 and 6") 
        
        dac = self.dac_list[port-1].voltage(4)
    
    def close_port(self, port: int):
        # for numbering of ports, see Radiall R583423141 datasheet
        # just shuts of the voltage really
        if port > 6 or port < 1:
            raise TypeError(f"Port should be inclusive between 1 and 6") 
        dac = self.dac_list[port-1].voltage(0)
        
    def close_all_ports(self):
        for port in range(1,7):
            self.close_port(port)
            
    def reset_all_ports(self):
        for port in range(1,7):
            self.open_port(port)
            sleep(0.25)
            self.close_port(port)
            sleep(0.25)
            
    def enable_ramp_for_all(self):
        for dac in self.dac_list:
            dac.ramping_enabled(False)
            
    def disable_ramp_for_all(self):
        for dac in self.dac_list:
            dac.ramping_enabled(False)
####################################################        
            #### ----- end ----- ####
####################################################    