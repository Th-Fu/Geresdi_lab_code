# ------------------------------------------------------------------------------
# Description    : D5a SPI module QCoDeS interface
# Git repository : https://gitlab.com/qblox/packages/software/qblox_instruments.git
# Copyright (C) Qblox BV (2021)
# ------------------------------------------------------------------------------
import numpy as np
import threading
from time import sleep

from qcodes.instrument.channel import InstrumentChannel
from qcodes.instrument.parameter import ManualParameter
from qcodes import validators

from functools import partial
from typing import List, Optional

import Geresdi_lab_code.AC
from Geresdi_lab_code.AC.spirack.D5a_module import D5a_module as D5a_api
from Geresdi_lab_code.AC.spi_rack.spi_module_base import spi_module_base


class dummy_d5a_api:
    """
    Mock implementation of spirack API (https://github.com/mtiggelman/SPI-rack/blob/master/spirack/D5a_module.py),
    for use with the dummy drivers.
    """

    SPAN_MAPPING = {"range_4V_uni": 0, "range_8V_uni": 1, "range_4V_bi": 2, "range_8V_bi": 3, "range_2V_bi": 4}
    DAC_RESOLUTION = 2**18

    def __init__(self, spi_rack, module: int, reset_voltages: bool =True, num_dacs: int = 16):
        """
        Instantiates the mock communication layer with identical parameters to the `spirack.D5a_module.D5a_module` constructor.

        Parameters
        ----------
        spi_rack : dummy_spi_api
            Mock SPI_rack class object via which the communication runs
        module : int
            module number set on the hardware
        reset_voltages : bool
            if True, then reset all voltages to zero and change the span
            to `range_4V_bi`. If a voltage jump would occur, then ramp
            to zero in steps of 10 mV
        num_dacs: int
            number of DAC channels available

        Returns
        ----------
        """
        self.parent = spi_rack
        self._voltages = [0.]*num_dacs
        init_span = self.SPAN_MAPPING["range_4V_bi"] if reset_voltages else 0
        self._num_dacs = num_dacs
        self._spans = [init_span] * num_dacs
        self.address = module

    def change_span_update(self, DAC: int, span: int):
        """
        Mocks the `change_span_update` function of the API.

        Parameters
        ----------
        DAC : int
            Current output of which to change the span
        span : int
            values for the span as mentioned in the datasheet, use from SPAN_MAPPING

        Returns
        ----------
        """
        self._spans[DAC] = span

    def change_span(self, DAC: int, span: int):
        """
        Mocks the `change_span` function of the API.

        Parameters
        ----------
        DAC : int
            Current output of which to change the span
        span : int
            values for the span as mentioned in the datasheet, use from SPAN_MAPPING

        Returns
        ----------
        """
        self.change_span_update(DAC, span)

    def set_voltage(self, DAC: int, voltage: float):
        """
        Mocks the `set_voltage` function of the API

        Parameters
        ----------
        DAC: int
            DAC inside the module of which to update the voltage
        voltage: float
            new DAC voltage

        Returns
        ----------
        """
        self._voltages[DAC] = voltage

    def get_settings(self, DAC: int):
        """
        Mocks the `get_settings` function of the API

        Parameters
        ----------
        DAC : int
            Current output of which the settings will be read

        Returns
        ----------
        float
            voltage
        int
            span
        """
        return self._voltages[DAC], self._spans[DAC]

    def get_stepsize(self, DAC: int):
        """
        Mocks the `get_stepsize` function of the API.

        Parameters
        ----------
        DAC : int
            DAC inside the module of which the stepsize is calculated

        Returns
        ----------
        float
            Smallest voltage step possible with DAC
        """
        if DAC not in range(self._num_dacs):
            raise ValueError('D5a module {} [get_stepsize]: DAC {} does not exist.'.format(self.address, DAC))

        voltage_ranges = {"range_4V_uni": 4.0, "range_8V_uni": 8.0, "range_4V_bi": 8.0, "range_8V_bi": 16.0, "range_2V_bi": 4}
        span_as_str = next(key for key, value in self.SPAN_MAPPING.items() if value == self._spans[DAC])  # reverse dict lookup
        return voltage_ranges[span_as_str] / self.DAC_RESOLUTION


class d5a_module(spi_module_base):
    """
    `QCoDeS <https://qcodes.github.io/Qcodes/>`- style instrument channel driver for the D5a SPI module.
    """

    NUMBER_OF_DACS = 16  # Set by hardware constraints

    def __init__(
        self,
        parent,
        name: str,
        address: int,
        reset_voltages: bool = True,
        dac_names: Optional[List[str]] = None,
        is_dummy: bool = False
    ):
        
        """
        Instantiates the driver object.

        Parameters
        ----------
        parent
            Reference to the spi_rack parent object. This is handled by the add_spi_module function.
        name : str
            Name given to the InstrumentChannel.
        address : int
            Module number set on the hardware.
        reset_voltages : bool
            If True, then reset all voltages to zero and change the span to `range_max_bi`.
        dac_names : List[str]
            List of all the names to use for the dac channels. If no list is given or is None, the default name "dac{i}"
            is used for the i-th dac channel.
        is_dummy : bool
            If true, do not connect to physical hardware, but use

        Returns
        ----------

        Raises
        ----------
        ValueError
            Length of the dac names list does not match the number of dacs.
        """
        super().__init__(parent, name, address)

        api = dummy_d5a_api if is_dummy else D5a_api

        self.api = api(
            parent.spi_rack,
            module=address,
            reset_voltages=reset_voltages,
            num_dacs=self.NUMBER_OF_DACS,
        )
        self._channels = []
        
        for dac in range(self.NUMBER_OF_DACS):
            if dac_names == None:
                ch_name = "dac{}".format(dac)
            elif len(dac_names) == self.NUMBER_OF_DACS:
                ch_name = dac_names[dac]
            else:
                raise ValueError(f"Length of dac_names must be {self.NUMBER_OF_DACS}")
            channel = d5a_dac_channel(self, ch_name, dac)
            self._channels.append(channel)
            self.add_submodule(ch_name, channel)
        
    def set_dacs_zero(self):
        """
        Sets all voltages of all outputs to 0.

        Returns
        ----------
        """
        for ch in self._channels:
            ch.voltage(0)
          
            

class d5a_dac_channel(InstrumentChannel):
    """
    `QCoDeS <https://qcodes.github.io/Qcodes/>`- style instrument channel driver for the dac channels of the D5a
    module. This class is used by the d5a_module to define the individual dac channels and should not be used
    directly.
    """

    def __init__(self, parent: d5a_module, name: str, dac: int):
        """
        Constructor for the dac channel instrument channel.

        Parameters
        ----------
        parent : d5a_module
            Reference to the parent s4g_module
        name : str
            Name for the instrument channel
        dac : int
            Number of the dac that this channel corresponds to

        Returns
        ----------

        Raises
        ----------
        """
        super().__init__(parent, name)
        self._api = parent.api
        self._is_ramping = False
        self._ramp_thread = None

        self.add_parameter(
            "voltage",
            get_cmd=partial(self._get_voltage, dac),
            set_cmd=partial(self._set_voltage, dac),
            unit="V",
            vals=validators.Numbers(min_value=-8.0, max_value=8.0),
            docstring="Sets the output voltage of the dac channel. Depending on the value of "
            "ramping_enabled, the output value is either achieved through slowly ramping, or "
            "instantaneously set.",
        )
        
        self.add_parameter(
            "span",
            val_mapping={
                "range_4V_uni": 0,
                "range_8V_uni": 1,
                "range_4V_bi": 2,
                "range_8V_bi": 3,
                "range_2V_bi": 4,
            },
            get_cmd=partial(self._get_span, dac),
            set_cmd=partial(self._api.change_span_update, dac),
            docstring="Sets the max range of the DACs. Possible values:"
            "\t'range_4V_uni':\t0 - 4 V,"
            "\t'range_8V_uni':\t0 - 8 V (only if non-standard 12 V power supply is present),"
            "\t'range_4V_bi':\t-4 - 4 V,"
            "\t'range_8V_bi':\t-8 - 8 V (only if non-standard 12 V power supply is present),"
            "\t'range_2V_bi':\t-2 - 2 V.",
        )
        self.add_parameter('ramp_rate',
                           unit='V/s',
                           initial_value=100e-3,  # 10 mV/s
                           docstring='Limits the rate at which currents can be changed. The size of of steps is still '
                                     'limited by `ramp_max_step`.',
                           parameter_class=ManualParameter
                           )
        self.add_parameter('ramp_max_step',
                           unit='V',
                           initial_value=100e-3,
                           docstring='Sets the maximum step size for voltage ramping. The rate at which it ramps is set'
                                     ' by `ramp_rate`.',
                           parameter_class=ManualParameter
                           )
        self.add_parameter('ramping_enabled',
                           initial_value=False,
                           vals=validators.Bool(),
                           parameter_class=ManualParameter,
                           docstring='Turns ramping on or off. Toggling `ramping_enabled` changed the behavior of the '
                                     'setter for the `current` parameter. If enabled, ramping is done at a rate set by '
                                     '`ramp_rate` and in steps specified by `ramp_max_step`.'
                           )
        self.add_parameter('is_ramping',
                           get_cmd=lambda: self._is_ramping,
                           set_cmd=False,
                           docstring="Returns whether the dac is currently in the process of ramping.")
        self.add_parameter('stepsize',
                           unit='V',
                           set_cmd=False,
                           get_cmd=partial(self._api.get_stepsize, dac),
                           docstring="Returns the smallest current step allowed by the dac for the current settings; usually 4V/(2^18)"
                           )
        self.add_parameter('dac_channel',
                           set_cmd=False,
                           get_cmd=lambda: dac,
                           docstring="Returns the dac number of this channel."
                           )

    def _get_span(self, dac: int):
        """"
        Gets the span set by the module.

        Parameters
        ----------
        dac : int
            the dac of which to get the span

        Returns
        ----------
        int
            The current span
        """
        _, span = self._api.get_settings(dac)
        return span

    def _get_voltage(self, dac: int):
        """"
        Gets the voltage set by the module.

        Parameters
        ----------
        dac : int
            the dac of which to get the voltage

        Returns
        ----------
        float
            The output voltage reported by the hardware
        """
        voltage, _ = self._api.get_settings(dac)
        return voltage

    def _set_voltage(self, dac: int, val: float):
        """
        Sets the voltage either through ramping or instantly.

        Parameters
        ----------
        dac : int
            the dac of which to set the voltage
        val : float
            The new value of the voltage

        Returns
        ----------
        """
        if self.ramping_enabled():
            self._set_voltage_ramp(dac, val)
        else:
            self.set_voltage_instant(dac, val)

    def set_voltage_instant(self, dac: int, val: float):
        """"
        Wrapper function around the set_voltage API call. Instantaneously sets the voltage.

        Parameters
        ----------
        dac : int
            the dac of which to set the voltage
        val : float
            The new value of the voltage

        Returns
        ----------
        """
        self._api.set_voltage(dac, val)

    def _set_voltage_ramp(self, dac: int, val: float):
        """
        Ramps the voltage in steps set by `ramp_max_step` with a rate set by `ramp_rate`. Ramping is non-blocking so the
        user should check `is_ramping() == False` to see if the final value is reached.

        Parameters
        ----------
        dac : int
            the dac of which to set the voltage
        val : float
            The new value of the voltage after ramping

        Returns
        ----------
        """
        class RampThread(threading.Thread):
            """"Inner class that defines a thread that can be safely killed."""
            def __init__(self, *args, **kwargs):
                super().__init__(target=self.worker_function, *args, **kwargs)
                self._stopped = False

            def worker_function(self, ch, dac, vals, dt):
                """"Conducts the actual ramping"""
                for val in vals:
                    if self._stopped:
                        return
                    sleep(dt)
                    ch.set_voltage_instant(dac, val)
                ch._is_ramping = False

            def exit(self):
                """"Stops the thread."""
                self._stopped = True

        # Start of function code

        if self._is_ramping:
            self._ramp_thread.exit()  # stop ramping to the value set previously, and ramp to new value.

        self._is_ramping = True

        val_begin = self.voltage()
        num_steps = int(np.ceil(np.abs(val_begin - val)/self.ramp_max_step()))
        vals = np.linspace(val_begin, val, num_steps)
        dt = np.abs(val_begin - val)/self.ramp_rate()/num_steps
        th = RampThread(args=(self, dac, vals, dt))
        self._ramp_thread = th
        th.start()

        
        
####################################################        
            #### Added some stuff ####
####################################################
    def get_settings(self, dac: int):
        """Reads current DAC settings

        Reads back the DAC registers of the given DAC for both the code
        and the span. Calculates the voltage set with the read out span.

        Args:
            DAC (int: 0-15): DAC inside the module of which the settings will be read
        Returns:
            List with voltage and span: [voltages, span] (int)
        """
  
        voltage = self._get_voltage(dac) 
        span = self._get_span(dac)
        
        return [voltage, span]

####################################################        
            #### ----- end ----- ####
####################################################    

