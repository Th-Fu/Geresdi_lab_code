# ------------------------------------------------------------------------------
# Description    : S4g SPI module QCoDeS interface
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

from spirack.S4g_module import S4g_module as S4g_api
from spi_rack.spi_module_base import spi_module_base


class dummy_s4g_api:
    """
    Mock implementation of spirack API (https://github.com/mtiggelman/SPI-rack/blob/master/spirack/S4g_module.py),
    for use with the dummy drivers.
    """

    SPAN_MAPPING = {"range_max_uni": 0, "range_max_bi": 2, "range_min_bi": 4}
    DAC_RESOLUTION = 2 ** 18
    NUMBER_OF_DACS = 4

    def __init__(self, spi_rack, module, max_current=50e-3, reset_currents=True):
        """
        Instantiates the mock communication layer with identical parameters to
        the `spirack.S4g_module.S4g_module` constructor.

        Parameters
        ----------
        spi_rack : dummy_spi_api
            Mock SPI_rack class object via which the communication runs
        module : int
            module number set on the hardware
        max_current : float
            maximum range of the S4g, configured in hardware
        reset_currents : bool
            if True, then reset all currents to zero and change the span to `range_max_bi`

        Returns
        ----------
        """
        self.parent = spi_rack
        self._currents = [0.] * self.NUMBER_OF_DACS
        init_span = self.SPAN_MAPPING["range_max_bi"] if reset_currents else 0
        self._num_dacs = self.NUMBER_OF_DACS
        self._spans = [init_span] * self.NUMBER_OF_DACS
        self.max_current = max_current
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

    def change_span(self, DAC, span):
        """
        Mocks the `change_span` function of the API

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

    def set_current(self, DAC: int, current: float):
        """
        Mocks the `set_voltage` function of the API

        Parameters
        ----------
        DAC: int
            Current output of which to update the current
        voltage: float
            new DAC current
        """
        self._currents[DAC] = current

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
            current
        int
            span
        """
        return self._currents[DAC], self._spans[DAC]

    def get_stepsize(self, DAC):
        """
        Mocks the `get_stepsize` function of the API

        Parameters
        ----------
        DAC : int
            Current output of which the stepsize is calculated

        Returns
        ----------
        float
            Smallest current step possible with DAC
        """
        if self._spans[DAC] == self.SPAN_MAPPING["range_max_bi"]:
            return 2 * self.max_current / self.DAC_RESOLUTION

        return self.max_current / self.DAC_RESOLUTION


class s4g_module(spi_module_base):
    """
    `QCoDeS <https://qcodes.github.io/Qcodes/>`- style instrument channel driver for the S4g SPI module.
    """

    NUMBER_OF_DACS = 4  # Set by hardware constraints

    def __init__(
        self,
        parent,
        name: str,
        address: int,
        reset_currents: bool = True,
        dac_names: Optional[List[str]] = None,
        is_dummy: bool = False,
    ):
        """
        Instantiates the driver object. This is the object that should be instantiated by the add_spi_module function.

        Parameters
        ----------
        parent
            Reference to the spi_rack parent object. This is handled by the add_spi_module function.
        name : str
            Name given to the InstrumentChannel.
        address : int
            Module number set on the hardware.
        reset_currents : bool
            If True, then reset all currents to zero and change the span to `range_max_bi`.
        dac_names : List[str]
            List of all the names to use for the dac channels. If no list is given or is None, the default name "dac{i}"
            is used for the i-th dac channel.

        Returns
        ----------

        Raises
        ----------
        ValueError
            Length of the dac names list does not match the number of dacs.
        """
        super().__init__(parent, name, address)

        api = dummy_s4g_api if is_dummy else S4g_api

        self.api = api(parent.spi_rack, module=address, reset_currents=reset_currents)
        self._channels = []

        for dac in range(self.NUMBER_OF_DACS):
            if dac_names == None:
                ch_name = "dac{}".format(dac)
            elif len(dac_names) == self.NUMBER_OF_DACS:
                ch_name = dac_names[dac]
            else:
                raise ValueError(f"Length of dac_names must be {self.NUMBER_OF_DACS}")
            channel = s4g_dac_channel(self, ch_name, dac)
            self._channels.append(channel)
            self.add_submodule(ch_name, channel)

    def set_dacs_zero(self):
        """
        Sets all currents of all outputs to 0.

        Returns
        ----------
        """
        for ch in self._channels:
            ch.current(0)


class s4g_dac_channel(InstrumentChannel):
    """
    `QCoDeS <https://qcodes.github.io/Qcodes/>`- style instrument channel driver for the dac channels of the S4g
    module. This class is used by the S4g_module to define the individual dac channels and should not be used
    directly.
    """

    def __init__(self, parent: s4g_module, name: str, dac: int):
        """
        Constructor for the dac channel instrument channel.

        Parameters
        ----------
        parent : s4g_module
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
            "current",
            get_cmd=partial(self._get_current, dac),
            set_cmd=partial(self._set_current, dac),
            unit="A",
            vals=validators.Numbers(min_value=-40e-3, max_value=40e-3),
            docstring="Sets the output current of the dac channel. Depending on the value of "
            "ramping_enabled, the output value is either achieved through slowly ramping, or "
            "instantaneously set.",
        )
        self.add_parameter(
            "span",
            val_mapping={"range_max_uni": 0, "range_max_bi": 2, "range_min_bi": 4},
            get_cmd=partial(self._get_span, dac),
            set_cmd=partial(self._api.change_span_update, dac),
            docstring="Sets the max range of the DACs. Possible values:"
            "\t'range_max_uni':\t0 - 40 mA,"
            "\t'range_max_bi':\t-40 - 40 mA,"
            "\t'range_min_bi':\t-20 - 20 mA.",
        )

        self.add_parameter(
            "ramp_rate",
            unit="A/s",
            initial_value=1e-3,  # 1 mA/s
            docstring="Limits the rate at which currents can be changed. The size of of steps is still "
            "limited by `ramp_max_step`.",
            parameter_class=ManualParameter,
        )
        self.add_parameter(
            "ramp_max_step",
            unit="A",
            initial_value=0.5e-3,
            docstring="Sets the maximum step size for current ramping. The rate at which it ramps is set"
            " by `ramp_rate`.",
            parameter_class=ManualParameter,
        )
        self.add_parameter(
            "ramping_enabled",
            initial_value=False,
            vals=validators.Bool(),
            parameter_class=ManualParameter,
            docstring="Turns ramping on or off. Toggling `ramping_enabled` changed the behavior of the "
            "setter for the `current` parameter. If enabled, ramping is done at a rate set by "
            "`ramp_rate` and in steps specified by `ramp_max_step`.",
        )
        self.add_parameter(
            "is_ramping",
            get_cmd=lambda: self._is_ramping,
            set_cmd=False,
            docstring="Returns whether the dac is currently in the process of ramping.",
        )
        self.add_parameter(
            "stepsize",
            unit="A",
            set_cmd=False,
            get_cmd=partial(self._api.get_stepsize, dac),
            docstring="Returns the smallest current step allowed by the dac for the current settings.",
        )
        self.add_parameter('dac_channel',
                           set_cmd=False,
                           get_cmd=lambda: dac,
                           docstring="Returns the dac number of this channel."
                           )

    def _get_span(self, dac: int):
        """
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

    def _get_current(self, dac: int):
        """
        Gets the current set by the module.

        Parameters
        ----------
        dac : int
            the dac of which to get the current

        Returns
        ----------
        float
            The output current reported by the hardware
        """
        current, _ = self._api.get_settings(dac)
        return current

    def _set_current(self, dac: int, val: float):
        """
        Sets the current either through ramping or instantly.

        Parameters
        ----------
        dac : int
            the dac of which to set the current
        val : float
            The new value of the current

        Returns
        ----------
        """
        if self.ramping_enabled():
            self._set_current_ramp(dac, val)
        else:
            self.set_current_instant(dac, val)

    def _set_current_ramp(self, dac, val):
        """
        Ramps the current in steps set by `ramp_max_step` with a rate set by `ramp_rate`. Ramping is non-blocking so the
        user should check `is_ramping() == False` to see if the final value is reached.

        Parameters
        ----------
        dac : int
            the dac of which to set the current
        val : float
            The new value of the current after ramping

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
                    ch.set_current_instant(dac, val)
                ch._is_ramping = False

            def exit(self):
                """"Stops the thread."""
                self._stopped = True

        # Start of function code

        if self._is_ramping:
            self._ramp_thread.exit()  # stop ramping to the value set previously, and ramp to new value.

        self._is_ramping = True

        val_begin = self.current()
        num_steps = int(np.ceil(np.abs(val_begin - val) / self.ramp_max_step()))
        vals = np.linspace(val_begin, val, num_steps)
        dt = np.abs(val_begin - val) / self.ramp_rate() / num_steps
        th = RampThread(args=(self, dac, vals, dt))
        self._ramp_thread = th
        th.start()

    def set_current_instant(self, dac, val):
        """ "
        Wrapper function around the set_current API call. Instantaneously sets the current.

        Parameters
        ----------
        dac : int
            the dac of which to set the current
        val : float
            The new value of the current

        Returns
        ----------
        """
        self._api.set_current(dac, val)
