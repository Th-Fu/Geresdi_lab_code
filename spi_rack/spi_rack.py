#------------------------------------------------------------------------------
# Description    : SPI Rack QCoDeS interface
# Git repository : https://gitlab.com/qblox/packages/software/qblox_instruments.git
# Copyright (C) Qblox BV (2021)
#------------------------------------------------------------------------------

from time import time
from typing import Type, List, Union, Optional

from qcodes.instrument.base import Instrument

from qblox_instruments import build

from spirack.spi_rack import SPI_rack as spi_api

from spi_rack.spi_module_base import spi_module_base, dummy_spi_module
from spi_rack.s4g_module import s4g_module
from spi_rack.d5a_module import d5a_module


class dummy_spi_api:
    """
    A dummy API that can be used to test the SPI rack driver.
    """
    TEMPERATURE = 25.0
    BATTERY_LVLS = [6., -6.]
    FIRMWARE_VERSION = "v1.0"

    def __init__(self, address, baud_rate, timeout=1.):
        """
        Instantiates the dummy API object.

        Parameters
        ----------
        address
            Mock value for the address on which the spi_rack is connected. This value is assigned to a member variable,
            but is not actually used.
        address
            Mock value for the baud_rate for the serial connection. This value is assigned to a member variable, but is
            not actually used.
        timeout
            Mock value for the timeout for the serial connection. This value is assigned to a member variable, but is
            not actually used.
        """
        self.address = address
        self.baud_rate = baud_rate
        self.timeout = timeout

        self.locked = True

    def get_temperature(self) -> float:
        """
        Return a mock temperature.

        Returns
        ----------
        float
            returns `dummy_spi_api.TEMPERATURE`
        """
        return self.TEMPERATURE

    def get_battery(self) -> List[float]:
        """
        Return a mock battery level list. In the actual API these are two values read from the battery ADCs. For the
        mock API simply constant values are returned. The expected values to be returned can also be gotten through
        `dummy_spi_api.BATTERY_LVLS`.

        Returns
        ----------
        List[float]
            returns `dummy_spi_api.BATTERY_LVLS`
        """
        return self.BATTERY_LVLS

    def get_firmware_version(self) -> str:
        """
        Returns a firmware version string. In the actual API this is returned by the microcontroller.

        Returns
        ----------
        str
            returns `dummy_spi_api.FIRMWARE_VERSION`
        """
        return self.FIRMWARE_VERSION

    def close(self):
        """Not relevant for dummy api but added to remove errors"""
        pass

    def unlock(self):
        """
        Unlocks the communication to the microcontroller. Since for the dummy implementation there is no communication
        to the microcontroller, this simply sets a self.locked to False.
        """
        self.locked = False


class spi_rack(Instrument):
    """
    SPI rack driver class based on `QCoDeS <https://qcodes.github.io/Qcodes/>`. This class relies on the `spirack API
    <https://github.com/mtiggelman/SPI-rack/>`.

    Example usage:

    .. highlight:: python
    .. code-block:: python

        from spi_rack.spi_rack import spi_rack
        from spi_rack.s4g_module import s4g_module

        spi = spi_rack("my_spi_rack", "COM4")     # connects to an SPI rack on COM port 4
        spi.add_spi_module(3, "D5a", "alice")     # adds an D5a module with address 3 named "alice"
        spi.add_spi_module(2, "S4g", "bob")       # adds an S4g module with address 2 named "bob"
        spi.add_spi_module(6, s4g_module)         # adds an S4g module with address 6 with the default name module6

        spi.bob.dac0.current(10e-3)               # sets the current of output 1 of the S4g module named "bob" to 10 mA
        spi.alice.dac6.voltage(-2)                # sets the voltage of output 7 of the D5a module named "alice" to -2 V
    """

    _MODULES_MAP = {"S4g": s4g_module,
                    "D5a": d5a_module,
                    "dummy": dummy_spi_module}  # If drivers are created for different modules they should be added here

    def __init__(self, name: str, address: str, baud_rate: int = 9600, timeout: float = 1,
                 is_dummy: bool = False):
        """
        Instantiates the driver object.

        Parameters
        ----------
        name : str
            Instrument name.
        address : str
            COM port used by SPI rack controller unit (e.g. "COM4")
        baud_rate : int
            Baud rate
        timeout : float
            Data receive timeout in seconds
        is_dummy : bool
            If true, the SPI rack driver is operating in "dummy" mode for testing purposes.

        Returns
        ----------
        """
        t0 = time()
        super().__init__(name)

        api = dummy_spi_api if is_dummy else spi_api
        self.spi_rack = api(address, baud_rate, timeout=timeout)

        self.spi_rack.unlock()

        self._modules = {}

        self._add_qcodes_params()

        self.connect_message(begin_time=t0)

    def _add_qcodes_params(self):
        """
        Function to add the QCoDeS parameters to the instrument
        """
        self.add_parameter('temperature',
                           unit='C',
                           set_cmd=False,
                           get_cmd=self.spi_rack.get_temperature,
                           docstring="Returns the temperature in the C1b module. Reads the temperature from the "
                                     "internal C1b temperature sensor. Accuracy is +- 0.5 degrees in 0-70 degree range."
                           )
        self.add_parameter('battery_voltages',
                           unit='V',
                           set_cmd=False,
                           get_cmd=self.spi_rack.get_battery,
                           docstring='Calculates the battery voltages from the ADC channel values. '
                                     'Returns: [VbatPlus, VbatMin]'
                           )

    def add_spi_module(self, address: int, module_type: Union[spi_module_base, str], name: str = None, **kwargs):
        """
        Add a module to the driver.

        Parameters
        ----------
        address : int
            Address that the module is set to (set internally on the module itself with a jumper).
        module_type : Union[str, spi_module_base]
             Either a str that is defined in _MODULES_MAP, or a reference to a class derived from SPI_Module_Base.
        name: str
            Optional name of the module. If no name is given or is None, a default name of "module{address}" is used.

        Returns
        ----------

        Raises
        ----------
        ValueError
            module_type is not a string or a subclass of SPI_Module_Base
        """
        if name is None:
            name = "module{}".format(address)

        if isinstance(module_type, str):
            module_obj = self._MODULES_MAP[module_type](self, name, address, **kwargs)
        elif issubclass(module_type, spi_module_base):
            module_obj = module_type(self, name, address, **kwargs)
        else:
            raise ValueError(f"{module_type} is not a valid SPI module.")

        self._modules[address] = module_obj
        self.add_submodule(name, module_obj)

    def get_idn(self):
        """
        Generates the IDN dict.

        Parameters
        ----------

        Returns
        ----------
        dict
            The QCoDeS style IDN dictionary. Currently only the firmware version is actually read from hardware.

        Raises
        ----------
        """
        return {"manufacturer": "Qblox",  # FIXME: This should all be read from microcontroller, not hardcoded.
                "model": "SPI Rack",
                "firmware": {"device": self.spi_rack.get_firmware_version(),
                             "driver": build.get_build_info()}}

    def close(self):
        """
        Closes connection to hardware and closes the Instrument.
        """
        self.spi_rack.close()
        super().close()

    def set_dacs_zero(self):
        """
        Calls the :meth:`set_dacs_zero` function on all the modules, which in turn should cause all output values to be
        set to 0.
        """
        for mod in self._modules.values():
            mod.set_dacs_zero()

    def connect_message(self, idn_param: str = 'IDN',
                        begin_time: Optional[float] = None) -> None:
        """
        Print a standard message on initial connection to an instrument. Overridden from superclass to accommodate
        IEEE488.2 for IDN.

        Parameters
        ----------
        idn_param: str
            Name of parameter that returns ID dict. Default ``IDN``.
        begin_time: Optional[float]
            ``time.time()`` when init started. Default is ``self._t0``, set at start of ``Instrument.__init__``.
        """
        # start with an empty dict, just in case an instrument doesn't
        # heed our request to return all 4 fields.
        idn = {'manufacturer': None, 'model': None,
               'serial': None, 'firmware': None}
        idn.update(self.get(idn_param))
        t = time() - (begin_time or self._t0)

        con_msg = ('Connected to: {manufacturer} {model} '
                   '(serial:{serial}, firmware:{firmware}) '
                   'in {t:.2f}s'.format(t=t, **idn))
        print(con_msg)
        self.log.info(f"Connected to instrument: {idn}")
