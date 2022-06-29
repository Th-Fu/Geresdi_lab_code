# ------------------------------------------------------------------------------
# Description    : SPI Module Base Class
# Git repository : https://gitlab.com/qblox/packages/software/qblox_instruments.git
# Copyright (C) Qblox BV (2021)
# ------------------------------------------------------------------------------
from qcodes.instrument.channel import InstrumentChannel
from qcodes.instrument.parameter import ManualParameter
from abc import ABCMeta, abstractmethod


class spi_module_base(InstrumentChannel, metaclass=ABCMeta):
    """
    Defines an abstract base class for SPI modules. All module drivers should inherit from this class.

    This class defines no actual functionality but rather serves to provide a common interface shared among all
    modules.

    Parameters
    ----------
    parent
        Reference to the spi_rack parent object. This is handled by the add_spi_module function.
    name : str
        Name given to the InstrumentChannel.
    address : int
        Module number set on the hardware (set internally on the module itself with a jumper).

    Returns
    ----------

    Raises
    ----------
    """
    def __init__(self, parent, name: str, address: int, **kwargs):
        super().__init__(parent, name)
        self._address = address

    @abstractmethod
    def set_dacs_zero(self):
        """ "
        Base method for set_dacs_zero. Should be overridden by subclass.
        """
        raise NotImplementedError(
            f"{type(self)} has no implementation of set_dacs_zero()"
        )


class dummy_spi_module(spi_module_base):
    """
    A dummy implementation of module driver, used for mock implementations.
    """

    def __init__(self, parent, name: str, address: int):
        super().__init__(parent, name, address)

        self.add_parameter(
            "output",
            unit="A",
            parameter_class=ManualParameter,
            docstring="This is a simple manual parameter that can be used during testing. The user is "
            "able to get and set values to/from it but no actual functionality is associated "
            "with it.",
        )

    def set_dacs_zero(self):
        self.output(0)
