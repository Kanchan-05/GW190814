import numpy as np
from pycbc.types import TimeSeries, FrequencySeries
from pyseobnr.generate_waveform import GenerateWaveform
from pycbc.waveform import get_waveform_filter_length_in_time, get_waveform_end_frequency
from pycbc.waveform import utils as wutils
from pycbc.conversions import get_final_from_initial, tau_from_final_mass_spin
from pycbc import pnutils


def teobresums_length_in_time(**kwargs):
    from pycbc.waveform import waveform
    """ Estimates the duration of SEOBNRv4HM waveforms that include higher modes.
    """
    # Default maximum node number for SEOBNRv4HM is 5
    # The relevant lower order approximant here is SEOBNRv4
    # NOTE: We can also use IMRPhenomD instead of SEOBNRv4 with e('IMRPhenomD', 2, **kwargs)
    from pycbc.waveform import waveform
    return waveform.get_hm_length_in_time('SEOBNRv4', 2, **kwargs)


def teobresums_get_end_frequency(**p):
    """Return the ringdown frequency of a template 
    """
    #f_end = get_waveform_end_frequency("IMRPhenomD", 2,**p)
    # Here we use SEOBNRv4PHM approximant
    import pycbc.conversions as conversions
    final_mass =  conversions.final_mass_from_initial(p['mass1'],p['mass2'],p['spin1z'],p['spin2z'], approximant='SEOBNRv4PHM', f_ref=-1)
    final_spin = conversions.final_spin_from_initial(p['mass1'],p['mass2'],p['spin1z'],p['spin2z'], approximant='SEOBNRv4PHM', f_ref=-1)
    ringdown_freq = conversions.freq_from_final_mass_spin(final_mass, final_spin, l=2, m=2, n=0)
    return ringdown_freq * 1.5