def teobresums_length_in_time(**kwargs):
    from pycbc.waveform import waveform
    """ Estimates the duration of SEOBNRv4HM waveforms that include higher modes.
    """
    # Default maximum node number for SEOBNRv4HM is 5
    # The relevant lower order approximant here is SEOBNRv4
    # NOTE: We can also use IMRPhenomD instead of SEOBNRv4 with e('IMRPhenomD', 2, **kwargs)
    return waveform.get_hm_length_in_time('IMRPhenomD', 2, **kwargs)

def teobresums_get_end_frequency(**p):
    """Estimate the ringdown frequency of a template 
    """
    # Here we use SEOBNRv4PHM approximant
    import pycbc.conversions as conversions
    final_mass =  conversions.final_mass_from_initial(p['mass1'],p['mass2'],p['spin1z'],p['spin2z'], approximant='SEOBNRv4PHM', f_ref=-1)
    final_spin = conversions.final_spin_from_initial(p['mass1'],p['mass2'],p['spin1z'],p['spin2z'], approximant='SEOBNRv4PHM', f_ref=-1)
    ringdown_freq = conversions.freq_from_final_mass_spin(final_mass, final_spin, l=2, m=2, n=0)
    # 1.5 is a fudge factor
    return ringdown_freq * 1.5


