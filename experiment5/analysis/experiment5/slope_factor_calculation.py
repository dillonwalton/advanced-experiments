import math

def slope_factor(m_bg, m_sig):
    return m_bg / m_sig

def slope_factor_percent(m_bg, m_sig):
    return (m_bg - m_sig) / m_bg

def slope_factor_uncertainty(m_bg, m_bg_err, m_sig, m_sig_err):
    bg = ((1 / m_sig) * m_bg_err)**2
    sig = ((-m_bg / (m_sig**2)) * m_sig_err)**2
    return math.sqrt(bg + sig)