"""Translate sampler math labels into the rich text supported by Qt plots."""
from html import escape
import re

from matplotlib._mathtext_data import tex2uni


def plot_label(label):
    """Keep ordinary text safe and format common TeX parameter notation."""
    text = escape(str(label).replace('$', ''))
    text = re.sub(r'\\([A-Za-z]+)',
                  lambda m: chr(tex2uni[m[1]]) if m[1] in tex2uni else m[0], text)
    # Font commands are unnecessary for Qt's default label font.
    text = re.sub(r'\\(?:mathrm|textrm|text|mathit|mathbf|operatorname)\s*', '', text)
    text = text.replace(r'\,', ' ').replace(r'\;', ' ').replace(r'\!', '')
    text = text.replace(r'\left', '').replace(r'\right', '')
    # Unwrap font groups before processing scripts (e.g. _{\\mathrm{RV}}).
    while re.search(r'\{([^{}]*)\}', text):
        unwrapped = re.sub(r'(?<![_^])\{([^{}]*)\}', r'\1', text)
        if unwrapped == text:
            break
        text = unwrapped
    # Resolve inner groups first, including nested subscripts/superscripts.
    while re.search(r'([_^])\{([^{}]*)\}', text):
        text = re.sub(r'([_^])\{([^{}]*)\}',
                      lambda m: '<{0}>{1}</{0}>'.format('sub' if m[1] == '_' else 'sup', m[2]), text)
    text = re.sub(r'([_^])([^\s<{}])',
                  lambda m: '<{0}>{1}</{0}>'.format('sub' if m[1] == '_' else 'sup', m[2]), text)
    return text.replace('{', '').replace('}', '')
