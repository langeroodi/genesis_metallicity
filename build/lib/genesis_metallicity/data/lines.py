import numpy as np

n = 1.000277

lines_dict = {}

lines_dict['N1238']  = {}
lines_dict['N1238']['lambda']  = 1238.821
lines_dict['N1238']['depend']  = None

lines_dict['N1242']  = {}
lines_dict['N1242']['lambda']  = 1242.804
lines_dict['N1242']['depend']  = 'N1238'

lines_dict['C1334']  = {}
lines_dict['C1334']['lambda']  = 1334.532
lines_dict['C1334']['depend']  = None

lines_dict['C1335']  = {}
lines_dict['C1335']['lambda']  = 1335.708
lines_dict['C1335']['depend']  = 'C1334'

lines_dict['O1397']  = {}
lines_dict['O1397']['lambda']  = 1397.232
lines_dict['O1397']['depend']  = None

lines_dict['O1399']  = {}
lines_dict['O1399']['lambda']  = 1399.780
lines_dict['O1399']['depend']  = 'O1397'

lines_dict['Si1402']  = {}
lines_dict['Si1402']['lambda'] = 1402.770
lines_dict['Si1402']['depend'] = None

lines_dict['C1548']  = {}
lines_dict['C1548']['lambda']  = 1548.187
lines_dict['C1548']['depend']  = None

lines_dict['C1550']  = {}
lines_dict['C1550']['lambda']  = 1550.772
lines_dict['C1550']['depend']  = 'C1548'

lines_dict['He1640']  = {}
lines_dict['He1640']['lambda'] = 1640.420
lines_dict['He1640']['depend'] = None

lines_dict['O1660']  = {}
lines_dict['O1660']['lambda']  = 1660.809
lines_dict['O1660']['depend']  = None

lines_dict['O1666']  = {}
lines_dict['O1666']['lambda']  = 1666.150
lines_dict['O1666']['depend']  = 'O1660'

lines_dict['C1908']  = {}
lines_dict['C1908']['lambda']  = 1908.734
lines_dict['C1908']['depend']  = None

lines_dict['O3727']  = {}
lines_dict['O3727']['lambda']  = 3726.032 * n
lines_dict['O3727']['depend']  = None

lines_dict['O3729']  = {}
lines_dict['O3729']['lambda']  = 3728.815 * n
lines_dict['O3729']['depend']  = 'O3727'

lines_dict['OII'] = {}
lines_dict['OII']['lambda']  = (lines_dict['O3729']['lambda']+lines_dict['O3727']['lambda'])/2
lines_dict['OII']['depend']  = None

lines_dict['Hdelta']  = {}
lines_dict['Hdelta']['lambda'] = 4101.742 * n
lines_dict['Hdelta']['depend'] = None

lines_dict['Hgamma']  = {}
lines_dict['Hgamma']['lambda'] = 4340.471 * n
lines_dict['Hgamma']['depend'] = None

lines_dict['O4363']  = {}
lines_dict['O4363']['lambda']  = 4363.210 * n
lines_dict['O4363']['depend']  = None

lines_dict['Hbeta']  = {}
lines_dict['Hbeta']['lambda']  = 4861.333 * n
lines_dict['Hbeta']['depend']  = None

lines_dict['Hbeta_EW'] = lines_dict['Hbeta']

lines_dict['O4959']  = {}
lines_dict['O4959']['lambda']  = 4958.911 * n
lines_dict['O4959']['depend']  = None

lines_dict['O5007']  = {}
lines_dict['O5007']['lambda']  = 5006.843 * n
lines_dict['O5007']['depend']  = None

lines_dict['OIII']  = {}
lines_dict['OIII']['lambda']  = (4958.911+5006.843)/2 * n
lines_dict['OIII']['depend']  = None

lines_dict['Halpha'] = {}
lines_dict['Halpha']['lambda'] = 6562.819 * n
lines_dict['Halpha']['depend'] = None

lines_dict['O7320']  = {}
lines_dict['O7320']['lambda']  = 7319.990 * n
lines_dict['O7320']['depend']  = None

lines_dict['O7330']  = {}
lines_dict['O7330']['lambda']  = 7330.730 * n
lines_dict['O7330']['depend']  = None

lines_dict['OII7320'] = {}
lines_dict['OII7320']['lambda']  = (lines_dict['O7320']['lambda']+lines_dict['O7330']['lambda'])/2
lines_dict['OII7320']['depend']  = None

#---- function to list the required lines and their metallicity ----#

required_lines = np.array(['OII', 'Hbeta', 'Hbeta_EW', 'O4959', 'O5007'])
optional_lines = np.array(['Hdelta', 'Hgamma', 'O4363', 'Halpha', 'OII7320'])
backend_lines  = np.concatenate((required_lines, optional_lines))

description_dict = {}
description_dict['OII']     = ': sum of the \'O3727\' and \'O3729\' lines (can accept the components separately)'
description_dict['O4959']   = ': can accept as \'OIII\' if not resolved from \'O5007\''
description_dict['O5007']   = ': can accept as \'OIII\' if not resolved from \'O4959\''
description_dict['OII7320'] = ': sum of the \'O7320\' and \'O7330\' lines (can accept the components separately)'

def print_lines():

    print('------------------------------------------------------------------------------------------')
    print(' -> these are the required lines:')
    print('(missing any of these lines will raise errors)')
    print('..............................................')
    for line in required_lines:
        try:
            description = description_dict[line]
        except:
            description = ' '
        print('    %s %s' %(line, description))

    print('------------------------------------------------------------------------------------------')
    print(' -> these are the optional lines:')
    print('(code can function without them)')
    print('..............................................')
    for line in optional_lines:
        try:
            description = description_dict[line]
        except:
            description = ' '
        print('    %s %s' %(line, description))
