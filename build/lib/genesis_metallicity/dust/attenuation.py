import numpy as np

#########################
# Kriek & Conroy (2013) #
#########################

def KC13(lambda_array, Av, delta=None, Eb=None, Rv=4.05, return_AxAv=False):

    #---- converting between delta & Eb if one is not provided (eq 3) ----#

    if Eb is None:
        Eb = 0.85-1.9*delta
        # print('converted')
    if delta is None:
        delta = (0.85-Eb)/1.9
        # print('converted')

    #---- calzetti ----#

    k_lambda = np.zeros(len(lambda_array))

    transition_lambda = 6300 # AA
    transition_index  = np.argmin(np.abs(lambda_array-transition_lambda))

    # 0.12 - 0.63 micron
    k_lambda[:transition_index] = 2.659 * (-2.156
                                        + 1.509/(lambda_array[:transition_index]*1e-4)
                                        - 0.198/np.power(lambda_array[:transition_index]*1e-4,2)
                                        + 0.011/np.power(lambda_array[:transition_index]*1e-4,3)) + Rv

    # 0.63 - 2.20 micron
    k_lambda[transition_index:] = 2.659 * (-1.857
                                        + 1.040/(lambda_array[transition_index:]*1e-4)) + Rv

    #---- druid ----#

    lambda_0 = 2175 # AA
    delta_lambda = 350 # AA
    D_lambda = Eb*np.power(lambda_array*delta_lambda,2) / (np.power(lambda_array**2-lambda_0**2,2)+np.power(lambda_array*delta_lambda,2))

    #---- kriek & conroy ----#

    # calculating A_lambda
    lambda_v = 5500 # AA
    A_lambda = (Av/Rv) * (k_lambda + D_lambda) * np.power(lambda_array/lambda_v,delta)

    if return_AxAv:
        return A_lambda/Av

    #---- outputing extinguish/fraction ----#

    extinguish = np.power(10,-0.4*A_lambda)
    return extinguish

###########
# Testing #
###########

#---- testing basic functionality ----#

if False:

    import matplotlib.pyplot as plt

    lambda_array = np.linspace(1000, 22000, 1000) # AA

    Eb = 1
    plt.plot(lambda_array/1e+4, KC13(lambda_array, 0.25, delta=-0.3, Eb=Eb, Rv=4.05), label='delta=-0.3')
    plt.plot(lambda_array/1e+4, KC13(lambda_array, 0.25, delta= 0.0, Eb=Eb, Rv=4.05), label='delta=0.0')
    plt.plot(lambda_array/1e+4, KC13(lambda_array, 0.25, delta=+0.3, Eb=Eb, Rv=4.05), label='delta=0.3')
    plt.legend()
    plt.show()

#---- testing AxAv shape ----#

if False:

    import matplotlib.pyplot as plt

    lambda_array = np.linspace(1000, 20000, 1000) # AA
    AxAv = KC13(lambda_array, 10, delta=0, Eb=0, return_AxAv=True)

    plt.plot(lambda_array/1e+4, AxAv)
    plt.xscale('log')
    plt.xlim(0.1,2.0)
    plt.ylim(0,8)
    plt.show()
