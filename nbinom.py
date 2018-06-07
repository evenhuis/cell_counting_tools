import scipy.stats as stats

def conv_mu_ab( mu, std_mu):
    ''' Converts the mean and stdev of the gamma distribution to
    alpha and beta'''
    alpha = 1./std_mu**2+1e-10
    beta  = alpha/mu
    return [alpha,beta]

def gamma_stdf_rvs( mu, std_mu, size=1 ):
    [alpha,beta]=conv_mu_ab(mu,std_mu)
    return stats.gamma.rvs(alpha,scale=1./beta,size=size)

def gamma_stdf_pdf( x, mu, std_mu ):
    [alpha,beta]=conv_mu_ab(mu,std_mu)
    return stats.gamma.pdf(x,alpha,scale=1./beta)


def conv_ab_rp( alpha, beta):
    ''' over-dispersed poission
     convolution of possion counts over a gamma dist 
        int_mu=0^inf ( gamma(mu, alpha,beta)  poission(mu) dmu )
    converts the alpha,beta for the gamma dist to the 
    r,p for the negative bionomiali
    see https://www.johndcook.com/negative_binomial.pdf
    '''
    r = alpha
    p = 1./(beta+1.)
    return [r,p]

def poisson_conv_gamma_stdf( mu, std_mu, s1=10,s2=25):
    ''' return counts from a poission rate mu  with from a gamm dist'''
    counts = np.zeros(s1*s2)
    for i in range(s1):
        mu_i = gamma_stdf_rvs( mu, std_mu ) # generate a rate
        counts[i*s2:(i+1)*s2] = stats.poisson.rvs(mu_i,size=s2)
    return counts

def nbinom_mu_stdf_rvs( mu, std_mu, size=1):
    [alpha,beta]=conv_mu_ab(mu,std_mu)
    [r,    p   ]=conv_ab_rp(alpha,beta)
    return stats.nbinom.rvs(r,1-p,size=size)     # need to convert the negative binom dist

def nbinom_mu_stdf_pmf( x, mu, std_mu):
    if( std_mu < 1e-7 ):
        return stats.poisson.pmf( x,mu)    
    else:
        [alpha,beta]=conv_mu_ab(mu,std_mu)
        [r,    p   ]=conv_ab_rp(alpha,beta)
        return stats.nbinom.pmf(x,r,1-p)             # need to convert the negative binom dist

def nbinom_mu_stdf_logpmf( x, mu, std_mu):
    if( std_mu < 1e-7 ):
        return stats.poisson.logpmf( x,mu)
    else:
        [alpha,beta]=conv_mu_ab(mu,std_mu)
        [r,    p   ]=conv_ab_rp(alpha,beta)
        return stats.nbinom.logpmf(x,r,1-p)         # need to convert the negative binom dist


if( __name__=="__main__"):
    import numpy as np
    import matplotlib.pyplot as plt

    mu   = 50.
    stdf = 0.25


    # generate the poission rate params
    bins = np.arange(0,100,1)

    plt.hist( gamma_stdf_rvs( mu, stdf/2., 10000 ),  normed=True, bins=bins)
    plt.hist( gamma_stdf_rvs( mu, stdf   , 10000 ),  normed=True, bins=bins)

    [alpha,beta] = conv_mu_ab(mu, stdf )
    print(alpha,beta)
    plt.plot( bins, stats.gamma.pdf( bins,alpha, scale=1./beta ))
    plt.show()

    plt.hist( poisson_conv_gamma_stdf( mu, stdf, s1=1000,s2=10), normed=True, bins=bins)
    plt.hist( nbinom_mu_stdf_rvs     ( mu, stdf,  10000      ),   normed=True, bins=bins)
    plt.plot( bins,nbinom_mu_stdf_pmf (bins, mu,stdf ))
    plt.plot( bins,stats.poisson.pmf(bins, mu) )
    plt.show()


