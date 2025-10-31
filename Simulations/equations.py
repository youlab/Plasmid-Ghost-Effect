import numpy as np

def func(x,time,para):
    x[x<0]=0

    Ns=para[0]#number of species
    mu=para[1]#basic reproduction rate //size=(Ns,)
    eta=para[2]#conjugation rate //size=(Ns,Ns)
    kappa=para[3]#plasmid loss rate//size=(Ns,)
    D=para[4]#dilution rate//size=1
    alpha=para[5]#fitness / cost by plasmid//size=(Ns,)
    gamma=para[6]#competition between two species//size=(Ns,Ns)
    c=para[7]#capacity//size=(Ns,)
    A=para[8]

    f=x[0:Ns]
    p=x[Ns:]

    reproduction_f=alpha*mu*(1-np.dot(gamma,f+p)/c)*f
    reproduction_p=mu*(1-np.dot(gamma,f+p)/c)*p

    conjugation=f*np.dot(eta,p)

    death_f=D*f+A*f
    death_p=D*p

    loss=kappa*p

    dfdt=reproduction_f-conjugation-death_f+loss
    dpdt=reproduction_p+conjugation-death_p-loss

    dydt=np.hstack((dfdt,dpdt))
    return dydt