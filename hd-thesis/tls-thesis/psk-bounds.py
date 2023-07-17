# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

from sage.all import *
#from math import log

#
# SIGMA bounds (ours and Canetti-Krawczyk 2002)
#

# Computes our SIGMA advantage (in bits)
def ourSIGMAbound(t,qn,qs,qRO,logp,nl,kl):
        sidcoll = 3*qs**2/(2**(nl+1+logp))  # sid collision bound: 3qs^2/(2^(nl+1)*p)
        tB1 = t + 2*qRO*(logp)              # stDH reduction time tB1: t + qRO exponentiations (costing 2 * log(p) group operations each)
        stDHadv = 4*tB1**2 / (2**logp)      # stDH adv for B1: GGM bound 4*t^2/p  (where t = tB1)
        PRFadv = qs*qRO/2**kl               # PRF, idealized ROM bound [q_New = qs]: qs * qRO / 2^kl
        Sigadv = qn*(t**2)/(2**logp)        # Signature mu-EUF-CMA: qN loss to single-user (the latter idealized as ECDlog: t^2/p)
        MACadv = qs/2**(kl) + qs*qRO/2**kl  # MAC, idealized ROM bound [q_Vf = q_New = qs, q_Corr = 0, ol = kl]: qs/2^kl + qs*qRO/2^kl
        tot = sidcoll + stDHadv + PRFadv + Sigadv + MACadv
        return float(log(tot)/log(2))		# compute bit advantage


# Computes the CK SIGMA advantage (in bits)
def theirSIGMAbound(t,qn,qs,qRO,logp,nl,kl):
        sidcoll = (2*qs**2)/2**(nl+logp)    # sid collision bound: 2qs^2/(2^nl * p)
        Sigadv1 = qn*(t**2)/(2**logp)       # Signature1 mu-EUF-CMA: qN loss to single-user idealized ECDlog (t^2/p)
        DDHadv = t**2/(2**logp)             # DDH GGM bound: t^2/p
        PRFadv = qRO/2**kl                  # PRF, single-user idealized ROM bound: qRO / 2^kl
        Sigadv2 = t**2/(2**logp)            # Signature2, single-user idealized ECDlog (t^2/p)
        MACadv = 2/2**kl + qRO/2**kl        # MAC, single-user idealized ROM bound [q_New = 1, q_Vf = 2, q_Corr = 0, ol = kl]: 2/2^kl + qRO/2^kl
        tot = sidcoll + Sigadv1 + qn * qs * (DDHadv + PRFadv + (qn+1)*(Sigadv2) + MACadv)
        return float(log(tot)/log(2))


#
# TLS 1.3 bounds (ours and MSKE model [DFGS JoC 2021 / ePrint 2020/1044])
#

# Computes our TLS 1.3 advantage (in bits)
def ourTLSbound(t,qn,qs,qRO,logp,nl,kl):
        sidcoll = 3*qs**2/(2**(nl+1+logp))  # sid collision bound: 3qs^2/(2^(nl+1)*p)
### HD 3-9-2021 -- changing hash back to RO bound idealized Hash CR: t^2/2^(kl+1)
        Hashadv = qRO**2/(2**(kl+1)) +1/2**(kl)         # idealized Hash CR: qRO^2/2^(kl+1) + 1/2^(kl)
        tB2 = t+2*qRO*logp                  # stDH reduction B2 time: t + qRO exponentiations
        stDHadv = 4*tB2**2/(2**logp)        # stDH adv for B2: GGM bound 4*t^2/p (where t = tB2)
        ROcoll = qRO*qs/(2**(kl-1))         # RO collision bound: qRO*qs / 2^(kl-1)
        Sigadv = qn*(t**2)/(2**logp)        # Signature mu-EUF-CMA: qN loss to single-user (the latter idealized as ECDlog: t^2/p)
        MACadv = qs/(2**kl)+qs*qRO/(2**kl)  # MAC, idealized ROM bound [q_Vf = q_New = qs, q_Corr = 0, ol = kl]: qs/2^kl + qs*qRO/2^kl
        tot = sidcoll + Hashadv + 2*stDHadv + ROcoll + Sigadv + MACadv
        return float(log(tot)/log(2))

def ourPSKDHEbound(t, qn, qs, qRO, logp, nl, kl):
        tmaul = t+qs*kl # maul bound: no advantage increase, but adv time increase.
        indiff2 = 2*(12*qs+qRO)**2/(2**kl) # HMAC bound, qPriv = 12qSend, qPub = qRO
        indiff3 = 2*qRO**2/(2**kl) + 8*(qRO+36*qs)**2/(2**kl) # indiff bound hop 3, qPriv = 6qSend, qPub = qRO
        psk_guess = qRO*qn/(2**kl) # probability that an RO query collides with an unrevealed PSK
        sidcoll = 2*qs**2/(2**(nl+logp)) # sid collision bound: 2qs^2/2^(nl)*p
        psk_coll = qn**2/2**(kl)    # collision between two psks
        coll_res = (qRO+qs)**2/(2**(kl)) +(qRO+6*qs)**2/(2**(kl))+ 2/(2**(kl)) #binder collisions and Thash collisions
        MACadv = qs/(2**kl) # probability of guessing a random MAC tag
        tB4 = tmaul + 4*qRO*logp #stDH reduction time: t + qRO exponentiations
        stDHadv = 4*tB4**2/(2**logp)
        tot = indiff2 + indiff3 + psk_guess + sidcoll + psk_coll + coll_res + MACadv + stDHadv
        return float(log(tot)/log(2))

def ourPSKbound(t, qn, qs, qRO, nl, kl):
        tmaul = t+ qs*kl # maul bound: no advantage increase, but adv time increase.
        indiff2 = 2*(12*qs+qRO)**2/(2**kl) # HMAC bound, qPriv = 12qSend, qPub = qRO
        indiff3 = 2*qRO**2/(2**kl) + 8*(qRO+36*qs)**2/(2**kl) # indiff bound hop 3, qPriv = 6qSend, qPub = qRO
        psk_guess = qRO*qn/(2**kl) # probability that an RO query collides with an unrevealed PSK
        sidcoll = qs**2/(2**nl) # sid collision bound: 2qs^2/2^(nl) (no group elements)
        psk_coll = qn**2/(2**kl)    # collision between two psks
        coll_res = (qRO+qs)**2/(2**(kl))+(qRO+6*qs)**2/(2**(kl)) + 2/(2**kl) # binder collisions and Thash collisions
        MACadv = qs/(2**kl) # probability of guessing a random MAC tag
        tot = indiff2 + indiff3 + psk_guess + sidcoll + psk_coll + coll_res + MACadv
        return float(log(tot)/log(2))

def theirPSKbound(t, qn, qs, qRO, nl, kl):
        match =  qRO**2/(2**kl) + qn**2/(2**kl) + qs**2/(2**nl) # Match bound: coll(HMAC) + psk collision + nonce collision
        Hashadv = qRO**2/(2**(kl+1))+1/(2**kl)
        PRFadv = qRO/(2**kl)
        tot = match + 8*qs * (Hashadv + qn*8*PRFadv)
        return float(log(tot)/log(2))

def theirPSKDHEbound(t,qn,qs,qRO, logp, nl,kl):
        match =  qRO**2/(2**kl) + qn**2/(2**kl) + qs**2/(2**(nl+logp)) # Match bound: coll(HMAC) + psk collision + sid collision
        Hashadv = qRO**2/(2**(kl+1))+1/(2**kl)
        PRFadv = qRO/(2**kl)
        MACadv = qs/(2**kl)+qs*qRO/(2**kl) 
        PRFODHadv = 4*t**2/(2**logp)        # PRF-ODH with qRO queries == stDH with qRO queries ==> GGM bound 4 * t^2/p, NOTE: we don't consider extra exponentiations here, so t
        tot = match + 8*qs * (Hashadv + qs* qn*8*PRFadv + qs*qn*2*MACadv + qs* (PRFODHadv + 5* PRFadv))
        return float(log(tot)/log(2))



# Computes the MSKE TLS 1.3 advantage (in bits)
def theirTLSbound(t,qn,qs,qRO,logp,nl,kl):
        match = qs**2 / (2**(nl + logp))    # Match bound: qs^2 / (2^nl * p)
### HD 3-9-2021 -- changing hash back to RO bound idealized Hash CR: t^2/2^(kl+1)
        Hashadv = qRO**2/2**(kl+1)+1/2**kl            # idealized Hash CR: qRO^2/2^(kl+1)+1/2**kl
        Sigadv = t**2/(2**logp)             # Signature, single-user: ECDlog (t^2/p)
        PRFODHadv = 4*t**2/(2**logp)        # PRF-ODH with qRO queries == stDH with qRO queries ==> GGM bound 4 * t^2/p, NOTE: we don't consider extra exponentiations here, so t
        PRFadv = qRO/(2**kl)                # PRF, idealized single-user ROM bound: qRO / 2^kl
        tot = match + qs * (Hashadv + qn*Sigadv + qs*(PRFODHadv + 5*PRFadv))
        return float(log(tot)/log(2))



# curve: [name, p, security level, symkey length]
# TLS 1.3: sym key length is either 256 bits (SHA-256, for curve bit-security 128) or 384 bits (SHA-384, for higher curve bit-security)
curves = [["secp256r1", 256, 128, 256], ["secp384r1", 384, 192, 384], ["secp521r1", 521, 256, 384], ["x25519", 252, 128, 256], ["x448", 446, 224, 384]]

highcolor = "red!25"
warncolor = "orange!25"
meetcolor = "green!25"



#
# Generate detailed LaTeX tables with numercial bounds.
#   fixedkeylength = 0   ---> Table 2: Numerical bounds for SIGMA and TLS 1.3, for symmetric lengths twice the curve’s bit-security.
#   fixedkeylength = 1   ---> Table 3: Numerical bounds for SIGMA and TLS 1.3, for symmetric lengths of 256 bits.
#   fixedkeylength = 2   ---> Table new: Numerical bounds for SIGMA and TLS 1.3, for symmetric lengths according to curves vector (256bits for sec-level 128, 384bits otherwise).
#
def genLatexbounds(fixedkeylength):

        curves = [["secp256r1", 256, 128, 256], ["x25519", 252, 128, 256], ["secp384r1", 384, 192, 384], ["x448", 446, 224, 384], ["secp521r1", 521, 256, 384], ]
        minDHEdiff = math.inf
        maxDHEdiff = 0
        minPSKdiff = math.inf
        maxPSKdiff = 0
        for t in [40, 60, 80]:
                for c in curves:
                        # symmetric key length: fixed / 2 * curve-bitsecurity / curve-dependend 256 or 384

                        symkl = c[3]
                        symnl = 256     # TLS 1.3 nonce length is fixed to 256 bit                                 

                        for qn in [25, 35]:
                                for qs in [35, 45,55]:
                                    if qs > t:
                                        continue
                                    print ("$2^{{{:d}}}$\t&".format(t),end='')
                                    print ("$2^{{{:d}}}$\t&".format(qn),end='')
                                    print ("$2^{{{:d}}}$\t&".format(qs),end='')
                                    print ("$2^{{{:d}}}$\t&".format(t-10),end='')
                                    print ("\\texttt{{{:s}}} ($b \\!=\\! {:d}$, \\! $p \\!\\approx\\! 2^{{{:d}}}$)\t&".format(c[0], c[2], c[1]),end='')
                                    print ("$2^{{{:d}}}$\t&".format(t-c[2]),end='')

            #t,    qn,    qs,       qRO, logp,    nl,    kl
      
                                    DHEtheir = int(round(theirPSKDHEbound(2**t, 2**qn, 2**qs, 2**(t-10), c[1], symnl, symkl)))
                                    DHEour     =   int(round(  ourPSKDHEbound(2**t, 2**qn, 2**qs, 2**(t-10), c[1], symnl, symkl)))
                                    if (DHEtheir - DHEour) < minDHEdiff:
                                        minDHEdiff = DHEtheir-DHEour 
                                    if (DHEtheir - DHEour) > maxDHEdiff:
                                        maxDHEdiff = DHEtheir-DHEour 
        
                                    col = "&"
                                    for adv in [DHEtheir, DHEour]:                                             
                                        #if adv >= 0:
                                        #        print("\\cellcolor{{{:s}}}1\t\t\t{:s}".format(highcolor,col),end="")
                                        #elif t - adv < c[2]:
                                        #        print ("\\cellcolor{{{:s}}}$\\approx 2^{{{:d}}}$\t{:s}".format(highcolor,adv,col),end="")
                                        if t - adv > c[2]:
                                                print ("\\cellcolor{{{:s}}}$\\approx 2^{{{:d}}}$\t{:s}".format(meetcolor,adv,col),end="")
                                        elif t - adv == c[2]:
                                                print ("\\cellcolor{{{:s}}}$\\approx 2^{{{:d}}}$\t{:s}".format(warncolor,adv,col),end="")
                                        else:
                                                print ("$\\approx 2^{{{:d}}}$\t{:s}".format(adv,col),end=" ")
                                        col = ""
                                    print( "\\\\")
                        print("\\midrule")
                print("\\midrule")

        print ("\\bottomrule")
        
        print("\n\n\n\n\n")

        for t in [40, 60, 80]:
                        print("\\midrule")
                        symkl = 256 
                        symnl = 256     # TLS 1.3 nonce length is fixed to 256 bit                                 
                        for qn in [25, 35]:
                                for qs in [35, 45,55]:
                                    if qs > t:
                                        continue
                                    print ("$2^{{{:d}}}$\t&".format(t),end='')
                                    print ("$2^{{{:d}}}$\t&".format(qn),end='')
                                    print ("$2^{{{:d}}}$\t&".format(qs),end='')
                                    print ("$2^{{{:d}}}$\t&".format(t-10), end='')
                                    print ("$2^{{{:d}}}$\t&".format(t-128),end='')

                                    #t,    qn,    qs,       qRO,    nl,    kl
                                    PSKtheir = int(round(theirPSKbound(2**t, 2**qn, 2**qs, 2**(t-10), symnl, symkl)))
                                    PSKour     =   int(round(  ourPSKbound(2**t, 2**qn, 2**qs, 2**(t-10), symnl, symkl))) 
                                    if (PSKtheir - PSKour) < minPSKdiff:
                                        minPSKdiff = PSKtheir-PSKour 
                                    if (PSKtheir - PSKour) > maxPSKdiff:
                                        maxPSKdiff = PSKtheir-PSKour 
                                    col = "&"
                                    for adv in [PSKtheir, PSKour]:                                             
                                        #if adv >= 0:
                                        #        print("\\cellcolor{{{:s}}}1\t\t\t{:s}".format(highcolor,col),end="")
                                        #elif t - adv < 128:
                                        #        print ("\\cellcolor{{{:s}}}$\\approx 2^{{{:d}}}$\t{:s}".format(highcolor,adv,col),end="")
                                        if t - adv > 128:
                                                print ("\\cellcolor{{{:s}}}$\\approx 2^{{{:d}}}$\t{:s}".format(meetcolor,adv,col),end="")
                                        elif t - adv == 128:
                                                print ("\\cellcolor{{{:s}}}$\\approx 2^{{{:d}}}$\t{:s}".format(warncolor,adv,col),end="")
                                        else:
                                                print ("$\\approx 2^{{{:d}}}$\t{:s}".format(adv,col),end=" ")
                                        col = ""
                                    print( "\\\\")                                
        print ("\\bottomrule")

        print("\n\n\n\n\n")
        print("min DHE:",minDHEdiff)
        print("max DHE:",maxDHEdiff)
        print("min PSK:",minPSKdiff)
        print("max PSK:",maxPSKdiff)

   
genLatexbounds(2)   ## new Table (sym key length 256 or 384 bits, nonce length 256 bits)

