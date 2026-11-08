#-------------------------------------------------------------------------------
#
# Error spending group sequential tests
#
# Three routines for different forms of 1-sided error spending test
#
#   ersp1a,  with no futility boundary at all
#
#   ersp1b,  with a NON-BINDING lower, futility boundary
#
#   ersp1c,  with a BINDING lower, futility boundary
#
# One routine for a 2-sided error spending test
#
#   ersp2,   with BINDING lower and upper boundaries (implemented
#            by a call to ersp1c)
#
# Routines to find:
#
#   the value of rho which gives a converging boundary for a given
#   sequence of information levels
#
#     find_rho_ersp1a, find_rho_ersp1b, find_rho_ersp1c, find_rho_ersp2
#
#   the value of rfac which gives a converging boundary for a given
#   value of rho
#
#     find_rfac_ersp1a, find_rfac_ersp1b, find_rfac_ersp1c, find_rfac_ersp2
#
#-------------------------------------------------------------------------------

ersp1a = function(r,na,inf,cume)
{
 # 1-sided error spending test of H0: theta =< 0 vs theta > 0, with
 # no stopping to accept H0 (lower, futility boundary is set at Z = -50).
 #
 # The same result should be obtained from ersp1b or ersp1c if all cumulative
 # type 2 error probabilities are set to be zero.
 #
 # IN:
 #
 #   r             multiplier for Simpson's rule grid, e.g., r=16
 #   na            maximum number of analyses
 #   inf[1:na]     information levels
 #   cume[2;1:na]  cumulative type 1 error probabilities
 #   cume[1;1:na]  cumulative type 2 error probabilities (all assumed = 0)
 #
 # OUT:
 #   [1] zbdy[1;1:na]  lower boundary for Z (an approximation to -infinity)
 #       zbdy[2;1:na]  upper boundary for Z
 #   [2] pu            P(Exit upper boundary under theta=0)
 #   [3] pl            P(Exit lower boundary under theta=0)
 #
 # NOTES:
 #
 # inf
 #
 # Information for theta at each analysis is given in inf[1:na],
 # so at analaysis k,
 #
 #  Z ~ N(theta*sqrt(inf(k)),inf(k)).
 #
 # Values in inf should be a strictly increasing positive sequence.
 # If this is not the case, an error will result.
 #
 # cume
 #
 # The labelling is because cume[2,1:na] are used primarily
 # in finding zbdy[2,1:na] and cume[1,1:na] are used primarily in
 # finding zbdy[1,1:na].
 #
 # Cumulative type I error probabilities are specified in cume[2;1:na].
 # These values should form an increasing positive sequence. Negative or
 # decreasing values will result in an upper boundary at a large positive
 # value to approximate +infinity (spending approximately zero type I error)
 #
 # If there were a futility boundary, cumulative type II error
 # probabilities would be specified in cume[1;1:na]. However, these
 # are assumed to be zero and lower boundary is set at a large
 # negative value to approximate -infinity.

 zbdy=array(-50,c(2,na))
 epsilon=0.00001

 theta=0   # Note: all calculations are under the null hypothesis, theta=0
 pl=0
 pu=0

 k1=0
 inf1=0
 m1=1
 z1=array(0, c(1))
 h1=array(1, c(1))
 k2=1

 stmesh=mesh(r)     # Standard mesh for numerical integration
 minmesh=stmesh[1]
 maxmesh=stmesh[length(stmesh)]

 while(k2 <= na)
 {
  inf2=inf[k2]
  dinf=inf2-inf1
  z2var=1-inf1/inf2
  z2sd=sqrt(z2var)
  ez2=array(0,c(m1))
  ez2=z1*sqrt(inf1/inf2)+theta*(dinf/sqrt(inf2))

  # Find upper boundary point

  fstar=cume[2,k2]-pu

  bhi=maxmesh
  fhi=sum((1-pnorm((bhi-ez2)/z2sd))*h1)
  if(fhi >= fstar)
  {
   zbdy[2,k2]=bhi
  }else
  {
   blo=minmesh
   flo=sum((1-pnorm((blo-ez2)/z2sd))*h1)
   if(flo <= fstar)
   {
    zbdy[2,k2]=blo
   }else
   {
    # print(c("A: blo, flo, bhi, fhi= "))
    # print(c(blo,flo,bhi,fhi))
    while(flo-fhi > epsilon)
    {
     bmid=0.5*(blo+bhi)
     fmid=sum((1-pnorm((bmid-ez2)/z2sd))*h1)
     if(fmid < fstar)
     {
      bhi=bmid
      fhi=fmid
     }else
     {
      blo=bmid
      flo=fmid
     }
    }
    bend=blo+(bhi-blo)*(fstar-flo)/(fhi-flo)
    fend=sum((1-pnorm((bend-ez2)/z2sd))*h1)
    zbdy[2,k2]=bend
   }
  } 
  if(k2 == na)
  {
   zbdy[1,k2]=zbdy[2,k2]
  }

  pupper=sum((1-pnorm((zbdy[2,k2]-ez2)/z2sd))*h1)
  pu=pu+pupper
  plower=sum(pnorm((zbdy[1,k2]-ez2)/z2sd)*h1)
  pl=pl+plower

  if(k2 < na)
  {
   zmesh=theta*sqrt(inf2)+stmesh                             # Mesh, analysis k2
   zmeshb=zmesh[zmesh >= zbdy[1,k2] & zmesh <= zbdy[2,k2]]   # Modify mesh
   if(min(zmesh) < zbdy[1,k2]) {zmeshb=c(zbdy[1,k2],zmeshb)}
   if(max(zmesh) > zbdy[2,k2]) {zmeshb=c(zmeshb,zbdy[2,k2])}
   len=length(zmeshb)
   m2=2*len-1
   z2=array(0,c(m2))
   for(i in 1:len) {z2[2*i-1]=zmeshb[i]}
   if(len > 1)
   {
    for(i in 1:(len-1)) {z2[2*i]=0.5*(z2[2*i-1]+z2[2*i+1])}
   }
   m2=2*len-1                                           # Simpson's rule weights
   w2=array(0, c(m2))
   w2[1]=0
   if(m2>1)
   {
    w2[1]=(z2[3]-z2[1])/6
    w2[m2]=(z2[m2]-z2[m2-2])/6
    if(m2>3)
    {
     for(i2 in seq(3,m2-2,2)) {w2[i2]=(z2[i2+2]-z2[i2-2])/6}
    }
    for(i2 in seq(2,m2-1,2)) {w2[i2]=(z2[i2+1]-z2[i2-1])*4/6}
   }
   h2=array(0, c(m2))
   for(i2 in 1:m2)
   {
    h2[i2]=sum(dnorm(z2[i2]-ez2,sd=z2sd)*w2[i2]*h1)
   }
   k1=k2
   inf1=inf2
   m1=m2
   z1=z2
   h1=h2
  }

  k2=k2+1
 }

 list(
      zbdy,
      pu,
      pl
     )
}

#-------------------------------------------------------------------------------

ersp1b = function(r,na,inf,delta,cume)
{
 # 1-sided error spending test of H0: theta =< 0 vs theta > 0
 # with an upper, efficacy boundary and a NON-BINDING lower,
 # futility boundary. 
 #
 # IN:
 #
 #   r             multiplier for Simpson's rule grid, e.g., r=16
 #   na            maximum number of analyses
 #   inf[1:na]     information levels
 #   delta         value of theta at which type II error is calculated
 #   cume[2;1:na]  cumulative type 1 error probabilities under theta=0
 #   cume[1;1:na]  cumulative type 2 error probabilities under theta=delta
 #
 # OUT:
 #   [1] zbdy[1;1:na]  lower boundary
 #       zbdy[2;1:na]  upper boundary
 #   [2] pu            P(Exit upper boundary under theta=0)
 #   [3] pl            P(Exit lower boundary under theta=delta)
 #
 # NOTES:
 #
 # inf
 #
 # Information for theta at each analysis is given in inf[1:na],
 # so at analaysis k,
 #
 #  Z ~ N(theta*sqrt(inf(k)),inf(k)).
 #
 # Values in inf should be a strictly increasing positive sequence.
 # If this is not the case, an error will result.
 #
 # cume
 #
 # The labelling is because cume[2,1:na] are used primarily
 # in finding zbdy[2,1:na] and cume[1,1:na] are used primarily in
 # finding zbdy[1,1:na].
 #
 # Cumulative type I error probabilities are specified in cume[2;1:na].
 # These values should form an increasing positive sequence. Negative or
 # decreasing values will result in an upper boundary at a large positive
 # value to approximate +infinity (spending approximately zero type I error)
 #
 # Cumulative type II error probabilities are specified in cume[2;1:na].
 # These values should form an increasing positive sequence. Negative or
 # decreasing values will result in an upper boundary at a large positive
 # value to approximate +infinity (spending approximately zero type I error)
 #
 # The binding upper boundary, zbdy[2,1:na] is calculated first.
 # Then, the non-binding lower boundary, zbdy[1,1:na] is found. 
 #
 # With too little information: The lower and upper boundaries do not
 # meet up. To use these in practice, take zbdy[na,2] as the final
 # critical value, and ignore zbdy[1,na].
 #
 # With excess information: The full type II error probability, under
 # theta=delta, cannot be spent. Assuming the futility boundary is
 # followed: at some analysis k*, zbdy[1,k*]=zbdy[2,k*] and the trial
 # terminates. For completeness, we set zbdy[1,k]=zbdy[2,k] for k > k*.

 zbdy=array(-50,c(2,na))
 epsilon=0.00001

 # Upper boundary -- calculate all of zbdy[2,1:na] first.

 theta=0
 pu=0

 k1=0
 inf1=0
 m1=1
 z1=array(0, c(1))
 h1=array(1, c(1))
 k2=1

 stmesh=mesh(r)     # Standard mesh for numerical integration
 minmesh=stmesh[1]
 maxmesh=stmesh[length(stmesh)]

 while(k2 <= na)
 {
  inf2=inf[k2]
  dinf=inf2-inf1
  z2var=1-inf1/inf2
  z2sd=sqrt(z2var)
  ez2=array(0,c(m1))
  ez2=z1*sqrt(inf1/inf2)+theta*(dinf/sqrt(inf2))

  # Find upper boundary point

  fstar=cume[2,k2]-pu
  bhi=maxmesh
  fhi=sum((1-pnorm((bhi-ez2)/z2sd))*h1)
  if(fhi >= fstar)
  {
   zbdy[2,k2]=bhi
  }else
  {
   blo=minmesh
   flo=sum((1-pnorm((blo-ez2)/z2sd))*h1)
   if(flo <= fstar)
   {
    zbdy[2,k2]=blo
   }else
   {
    while(flo-fhi > epsilon)
    {
     bmid=0.5*(blo+bhi)
     fmid=sum((1-pnorm((bmid-ez2)/z2sd))*h1)
     if(fmid < fstar)
     {
      bhi=bmid
      fhi=fmid
     }else
     {
      blo=bmid
      flo=fmid
     }
    }
    bend=blo+(bhi-blo)*(fstar-flo)/(fhi-flo)
    zbdy[2,k2]=bend
   }
  }

  pupper=sum((1-pnorm((zbdy[2,k2]-ez2)/z2sd))*h1)
  pu=pu+pupper

  if(k2 < na)
  {
   zmesh=theta*sqrt(inf2)+stmesh                             # Mesh, analysis k2
   zmeshb=zmesh[zmesh >= zbdy[1,k2] & zmesh <= zbdy[2,k2]]   # Modify mesh
   if(min(zmesh) < zbdy[1,k2]) {zmeshb=c(zbdy[1,k2],zmeshb)}
   if(max(zmesh) > zbdy[2,k2]) {zmeshb=c(zmeshb,zbdy[2,k2])}
   len=length(zmeshb)
   m2=2*len-1
   z2=array(0,c(m2))
   for(i in 1:len) {z2[2*i-1]=zmeshb[i]}
   if(len > 1)
   {
    for(i in 1:(len-1)) {z2[2*i]=0.5*(z2[2*i-1]+z2[2*i+1])}
   }
   m2=2*len-1                                           # Simpson's rule weights
   w2=array(0,c(m2))
   w2[1]=0
   if(m2>1)
   {
    w2[1]=(z2[3]-z2[1])/6
    w2[m2]=(z2[m2]-z2[m2-2])/6
    if(m2>3)
    {
     for(i2 in seq(3,m2-2,2)) {w2[i2]=(z2[i2+2]-z2[i2-2])/6}
    }
    for(i2 in seq(2,m2-1,2)) {w2[i2]=(z2[i2+1]-z2[i2-1])*4/6}
   }
   h2=array(0,c(m2))
   for(i2 in 1:m2)
   {
    h2[i2]=sum(dnorm(z2[i2]-ez2,sd=z2sd)*w2[i2]*h1)
   }
   k1=k2
   inf1=inf2
   m1=m2
   z1=z2
   h1=h2
  }

  k2=k2+1
 }

 # Lower boundary

 theta=delta
 pl=0

 k1=0
 inf1=0
 m1=1
 z1=array(0, c(1))
 h1=array(1, c(1))
 k2=1

 while(k2 <= na)
 {
  inf2=inf[k2]
  dinf=inf2-inf1
  z2var=1-inf1/inf2
  z2sd=sqrt(z2var)
  ez2=array(0,c(m1))
  ez2=z1*sqrt(inf1/inf2)+theta*(dinf/sqrt(inf2))

  # Find lower boundary point

  fstar=cume[1,k2]-pl
  # print(c("fstar= ",fstar))

  bhi=zbdy[2,k2]
  fhi=sum(pnorm((bhi-ez2)/z2sd)*h1)
  if(fhi <= fstar)
  {
   zbdy[1,k2]=bhi
  }else
  {
   blo=minmesh
   flo=sum(pnorm((blo-ez2)/z2sd)*h1)
   if(flo >= fstar)
   {
    zbdy[1,k2]=blo
   }else
   {
    while(fhi-flo > epsilon)
    {
     bmid=0.5*(blo+bhi)
     fmid=sum(pnorm((bmid-ez2)/z2sd)*h1)
     if(fmid < fstar)
     {
      blo=bmid
      flo=fmid
     }else
     {
      bhi=bmid
      fhi=fmid
     }
    }
    bend=blo+(bhi-blo)*(fstar-flo)/(fhi-flo)
    fend=sum(pnorm((bend-ez2)/z2sd)*h1)
    zbdy[1,k2]=bend
   }
  }

  plower=sum(pnorm((zbdy[1,k2]-ez2)/z2sd)*h1)
  pl=pl+plower

  if(k2 < na)
  {
   zmesh=theta*sqrt(inf2)+stmesh                             # Mesh, analysis k2
   zmeshb=zmesh[zmesh >= zbdy[1,k2] & zmesh <= zbdy[2,k2]]   # Modify mesh
   if(min(zmesh) < zbdy[1,k2]) {zmeshb=c(zbdy[1,k2],zmeshb)}
   if(max(zmesh) > zbdy[2,k2]) {zmeshb=c(zmeshb,zbdy[2,k2])}
   len=length(zmeshb)
   m2=2*len-1
   z2=array(0,c(m2))
   for(i in 1:len) {z2[2*i-1]=zmeshb[i]}
   if(len > 1)
   {
    for(i in 1:(len-1)) {z2[2*i]=0.5*(z2[2*i-1]+z2[2*i+1])}
   }
   m2=2*len-1                                           # Simpson's rule weights
   w2=array(0, c(m2))
   w2[1]=0
   if(m2>1)
   {
    w2[1]=(z2[3]-z2[1])/6
    w2[m2]=(z2[m2]-z2[m2-2])/6
    if(m2>3)
    {
     for(i2 in seq(3,m2-2,2)) {w2[i2]=(z2[i2+2]-z2[i2-2])/6}
    }
    for(i2 in seq(2,m2-1,2)) {w2[i2]=(z2[i2+1]-z2[i2-1])*4/6}
   }
   h2=array(0, c(m2))
   for(i2 in 1:m2)
   {
    h2[i2]=sum(dnorm(z2[i2]-ez2,sd=z2sd)*w2[i2]*h1)
   }
   k1=k2
   inf1=inf2
   m1=m2
   z1=z2
   h1=h2
  }

  k2=k2+1
 }

list(
     zbdy,
     pu,
     pl
     )
}

#-------------------------------------------------------------------------------

ersp1c = function(r,na,inf,delta,cume)
{
 # 1-sided error spending test of H0: theta =< 0 vs theta > 0
 # with an upper, efficacy boundary and a BINDING lower, futility
 # boundary.
 #
 # IN:
 #
 #   r             multiplier for Simpson's rule grid, e.g., r=16
 #   na            maximum number of analyses
 #   inf[1:na]     information levels
 #   delta         value of theta at which type II error is calculated
 #   cume[2;1:na]  cumulative type 1 error probabilities under theta=0 
 #   cume[1;1:na]  cumulative type 2 error probabilities under theta=delta
 #
 # OUT:
 #
 #   [1] zbdy[1;1:na]  lower boundary
 #       zbdy[2;1:na]  upper boundary
 #   [2] pu            P(Exit upper boundary under theta=0)
 #   [3] pl            P(Exit lower boundary under theta=delta)
 #
 # NOTES:
 #
 # inf
 #
 # Information for theta at each analysis is given in inf[1:na],
 # so at analaysis k,
 #
 #  Z ~ N(theta*sqrt(inf(k)),inf(k)).
 #
 # Values in inf should be a strictly increasing positive sequence.
 # If this is not the case, an error will result.
 #
 # cume
 #
 # The labelling is because cume[2,1:na] are used primarily
 # in finding zbdy[2,1:na] and cume[1,1:na] are used primarily in
 # finding zbdy[1,1:na].
 #
 # Cumulative type I error probabilities are specified in cume[2;1:na].
 # These values should form an increasing positive sequence with
 # cume[2,na]=alpha and cume[1,na]=beta. Negative or decreasing values
 # will result in an upper boundary at a large positive value to
 # approximate +infinity (spending approximately zero type I error).
 #
 # Cumulative type II error probabilities are specified in cume[2;1:na].
 # These values should form an increasing positive sequence. Negative or
 # decreasing values will result in a lower boundary at a large negative
 # value to approximate -infinity (spending approximately zero type II error).
 #
 # Upper and lower boundaries are computed together.
 # For each k, the upper boundary, zbdy[2,1:na] is calculated first.
 # Then, the lower boundary, zbdy[1,1:na] is found. 
 #
 # With too little information: The lower and upper boundaries do not
 # meet up. To use these in practice, take zbdy[2,na] as the final
 # critical value, and ignore zbdy[1,na].
 #
 # With excess information: The full type II error probability, under
 # theta=delta, cannot be spent. At some analysis k* < na,
 # zbdy[1,k*]=zbdy[2,k*] and the trial will terminate. For completeness,
 # we assign large negative values to zbdy[1,k] and zbdy[2,k] for k > k*
 # -- these follow from trying to spend more type I error probability in
 # setting zbdy[2,k] and not managing to do so as there is no probabiity
 # associated with continuing paths -- and then zbdy[1,k] is set equal
 # to zbdy[2,k].
 #
 # A technical point:
 #
 # Define alpha=cume[2,na]. Suppose we set zbdy[1,k] and zbdy[2,k]
 # to spend the specified cumulative error probabilities, but the
 # probability of continuing paths under theta=0 is less than
 # alpha-cume[2,k]: then it will not be possible to bring the total
 # type I probability up to alpha in the remainder of the trial.
 # This situation occurs precisely if the trial can be stopped at
 # analysis k with cumulative type I error equal to alpha and
 # cumulative type II error probability less than or equal to cume[1,k].
 # In this situation it seems reasonable to insist on termination at
 # analysis k, setting the upper boundary zbdy[2,k] to make the cumulative
 # type I error probability equal to alpha, then setting zbdy[1,k]=zbdy[2,k].
 # This approach protects overall type I and II error rates. Stopping to
 # reject H0, rather than continue sampling, can only increase power;
 # if k=na-1, continuing to analysis k would give zbdy[2,na]=-infinity
 # so one always rejects H0 and the extra stage has no effect on power.
 # We have done this.
 #
 # One could take the above argument a step further and terminate at
 # analysis k if it is possible to stop with cumulative type I error
 # equal to alpha and cumulative type II error probability less than or
 # equal to beta. However, in the extra cases where this criterion leads
 # to earlier stopping, additional data would provide additional power.
 # We have not implemented this further step.

 zbdy=array(-50,c(2,na))
 epsilon=0.00001

 alpha=cume[2,na]
 istop=0

 pu=0
 pl=0

 k1=0
 inf1=0
 m1=1
 z1=array(0, c(1))
 hu1=array(1, c(1))
 hl1=array(1, c(1))
 k2=1

 stmesh=mesh(r)     # Standard mesh for numerical integration
 minmesh=stmesh[1]
 maxmesh=stmesh[length(stmesh)]

 while(k2 <= na & istop==0)
 {
  inf2=inf[k2]
  dinf=inf2-inf1
  z2var=1-inf1/inf2
  z2sd=sqrt(z2var)
  ezu2=z1*sqrt(inf1/inf2)
  ezl2=z1*sqrt(inf1/inf2)+delta*(dinf/sqrt(inf2))

  # Find upper boundary point, under theta=0

  # (i) Check if it is possible to spend all remaining type I error
  # probability now and, if so, do that here.

  fstar=alpha-pu
  bhi=maxmesh
  fhi=sum((1-pnorm((bhi-ezu2)/z2sd))*hu1)
  if(fhi >= fstar)
  {
   zbdyall=bhi
  }else
  {
   blo=minmesh
   flo=sum((1-pnorm((blo-ezu2)/z2sd))*hu1)
   if(flo <= fstar)
   {
    zbdyall=blo
   }else
   {
    while(flo-fhi > epsilon)
    {
     bmid=0.5*(blo+bhi)
     fmid=sum((1-pnorm((bmid-ezu2)/z2sd))*hu1)
     if(fmid < fstar)
     {
      bhi=bmid
      fhi=fmid
     }else
     {
      blo=bmid
      flo=fmid
     }
    }
    bend=blo+(bhi-blo)*(fstar-flo)/(fhi-flo)
    zbdyall=bend
   }
  }
  ff=sum(pnorm((zbdyall-ezl2)/z2sd)*hl1)
  if(ff <= cume[1,k2]-pl)
  {

   # Calclations shows it is possible to terminate at analysis k2,
   # spending all remaining type I error probability and cumulative
   # type II error probability no greater than cume[1,k2].
   # So, do this:

   zbdy[2,k2]=zbdyall
   zbdy[1,k2]=zbdy[2,k2]
   pupper=sum((1-pnorm((zbdy[2,k2]-ezu2)/z2sd))*hu1)
   pu=pu+pupper
   plower=sum(pnorm((zbdy[1,k2]-ezl2)/z2sd)*hl1)
   pl=pl+plower  
   istop=1
  }

  if(istop==0)
  {
   # (ii) Spend type I error probability cume[2,k2] knowing that, under
   # theta=0, the probability of continuing is > alpha-cume[2,k2].

   fstar=cume[2,k2]-pu
   bhi=maxmesh
   fhi=sum((1-pnorm((bhi-ezu2)/z2sd))*hu1)
   if(fhi >= fstar)
   {
    zbdy[2,k2]=bhi
   }else
   {
    blo=minmesh
    flo=sum((1-pnorm((blo-ezu2)/z2sd))*hu1)
    if(flo <= fstar)
    {
     zbdy[2,k2]=blo
   }else
    {
     while(flo-fhi > epsilon)
     {
      bmid=0.5*(blo+bhi)
      fmid=sum((1-pnorm((bmid-ezu2)/z2sd))*hu1)
      if(fmid < fstar)
      {
       bhi=bmid
       fhi=fmid
      }else
      {
       blo=bmid
       flo=fmid
      }
     }
     bend=blo+(bhi-blo)*(fstar-flo)/(fhi-flo)
     zbdy[2,k2]=bend
    }
   }

   pupper=sum((1-pnorm((zbdy[2,k2]-ezu2)/z2sd))*hu1)
   pu=pu+pupper

   # Find lower boundary point

   fstar=cume[1,k2]-pl
   # print(c("fstar= ",fstar))

   bhi=zbdy[2,k2]
   fhi=sum(pnorm((bhi-ezl2)/z2sd)*hl1)
   if(fhi <= fstar)
   {
    zbdy[1,k2]=bhi
   }else
   {
    blo=minmesh
    flo=sum(pnorm((blo-ezl2)/z2sd)*hl1)
    if(flo >= fstar)
    {
     zbdy[1,k2]=blo
    }else
    {
     while(fhi-flo > epsilon)
     {
      bmid=0.5*(blo+bhi)
      fmid=sum(pnorm((bmid-ezl2)/z2sd)*hl1)
      if(fmid < fstar)
      {
       blo=bmid
       flo=fmid
      }else
      {
       bhi=bmid
       fhi=fmid
      }
     }
     bend=blo+(bhi-blo)*(fstar-flo)/(fhi-flo)
     fend=sum(pnorm((bend-ezl2)/z2sd)*hl1)
     zbdy[1,k2]=bend
    }
   }

   plower=sum(pnorm((zbdy[1,k2]-ezl2)/z2sd)*hl1)
   pl=pl+plower
  }

  if(k2 < na & istop==0)
  {
   zmesh=stmesh                                              # Mesh, analysis k2
   zmeshb=zmesh[zmesh >= zbdy[1,k2] & zmesh <= zbdy[2,k2]]   # Modify mesh
   if(min(zmesh) < zbdy[1,k2]) {zmeshb=c(zbdy[1,k2],zmeshb)}
   if(max(zmesh) > zbdy[2,k2]) {zmeshb=c(zmeshb,zbdy[2,k2])}
   len=length(zmeshb)
   m2=2*len-1
   z2=array(0,c(m2))
   for(i in 1:len) {z2[2*i-1]=zmeshb[i]}
   if(len > 1)
   {
    for(i in 1:(len-1)) {z2[2*i]=0.5*(z2[2*i-1]+z2[2*i+1])}
   }
   m2=2*len-1                                           # Simpson's rule weights
   w2=array(0, c(m2))
   w2[1]=0
   if(m2>1)
   {
    w2[1]=(z2[3]-z2[1])/6
    w2[m2]=(z2[m2]-z2[m2-2])/6
    if(m2>3)
    {
     for(i2 in seq(3,m2-2,2)) {w2[i2]=(z2[i2+2]-z2[i2-2])/6}
    }
    for(i2 in seq(2,m2-1,2)) {w2[i2]=(z2[i2+1]-z2[i2-1])*4/6}
   }
   hu2=array(0, c(m2))
   hl2=array(0, c(m2))
   for(i2 in 1:m2)
   {
    hu2[i2]=sum(dnorm(z2[i2]-ezu2,sd=z2sd)*w2[i2]*hu1)
    hl2[i2]=sum(dnorm(z2[i2]-ezl2,sd=z2sd)*w2[i2]*hl1)
   }
   k1=k2
   inf1=inf2
   m1=m2
   z1=z2
   hu1=hu2
   hl1=hl2
  }

  k2=k2+1
}

list(
     zbdy,
     pu,
     pl
     )
}

#-------------------------------------------------------------------------------

ersp2 = function(r,na,inf,cume)
{
 # 2-sided error spending test of H0: theta = 0 vs (theta < 0 or theta > 0).
 # with BINDING upper and lower boundaries.
 #
 # IN:
 #
 #   r             multiplier for Simpson's rule grid, e.g., r=16
 #   na            maximum number of analyses
 #   inf[1:na]     information levels
 #   cume[2;1:na]  cumulative probability under theta=0 of rejecting H0
 #                 in favour of theta > 0 
 #   cume[1;1:na]  cumulative probability under theta=0 of rejecting H0
 #                 in favour of theta < 0 
 #
 # OUT:
 #
 #   [1] zbdy[1;1:na]  lower boundary
 #       zbdy[2;1:na]  upper boundary
 #
 # NOTES:
 #
 # inf
 #
 # Information for theta at each analysis is given in inf[1:na],
 # so at analaysis k,
 #
 #  Z ~ N(theta*sqrt(inf(k)),inf(k)).
 #
 # Values in inf should be a strictly increasing positive sequence.
 # If this is not the case, an error will result.
 #
 # cume
 #
 # Cumulative type I error probabilities are specified in cume[2;1:na]
 # and cume[1;1:na]. For a symmetric test cume[1;1:na] = cume[2;1:na].
 # Each set of values should form an increasing positive sequence with
 # cume[1,na]+cume[2,na]=alpha.
 # Negative or decreasing values will result in an upper boundary at a
 # large positive value to approximate +infinity or a lower boundary at
 # a large negative value to approximate -infinity (spending ~ zero
 # type I error).
 #
 # Upper and lower boundaries are computed together.
 # For each k, the upper boundary, zbdy[2,1:na] is calculated first.
 # Then, the lower boundary, zbdy[1,1:na] is found. 
 #
 # Calculation uses a neat trick to use the one-sided test produced
 # by ersp1c as a two-sided test. This is done by setting the type I
 # error probability arising from crossing the lower boundary to be
 # the type II error probability of the one-sided test at the nominal
 # "alternative" theta=delta=0.

 delta=0
 out=ersp1c(r,na,inf,delta,cume)
 zbdy=out[[1]]
 pu=out[[2]]
 pl=out[[3]]

 list(
      zbdy,
      pu,
      pl
      )
}

#-------------------------------------------------------------------------------

find_rho_ersp1a=function(r,na,inf,alpha,beta,delta)
{
 # Finds the value rho such that a one-side error spending test with no
 # stopping to accept H0 (lower, futility boundary is set at Z = -50) has
 # type 1 and type 2 error rates precisely equal to alpha (at theta=0)
 # and beta (at theta=delta).

 fof_rho_1a=function(rho,par)
 {
  r=par[[1]]
  na=par[[2]]
  inf=par[[3]]
  alpha=par[[4]]
  beta=par[[5]]
  delta=par[[6]]
  infmax=inf[na]
  cume=array(0,c(2,na))
  cume[1,]=0
  cume[2,]=alpha*(inf/infmax)^rho
  out=ersp1a(r,na,inf,cume)
  zbdy=out[[1]]
  zbdy[1,na]=zbdy[2,na]
  out=gst1(r,na,inf,zbdy,delta)
  fof_rho=out[[2]]

  fof_rho
 }

 par=list(r,na,inf,alpha,beta,delta)
 rho0=1
 eps=10^(-5)
 rho.root=findroot(fof_rho_1a,1-beta,rho0,0.1,eps,par)

 # par=list(r,na,inf,alpha,beta,delta)
 # check=fof_rho_1a(rho.root,par)
 # print(c("Check: ",rho.root,check,beta))

 rho.root
}

#-------------------------------------------------------------------------------

find_rho_ersp1b=function(r,na,inf,alpha,beta,delta)
{
 # Finds the value rho such that a one-side error spending test
 # with an upper, efficacy boundary and a NON-BINDING lower,
 # futility boundary has type 1 and type 2 error rates precisely
 # equal to alpha (at theta=0) and beta (at theta=delta).

 fof_rho_1b=function(rho,par)
 {
  r=par[[1]]
  na=par[[2]]
  inf=par[[3]]
  alpha=par[[4]]
  beta=par[[5]]
  delta=par[[6]]
  infmax=inf[na]
  cume=array(0,c(2,na))
  cume[1,]=beta*(inf/infmax)^rho
  cume[2,]=alpha*(inf/infmax)^rho
  out=ersp1b(r,na,inf,delta,cume)
  zbdy=out[[1]]
  zbdy[1,na]=zbdy[2,na]
  out=gst1(r,na,inf,zbdy,delta)
  fof_rho=out[[2]]

  fof_rho
 }

 par=list(r,na,inf,alpha,beta,delta)
 rho0=1
 eps=10^(-5)
 rho.root=findroot(fof_rho_1b,1-beta,rho0,0.1,eps,par)

 # par=list(r,na,inf,alpha,beta,delta)
 # check=fof_rho_1b(rho.root,par)
 # print(c("Check: ",rho.root,check,beta))

 rho.root
}

#-------------------------------------------------------------------------------

find_rho_ersp1c=function(r,na,inf,alpha,beta,delta)
{
 # Finds the value rho such that a one-side error spending test
 # with an upper, efficacy boundary and a BINDING lower,
 # futility boundary has type 1 and type 2 error rates precisely
 # equal to alpha and beta

 fof_rho_1c=function(rho,par)
 {
  r=par[[1]]
  na=par[[2]]
  inf=par[[3]]
  alpha=par[[4]]
  beta=par[[5]]
  delta=par[[6]]
  infmax=inf[na]
  cume=array(0,c(2,na))
  cume[1,]=beta*(inf/infmax)^rho
  cume[2,]=alpha*(inf/infmax)^rho
  out=ersp1c(r,na,inf,delta,cume)
  zbdy=out[[1]]
  zbdy[1,na]=zbdy[2,na]
  out=gst1(r,na,inf,zbdy,delta)
  fof_rho=out[[2]]

  fof_rho
 }

 par=list(r,na,inf,alpha,beta,delta)
 rho0=1
 eps=10^(-5)
 rho.root=findroot(fof_rho_1c,1-beta,rho0,0.1,eps,par)

 # par=list(r,na,inf,alpha,beta,delta)
 # check=fof_rho_1c(rho.root,par)
 # print(c("Check: ",rho.root,check,beta))

 rho.root
}

#-------------------------------------------------------------------------------

find_rho_ersp2=function(r,na,inf,alpha,beta,delta)
{
 # Finds the value rho such that a two-side error spending test
 # with BINDING upper and lower boundaries has type 1 and type 2
 # error rates precisely equal to alpha and beta

 fof_rho_2=function(rho,par)
 {
  r=par[[1]]
  na=par[[2]]
  inf=par[[3]]
  alpha=par[[4]]
  beta=par[[5]]
  delta=par[[6]]
  infmax=inf[na]
  cume=array(0,c(2,na))
  cume[1,]=(alpha/2)*(inf/infmax)^rho
  cume[2,]=(alpha/2)*(inf/infmax)^rho
  out=ersp2(r,na,inf,cume)
  zbdy=out[[1]]
  out=gst1(r,na,inf,zbdy,delta)
  fof_rho=out[[2]]

  fof_rho
 }

 par=list(r,na,inf,alpha,beta,delta)
 rho0=1
 eps=10^(-5)
 rho.root=findroot(fof_rho_2,1-beta,rho0,0.1,eps,par)

 # par=list(r,na,inf,alpha,beta,delta)
 # check=fof_rho_2(rho.root,par)
 # print(c("Check: ",rho.root,check,beta))

 rho.root
}

#-------------------------------------------------------------------------------

find_rfac_ersp1a=function(r,na,alpha,beta,rho)
{
 # Finds the value rfac such that a one-side error spending test with
 # NO STOPPING TO ACCEPT H0 (lower, futility boundary is set at Z = -50)
 # and immediate response has type 1 and type 2 error rates precisely
 # equal to alpha and beta.
 #
 # Information level are equally spaced up to the maximum rfac*ifix.
 #
 # Calculations assume delta=1 but the result is invariant to this choice.

 fof_rfac_1a=function(rfac,par)
 {
  r=par[1]
  na=par[2]
  alpha=par[3]
  beta=par[4]
  rho=par[5]
  delta=1
  ifix=(qnorm(1-alpha)+qnorm(1-beta))^2
  infmax=rfac*ifix
  inf=infmax*c(1:na)/na
  cume=array(0,c(2,na))
  cume[1,]=0
  cume[2,]=alpha*(inf/infmax)^rho
  out=ersp1a(r,na,inf,cume)
  zbdy=out[[1]]
  out=gst1(r,na,inf,zbdy,delta)
  fof_rfac=out[[2]]

  fof_rfac
 }

 par=c(r,na,alpha,beta,rho)
 rfac0=1
 eps=10^-5
 rfac.root=findroot(fof_rfac_1a,1-beta,rfac0,0.1,eps,par)

 # par=c(r,na,alpha,beta,rho)
 # check=fof_rfac_1a(rfac.root,par)
 # print(c("Check: ",rfac.root,check,beta))

 rfac.root
}

#-------------------------------------------------------------------------------

find_rfac_ersp1b=function(r,na,alpha,beta,rho)
{
 # Finds the value rfac such that a one-side error spending test with
 # a NON-BINDING lower, futility boundary has type 1 and type 2 error
 # rates precisely equal to alpha and beta.
 #
 # Information level are equally spaced up to the maximum rfac*ifix.
 #
 # Calculations assume delta=1 but the result is invariant to this choice.

 fof_rfac_1b=function(rfac,par)
 {
  r=par[1]
  na=par[2]
  alpha=par[3]
  beta=par[4]
  rho=par[5]
  delta=1
  ifix=(qnorm(1-alpha)+qnorm(1-beta))^2
  infmax=rfac*ifix
  inf=infmax*c(1:na)/na
  cume=array(0,c(2,na))
  cume[1,]=beta*(inf/infmax)^rho
  cume[2,]=alpha*(inf/infmax)^rho
  out=ersp1b(r,na,inf,delta,cume)
  zbdy=out[[1]]
  out=gst1(r,na,inf,zbdy,delta)
  fof_rfac=out[[2]]

  fof_rfac
 }

 par=c(r,na,alpha,beta,rho)
 rfac0=1
 eps=10^-5
 rfac.root=findroot(fof_rfac_1b,1-beta,rfac0,0.1,eps,par)

 # par=c(r,na,alpha,beta,rho)
 # check=fof_rfac_1b(rfac.root,par)
 # print(c("Check: ",rfac.root,check,beta))

 rfac.root
}

#-------------------------------------------------------------------------------

find_rfac_ersp1c=function(r,na,alpha,beta,rho)
{
 # Finds the value rfac such that a one-side error spending test with
 # immediate response has type 1 and type 2 error rates precisely
 # equal to alpha and beta.
 #
 # Information level are equally spaced up to the maximum rfac*ifix.
 #
 # Calculations assume delta=1 but the result is invariant to this choice.

 fof_rfac_1c=function(rfac,par)
 {
  r=par[1]
  na=par[2]
  alpha=par[3]
  beta=par[4]
  rho=par[5]
  delta=1
  ifix=(qnorm(1-alpha)+qnorm(1-beta))^2
  infmax=rfac*ifix
  inf=infmax*c(1:na)/na
  cume=array(0,c(2,na))
  cume[1,]=beta*(inf/infmax)^rho
  cume[2,]=alpha*(inf/infmax)^rho
  out=ersp1c(r,na,inf,delta,cume)
  zbdy=out[[1]]
  out=gst1(r,na,inf,zbdy,delta)
  fof_rfac=out[[2]]

  fof_rfac
 }

 par=c(r,na,alpha,beta,rho)
 rfac0=1
 eps=10^-5
 rfac.root=findroot(fof_rfac_1c,1-beta,rfac0,0.1,eps,par)

 # par=c(r,na,alpha,beta,rho)
 # check=fof_rfac_1c(rfac.root,par)
 # print(c("Check: ",rfac.root,check,beta))

 rfac.root
}

#-------------------------------------------------------------------------------

find_rfac_ersp2=function(r,na,alpha,beta,rho)
{
 # Finds the value of the inflation factor rfac such that a two-side
 # error spending test with binding upper and lower boundaries has
 # two-sided type 1 error probability alpha at theta=0 and power 1-beta
 # at theta=delta.
 #
 # Information level are equally spaced up to the maximum rfac*ifix.
 #
 # Calculations assume delta=1 but the result is invariant to this choice.

 fof_rfac_2=function(rfac,par)
 {
  r=par[1]
  na=par[2]
  alpha=par[3]
  beta=par[4]
  rho=par[5]
  delta=1
  ifix=(qnorm(1-alpha/2)+qnorm(1-beta))^2
  infmax=rfac*ifix
  inf=infmax*c(1:na)/na
  cume=array(0,c(2,na))
  cume[1,]=(alpha/2)*(inf/infmax)^rho
  cume[2,]=(alpha/2)*(inf/infmax)^rho
  out=ersp2(r,na,inf,cume)
  zbdy=out[[1]]
  out=gst1(r,na,inf,zbdy,delta)
  fof_rfac=out[[2]]

  fof_rfac
 }

 par=c(r,na,alpha,beta,rho)
 rf20=1
 eps=10^-5
 rfac.root=findroot(fof_rfac_2,1-beta,rf20,0.1,eps,par)

 rfac.root
}

#-------------------------------------------------------------------------------