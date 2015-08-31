# Draft R code for BUGS model to run hierachical Bayesian model
#  model { phi[y*k] p[y*k] Psi[y*k] }  for k locations scenario: k=1, k=2, k = 1...n observable, k= n+1 unobserved/dead

# data file requirements: K locations, N[y], T, Y, z[y,i,t], w[y,i,t], c[y,i]
# initial values file requirements
# z[y,i,t],
# mu.p[k]
# tau.p[k]
# mu.phi[k]
# tau.phi[k]                                     
# mu.Psi[k]                                      
# tau.Psi[k]                                     
# etaBphi[k]
# etaBp[k]
# etaBphi[k]
# etap[y,k]
# etaphi[y,k]
#etaPsi[y,k]


{
	mu.p[1] ~ dnorm(0,0.001) 		
	tau.p[1] ~ dgamma(0.001,0.001) 		
	mu.p[2] ~ dnorm(0,0.001) 
	tau.p[2] ~ dgamma(0.001,0.001) 

	mu.phi[1] ~ dnorm(0,0.001) 
	tau.phi[1] ~ dgamma(0.001,0.001) 
	mu.phi[2] ~ dnorm(0,0.001) 
	tau.phi[2] ~ dgamma(0.001,0.001) 

	mu.Psi[1] ~ dnorm(0,0.001) 
	tau.Psi[1] ~ dgamma(0.001,0.001) 
	mu.Psi[2] ~ dnorm(0,0.001) 
	tau.Psi[2] ~ dgamma(0.001,0.001)

mean.p[1] <- exp(mu.p[1]) / (1+exp(mu.p[1]))
mean.p[2] <- exp(mu.p[2]) / (1+exp(mu.p[2]))
mean.phi[1] <- exp(mu.phi[1]) / (1+exp(mu.phi[1]))
mean.phi[2] <- exp(mu.phi[2]) / (1+exp(mu.phi[2]))
mean.Psi[1] <- exp(mu.Psi[1]) / (1+exp(mu.Psi[1]))
mean.Psi[2] <- exp(mu.Psi[2]) / (1+exp(mu.Psi[2]))

	logit(Bpd.phi[1]) <- etaBphi[1] 
	etaBphi[1] ~ dnorm(mu.phi[1], tau.phi[1])
	logit(Bpd.phi[2]) <- etaBphi[2] 
	etaBphi[2] ~ dnorm(mu.phi[2], tau.phi[2])

	logit(Bpd.p[1]) <- etaBp[1] 
	etaBp[1] ~ dnorm(mu.p[1], tau.p[1])
	logit(Bpd.p[2]) <- etaBp[2] 
	etaBp[2] ~ dnorm(mu.p[2], tau.p[2])

	logit(Bpd.Psi[1]) <- etaBPsi[1] 
	etaBPsi[1] ~ dnorm(mu.Psi[1], tau.Psi[1])
	logit(Bpd.Psi[2]) <- etaBPsi[2] 
	etaBPsi[2] ~ dnorm(mu.Psi[2], tau.Psi[2])
        
for (y in 1:Y){
	logit(p[y,1]) <- etap[y,1] 
	etap[y,1] ~ dnorm(mu.p[1],tau.p[1]) 
	logit(p[y,2]) <- etap[y,2] 
	etap[y,2] ~ dnorm(mu.p[2],tau.p[2]) 
	p[y,K+1] <- 0		
	logit(phi[y,1]) <- etaphi[y,1] 
	etaphi[y,1] ~ dnorm(mu.phi[1],tau.phi[1]) 
	logit(phi[y,2]) <- etaphi[y,2] 
	etaphi[y,2] ~ dnorm(mu.phi[2],tau.phi[2]) 

	logit(Psi[y,1,1]) <- etaPsi[y,1] 
	etaPsi[y,1] ~ dnorm(mu.Psi[1],tau.Psi[1]) 
	logit(Psi[y,2,2]) <- etaPsi[y,2] 
	etaPsi[y,2] ~ dnorm(mu.Psi[2],tau.Psi[2]) 

	Psi[y,1,2] <- 1 - Psi[y,1,1]
	Psi[y,2,1] <- 1 - Psi[y,2,2]

	q[y, K+1, K+1] <- 1					
	for (k in 1:K){
		q[y, K+1, k] <- 0				
		q[y, k, K+1] <- 1 - phi[y, k]			

		for (j in 1:K){
			q[y,k,j] <- phi[y,k] * Psi[y,k,j] 			
			}
		}			
	for (i in 1:N[y]){				
		for (t in (c[y,i]+1):T){				

		       z[y,i,t] ~ dcat(q[y,z[y,i,t-1],1:(K+1)])
 
		       w[y,i,t] ~ dbern(p[y,z[y,i,t]])	      
		       }
		}
	}

