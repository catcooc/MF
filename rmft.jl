using LinearAlgebra
#using GR
const t=0.2
const tt=-0.05
const tc=0.28
const J=0.1
const K=1
const L=50
const eta=10^-5
const N=L*L
#const KK=[((i*2*pi)/L,(j*2*pi)/L)  for i in 1:L for j in 1:L]
const KK=[((i*2*pi)/L,(j*2*pi)/L)  for i in -L/2+1:L/2 for j in -L/2+1:L/2]
nrand=1
#const KKK=[map(x->(x[1]+rand()*(2*pi/L),x[2]+rand()*(2*pi/L)),KK) for i in 1:nrand]
const KKK= nrand==1 ? [ KK for i in 1:nrand] : [map(x->(x[1]+rand()*(2*pi/L),x[2]+rand()*(2*pi/L)),KK) for i in 1:nrand]
const nc=0.1

function theta(x,T=0)
	if T >0 
		fek=1/(1+exp(x/T))
	else 
		fek=x>0 ? 0 : 1
	end
return  fek
end


function get_order(U,ek,x,y,fz=1,T=0)
	
	if fz == 1
	result=sum([conj(U[x,i])*U[y,i]*theta(ek[i],T) for i in 1:length(ek)])
	
	else

	result=sum([U[x,i]*conj(U[y,i])*(1-theta(ek[i],T)) for i in 1:length(ek)])
 	end
 	return result

end


function hz(kx,ky,
	B=0.
	,D=0.
	,chix=0+0*im
	,chiy=0+0*im
    ,chix_=0+0*im
    ,chiy_=0+0*im
	,deltax=0+0*im
	,deltay=0+0*im
    ,deltax_=0+0*im
    ,deltay_=0+0*im
	,dop=0.
	,mu1=0.,mu2=0.,T=0.)

    # B=0.
    # D=0.
    # deltax=0+0*im
    # deltay=0+0*im
  

    nh=dop+nc
    gt=nh/(1+nh)
    gJ=4/(1+nh)^2
    gk=2/(1+nh)
      xs=3/8
    Chi= -(t*gt + xs*J*gJ*chix)*cos(kx)-(t*gt + xs*J*gJ*chix_)*cos(-kx)-(t*gt + xs*J*gJ*chiy)*cos(ky)-(t*gt +xs*J*gJ*chiy_)*cos(-ky)-tt*gt*2*(cos(kx+ky)+cos(kx-ky))+mu1
    ek_c=-tc*2*(cos(kx)+cos(ky))+mu2
    delta=-xs * J *gJ*((deltax+deltax_)*cos(kx)+(deltay_+deltay)*cos(ky))
    KD=-3/4 * gk * K * D /(2^.5)
    KB=-3/4 * gk * K * B /(2^.5)
    
    

    h= [Chi  KD  conj(delta)  conj(KB);
       conj(KD)  ek_c  conj(KB)  0 ;
       delta  KB  -Chi  -conj(KD) ;
       KB  0  -KD  -ek_c]
   
    ek,U=eigen(Hermitian(h))
    # ek2,U2=eigen(h)
    #   println(ek)
    #   println(ek2)
    chi1=get_order(U,ek,1,1,1,T)
    chi2=get_order(U,ek,3,3,0,T)
    delta1=get_order(U,ek,1,3,1,T)
    delta2=get_order(U,ek,3,1,0,T)
    B=(get_order(U,ek,1,4,1,T)-get_order(U,ek,3,2,0,T))/(N*2^.5)
    D=(get_order(U,ek,2,1,1,T)+get_order(U,ek,4,3,0,T))/(N*2^.5)
    deltax=(delta1*exp(im*kx)-delta2*exp(-im*kx))/N
    deltay=(delta1*exp(im*ky)-delta2*exp(-im*ky))/N
    deltax_=(delta1*exp(-im*kx)-delta2*exp(im*kx))/N
    deltay_=(delta1*exp(-im*ky)-delta2*exp(im*ky))/N
    nc_get=(get_order(U,ek,2,2,1,T)+get_order(U,ek,4,4,0,T))/N
    nd=((chi1+chi2))/N
    chix=((chi1*exp(-im*kx)+chi2*exp(im*kx)))/N
    chiy=((chi1*exp(-im*ky)+chi2*exp(im*ky)))/N
    chix_=((chi1*exp(im*kx)+chi2*exp(-im*kx)))/N
    chiy_=((chi1*exp(im*ky)+chi2*exp(-im*ky)))/N
   
    # ek2,U2=eigen(Hermitian(h))
    # println(Hermitian(h))
    # println(h)
    # println(inv(U)*h*(U))

    return ([B D chix chiy chix_ chiy_ deltax  deltay deltax_  deltay_ nc_get  nd])
end
#     dop=0
#     mu1=0
#     mu2=0
#     T=0
#     B,D=rand(2).+0.001
#     chix,chiy,deltax,deltay =rand(ComplexF64,4).+0.001
#   result=sum(map(xx->sum(map(x->hz(x[1],x[2],B,D,chix,chiy,deltax,deltay,dop,mu1,mu2,T),xx)),KKK))/nrand
# println(result)



function find_mu1(mu1,mu2,B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_,n,dop=0.,T=0.,eta=10^-6)
    mumax=5
    mumin=-5
    mu=mu1
    nd_get=abs(start(0,mu,mu2,dop,T,B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_)[end])
    n_get=1-nd_get
    errors=n-n_get
    convergence=0
    for i in 1:3000
        if abs(errors) >eta

            #mumin,mumax= (-1*errors) > 0 ? (mu,mumax) : (mumin,mu)
             mumin,mumax= errors > 0 ? (mu,mumax) : (mumin,mu)
        else 
            
            #println("over in : ",i)
            convergence=i
            break

        end
       
         mu=(mumin+mumax)/2
        nd_get=abs(start(0,mu,mu2,dop,T,B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_)[end])
        n_get=1-nd_get
        errors=n-n_get
    end
    if mumax == mumin
        println("nm")
        convergence="nm"
    end
    return [mu convergence]
end



function find_mu2(mu1,mu2,B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_,n,dop=0.,T=0.,eta=10^-6)
    mumax=5
    mumin=-5
    mu=mu2
    n_get=abs(start(0,mu1,mu,dop,T,B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_)[end-1])

    errors=n-n_get
    convergence=0
    for i in 1:3000
       
        if abs(errors) >eta
            #mumin,mumax= errors > 0 ? (mu,mumax) : (mumin,mu)
            mumin,mumax= (-1*errors) > 0 ? (mu,mumax) : (mumin,mu)
           
        else 
            #println("over in : ",i)
            convergence=i
            break
        end


        mu=(mumin+mumax)/2
        n_get=abs(start(0,mu1,mu,dop,T,B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_)[end-1])
        errors=n-n_get
      
        if mumax == mumin
            convergence=string("nm: ",i)
        end
    end
    return [mu convergence]
end

function start(r=1,mu1=0.
    ,mu2=0.,dop=0.,T=0.,B=0.
    ,D=0.,chix=0+0*im
    ,chiy=0+0*im
    ,chix_=0+0*im
    ,chiy_=0+0*im
    ,deltax=0+0*im
    ,deltay=0+0*im
    ,deltax_=0+0*im
    ,deltay_=0+0*im
    )
    nh=nc+dop


    if r == 0 
      
       

        result=sum(map(xx->sum(map(x->hz(x[1],x[2],B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_,dop,mu1,mu2,T),xx)),KKK))/nrand

        
       
      
    else
        println("start")
        
        for i in 1:500
        r2=find_mu2(mu1,mu2,B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_,nc,dop,T)
        
        r1=find_mu1(mu1,r2[1],B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_,nh,dop,T)
        if abs(mu2-r2[1])<eta && abs(mu1-r1[1])< eta
            if i%100==0
                println("ddfdd")
            end
            break
        else
            if i%100==0
            println("dddd")
            println(mu2 ," ",r2[1]," ", abs(mu2-r2[1]))
            end
           mu2=r2[1]
           mu1=r1[1]
        end
        end
        result=sum(map(xx->sum(map(x->hz(x[1],x[2],B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_,dop,mu1,mu2,T),xx)),KKK))/nrand
    end
return result
end
#start(1,0.,0.,0.,0.)

function MF(B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_,mu1=0.,mu2=0.,dop=0.,T=0.)
	convergence=0
	nd=1-nc-dop
    nh=nc+dop
	B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_,nc_get,nd_get=start(1,mu1,mu2,dop,T,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_)
	error_nc=abs(nc - nc_get)
    error_nd=abs(nd - nd_get)
	for i  in 1:1500
            
         for jj in 1:500
            global r1
            global r2
            r2=find_mu2(mu1,mu2,B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_,nc,dop,T)
            r1=find_mu1(mu1,r2[1],B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_,nh,dop,T)
            if abs(mu2-r2[1])<eta && abs(mu1-r1[1])< eta
                if jj%100==0
                println("ddfdd")
                end
                break
            else
                if i%100==0
                println("dddd")
                end
                mu2=r2[1]
                mu1=r1[1]
            end
            end
            result=start(0,mu1,mu2,dop,T,B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_)
		    errors=hcat((result[1]-B)/(B+eta),(D-result[2])/(D+eta),(chix-result[3])/(chix+eta)
		          ,(chiy-result[4])/(chiy+eta),(chix_-result[5])/(chix_+eta)
                  ,(chiy_-result[6])/(chiy_+eta),(deltax-result[7])/(deltax+eta),(deltay-result[8])/(deltay+eta)
                  ,(deltax_-result[9])/(deltax_+eta),(deltay_-result[10])/(deltay_+eta)
		         ,(nc_get-result[11])/(nc_get+eta),(nd_get-result[12])/(nd_get+eta))
	
		errors=map(abs,errors)
		error_nc=abs(nc - nc_get)
    	error_nd=abs(nd - nd_get)
    	B,D,chix,chiy,chix_,chiy_,deltax,deltay,deltax_,deltay_,nc_get,nd_get=result
        nc_get=abs(nc_get)
        nd_get=abs(nd_get)
	#if (sum([(10>=i>=4 ? errors[i]>0.001 : errors[i] >eta), [i for i in 1:12]) == 0  &&   error_nc <=eta && error_nd <=eta )
        if (sum(errors.>eta) == 0  &&   error_nc <=eta && error_nd <=eta )
			convergence=1
			println("over in ",i)
            println(result)
            println(nd_get)
            println(nc_get)
			break

		end
		if (i<=20)|| i%20==0 || i%21==0
        #if  i >10
			println("tc:",tc," i:",i," T:",T," eta:",eta," nrand:",nrand," L:",L)
			println("-----------------------------------------------------------------------------")
			println("errors")
            println(errors)
            println("B D nc nd ")
            println(map(abs,result[1:2])," ",map(abs,result[11:end]))
            println("chix chiy chix_ chiy_")
            println(result[3:6])
            println(result[7:10])
            println("ds: ",abs(nh/(1+nh)*(deltax+deltay)/2)," dd:",abs(nh/(1+nh)*(deltax-deltay)/2))
			println("error_nc:",error_nc ,"  error_nd:",error_nd )
			println("nc:", nc,"  nc_get:" ,nc_get,"   nh:",nh,"  nh_get:",1 - nd_get)
			println(mu1 ,"  ",mu2)
            println(r1[2],"  ", r2[2])
			#println(mu11 ,"  ",mu22)
			println("-----------------------------------------------------------------------------")

		end
	end
    convergence=abs(convergence)
   return([B D chix chiy chix_ chiy_ deltax  deltay deltax_  deltay_  nc_get  1-nd_get mu1 mu2 convergence])


end
	

#println(MF(0.03263248429615925,0.980440315290581))



function caculata(T=0)

	 fnc=open("Nc.txt","w")
     fnh=open("Nh.txt","w")
     fchix=open("chix.txt","w")
     fchiy=open("chiy.txt","w")
     fdeltax=open("deltax.txt","w")
     fdeltay=open("deltay.txt","w")
     fB=open("B.txt","w")
     fD=open("D.txt","w")
     fmu1=open("mu1.txt","w")
     fmu2=open("mu2.txt","w")
     fconvergence=open("convergence.txt","w")
     close(fnc)
     close(fnh)
     close(fchix)
     close(fchiy)
     close(fdeltay)
     close(fdeltax)
     close(fB)
     close(fD)
     close(fmu1)
     close(fmu2)
     close(fconvergence)
	dop=[i for i in -0.1:0.1:0.3]
	nn=length(dop)
	mu1=0
	mu2=0
	B=zeros(ComplexF64,nn)
    D=zeros(ComplexF64,nn)
    chix=zeros(ComplexF64,nn)
    chiy=zeros(ComplexF64,nn)
    deltax=zeros(ComplexF64,nn)
    deltay=zeros(ComplexF64,nn)
    Nc=zeros(nn)
    Nh=zeros(nn)
    sl=zeros(nn)
    Bs,Ds=rand(2).+0.001
        chixs,chiys,deltaxs,deltays =rand(ComplexF64,4).+0.001
        chixs_,chiys_,deltaxs_,deltays_ =rand(ComplexF64,4).+0.001
 	for i in 1:nn
 		println("i:",i)
 		fnc=open("Nc.txt","a")
     	fnh=open("Nh.txt","a")
     	fchix=open("chix.txt","a")
     	fchiy=open("chiy.txt","a")
     	fdeltax=open("deltax.txt","a")
     	fdeltay=open("deltay.txt","a")
        fmu1=open("mu1.txt","a")
        fmu2=open("mu2.txt","a")
     	fB=open("B.txt","a")
     	fD=open("D.txt","a")
        fconvergence=open("convergence.txt","a")

 		result=MF(Bs,Ds,chixs,chiys,chixs_,chiys_,deltaxs,deltays,deltaxs_,deltays_,mu1,mu2,dop[i],T)
        println(result)
          Bs,Ds,chixs,chiys,chixs_,chiys_,deltaxs,deltays,deltaxs_,deltays_=result[1:10].+0.001
 		B[i],D[i],chix[i],chiy[i]=result[1:4]
        deltax[i],deltay[i]=result[7:8]
 		Nc[i],Nh[i],mu1,mu2,sl[i]=result[11:end]
 		write(fnc,string(Nc[i]))
 		write(fnc,"\n")
 		write(fnh,string(Nh[i]))
 		write(fnh,"\n")
 		write(fconvergence,string(sl[i]))
 		write(fconvergence,"\n")
 		write(fB,string(B[i]))
 		write(fB,"\n")
 		write(fD,string(D[i]))
 		write(fD,"\n")
 		write(fchix,string(chix[i]))
 		write(fchix,"\n")
 		write(fchiy,string(chiy[i]))
 		write(fchiy,"\n")
 		write(fdeltax,string(deltax[i]))
 		write(fdeltax,"\n")
 		write(fdeltay,string(deltay[i]))
 		write(fdeltay,"\n")
        write(fmu1,string(mu1))
        write(fmu1,"\n")
        write(fmu2,string(mu2))
        write(fmu2,"\n")
		close(fnc)
     	close(fnh)
     	close(fchix)
     	close(fchiy)
     	close(fdeltay)
     	close(fdeltax)
     	close(fB)
     	close(fD)
     	close(fconvergence)
        close(fmu1)
        close(fmu2)
        mu1=0
        mu2=0
     end
     println(sl)
 end

caculata(0)


