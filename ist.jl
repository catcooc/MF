using LinearAlgebra
const t=-1 #-1
const tt=0.32
const ttt=-0.16 #-0.16
const J=0.3
const eta=10^-4
const L=200
const N=L*L*2
 KK=[((i*2*pi)/L,(j*2*pi)/L)  for i in -L/2:L/2-1 for j in -L/2:L/2-1]
const MK=map(x->((x[1]+x[2])/(2^.5),(x[2]-x[1])/(2^.5)),KK)

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
	jg=sum([conj(U[x,i])*U[y,i]*(theta(ek[i],T)) for i in 1:length(ek)])
	
	else

	jg=sum([U[x,i]*conj(U[y,i])*(1-theta(ek[i],T)) for i in 1:length(ek)])
 	end
 	return jg

end



function hz(kx,ky
 	,Chi=0
 	,de=0
  	,m=0
    ,dop=0
    ,mu=0,
    T=0
	)
    x=dop
    Q=pi

	ek1=-(2*x*t + 3/4 *J * Chi)*(cos(kx)+cos(ky))-4*x*tt*cos(kx)*cos(ky)-2*x*ttt*(cos(2*kx)+cos(2*ky))-mu
    ekQ=-(2*x*t + 3/4 *J * Chi)*(cos(kx+Q)+cos(ky+Q))-4*x*tt*cos(kx+Q)*cos(ky+Q)-2*x*ttt*(cos(2*kx+2*Q)+cos(2*ky+2*Q))-mu
	Deltak=-3/4 *J *de*(cos(kx)-cos(ky))
	DeltakQ=-3/4 *J *de*(cos(kx+Q)-cos(ky+Q))

  	MM=-4*J*m
	h=[ek1 Deltak MM 0 ; Deltak -ek1 0 MM ; MM 0 ekQ DeltakQ ; 0 MM DeltakQ -ekQ]

    ek,U=eigen(h) 
   Chi=((get_order(U,ek,1,1,1,T)+get_order(U,ek,2,2,0,T))*(cos(kx)+cos(ky))+(get_order(U,ek,3,3,1,T)+get_order(U,ek,4,4,0,T))*(cos(kx+Q)+cos(ky+Q)))/(2*N)
   de=((get_order(U,ek,1,2,1,T)-get_order(U,ek,2,1,0,T))*(cos(kx)-cos(ky))+(get_order(U,ek,3,4,1,T)-get_order(U,ek,4,3,0,T))*(cos(kx+Q)-cos(ky+Q)))/(2*N)

   m=(get_order(U,ek,1,3,1,T)-get_order(U,ek,4,2,0,T)+get_order(U,ek,3,1,1,T)-get_order(U,ek,2,4,0,T))/(2*N)

	n=(get_order(U,ek,1,1,1,T)+get_order(U,ek,2,2,0,T)+get_order(U,ek,3,3,1,T)+get_order(U,ek,4,4,0,T))/N

	return  ([Chi de m n])
 end

function find_mu(mu,n,Chi,de,m,dop,T=0,j=0,eta=10^-6)
	mumax=10
	mumin=-10
	x=sum(map(x->hz(x[1],x[2],Chi,de,m,dop,mu,T),MK))[end]
 for i in 1:2000
	if abs(n-x) > eta
   		(mumin,mumax)= n-x>0  ?  (mu,mumax) : (mumin,mu)
   		mu=(mumin+mumax)/2
   		x=sum(map(x->hz(x[1],x[2],Chi,de,m,dop,mu,T),MK))[end]
   	else 
   		if  (j<20&&j>5) || j%50 ==0
   			println("over in :", i," x:",x)
   		end
   		break
   	end
  end
  		if  (j<20&&j>5) || j%50 ==0
  			if mumax==mumin
   				println("nm")
   			end
   		end
   return  mu
end

 function start(dop=0,mu=0,T=0)
 	println(start)
	Chi,de,m=rand(3).+0.0001

	jg=sum(map(x->hz(x[1],x[2],Chi,de,m,dop,mu,T),MK))
return jg
end


function MF(Chi,de,m,dop=0.
    ,mu=0.
    ,T=0.)

wc=zeros(4)
jg=zeros(8)
n=1-dop
sl=0

mu=find_mu(mu,n,Chi,de,m,dop,T,0,eta)
nz=sum(map(x->hz(x[1],x[2],Chi,de,m,dop,mu,T),MK))[end]
println(Chi,"  ",de," ",m)

for i in 1:8000
	mu=find_mu(mu,n,Chi,de,m,dop,T,i,eta)
	jg=sum(map(x->hz(x[1],x[2],Chi,de,m,dop,mu,T),MK))
	#jg=sum(map(x->hzcs(x[1],x[2],B,Dz,Chix,Chiy,dex,dey,nc,nh,mu1,mu2),KK))
	
   
    	
        	

   

	wc=hcat((jg[1]-Chi)/(Chi+im*eta),(de-jg[2])/(de+im*eta),(m-jg[3])/(m+im*eta)
		,(nz-jg[4])/(nz+im*eta))
	
	wc=map(abs,wc)
	
    Chi,de,m,nz=jg
	 # Chi,de,m=(map(x->x*0.7,jg[1:3]) + transpose(map(x->0.3*x, [Chi de m])))
	
	 # nz=jg[4] 

   	wnc=abs(n-nz)

	

	if (sum(wc.>eta) == 0  &&   abs(n-nz) <=eta )
		sl=1
		println("ddd",i)
		break

	end
	if i%50 == 0
		println(" ")
		println("-----------------------------------------------------------------------------")
		println("t:",t," tt:",tt," ttt:",ttt," J:",J," T:",T,"eta:",eta)
		println(wc)
		println(jg)
		println("wnc:",wnc)

		println("nn:",n , " n:" ,nz)
		println("mu:",mu)

		println("-----------------------------------------------------------------------------")

	end
end

return (Chi,de,m,nz,sl,mu)
end

#println(start())
#println(MF(0.01))


function caculate(T=0)
	println("start")
fCx=open(string("Chi3",".txt"),"w")
fM=open(string("M3",".txt"),"w")
fDe=open(string("de3",".txt"),"w")
fn=open(string("n3",".txt"),"w")
fsl=open(string("sl3",".txt"),"w")
mu=0
Chi,de,m=rand(3).+0.001
for i in  0:10
println(i) 
dop=(i*0.03)
jg=MF(Chi,de,m,dop,mu,T)

println(jg)
write(fCx,string(jg[1]))
write(fCx,"\n")
write(fDe,string(jg[2]))
write(fDe,"\n")
write(fM,string(jg[3]))
write(fM,"\n")
write(fsl,string(jg[5]))
write(fsl,"\n")
write(fn,string(jg[4]))
write(fn,"\n")
mu=jg[6]
Chi,de,m=jg
Chi+=0.001
de+=0.001
m+=0.001
end
close(fDe)
close(fCx)
close(fsl)
close(fn)
close(fM)
println("t:",t,"tt:",tt,"ttt:",ttt,"J:",J)

end

caculate(eta)
