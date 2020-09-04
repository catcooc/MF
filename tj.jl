using LinearAlgebra
const t=-1 #-1
const tt=-0.3*t
const ttt=0.2*t #-0.16
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


function get_order(U,ek,x,y,fz=1,T=0.)
	
	if fz == 1
	jg=sum([conj(U[x,i])*U[y,i]*(theta(ek[i],T)) for i in 1:length(ek)])
	
	else

	jg=sum([U[x,i]*conj(U[y,i])*(1-theta(ek[i],T)) for i in 1:length(ek)])
 	end
 	return jg

end







function hz(kx,ky
 	,Chi=0.
  	,m=0.
    ,dop=0.
    ,mu=0.,
    T=0.
	)
    x=dop
    Q=pi 
    xs=1
	ek1=-(2*x*t +  xs*J * Chi)*(cos(kx)+cos(ky))-4*x*tt*cos(kx)*cos(ky)-2*x*ttt*(cos(2*kx)+cos(2*ky))-mu
    ekQ=-(2*x*t + xs *J * Chi)*(cos(kx+Q)+cos(ky+Q))-4*x*tt*cos(kx+Q)*cos(ky+Q)-2*x*ttt*(cos(2*kx+2*Q)+cos(2*ky+2*Q))-mu
	
  	MM=2*J*m
	h=[ek1 -MM 0 0;-MM ekQ 0 0 ; 0 0 ek1 MM ; 0 0 MM ekQ]
    ek,U=eigen(h) 
   	Chi=((get_order(U,ek,1,1,1,T)+get_order(U,ek,3,3,1,T))*(cos(kx)+cos(ky))
		+(get_order(U,ek,2,2,1,T)+get_order(U,ek,4,4,1,T))*(cos(kx+Q)+cos(ky+Q)))/(2*N)
   
    m=(get_order(U,ek,1,2,1,T)-get_order(U,ek,3,4,1,T)
		+get_order(U,ek,2,1,1,T)-get_order(U,ek,4,3,1,T))/(2*N)

	n=(get_order(U,ek,1,1,1,T)+get_order(U,ek,2,2,1,T)
		+get_order(U,ek,3,3,1,T)+get_order(U,ek,4,4,1,T))/N

	return  ([Chi m n])
 end


function find_mu(mu,n,Chi,m,dop,T=0.,j=100,eta=10^-6)
	mumax=10
	mumin=-10

	x=sum(map(x->hz(x[1],x[2],Chi,m,dop,mu,T),MK))[end]
 for i in 1:2000
	if abs(n-x) > eta
   		(mumin,mumax) = n-x>0 ?  (mu,mumax) : (mumin,mu)
   		 # if i%2 == 0 
   			# println("mumax:",mumax ," mumin:",mumin, " mu:",mu)

   		 # 	println("n:",n," n_get:",x)
   		 # end
   		mu=(mumin+mumax)/2
   		x=sum(map(x->hz(x[1],x[2],Chi,m,dop,mu,T),MK))[end]
   	else 
   		if (j<20&&j>5) || j%50 ==0
   			println("over in :" ,i)
   		end
   		break

   	end

  end
  	if mumax == mumin
   			println("nm " ,j )
   	end
   return  mu
end

# Chi,m=rand(2).+0.0001
# dop=0.1
# T=0.001
# Chi,m,nz=sum(map(x->hz(x[1],x[2],Chi,m,dop,0,T),MK))
# println(find_mu(0,1-dop,nz,Chi,m,dop,T))

 function MF(dop=0.
    ,mu=0.
    ,T=0.)

wc=zeros(3)
jg=zeros(3)
n=1-dop
sl=0
Chi,m=rand(2).+0.001
mu=find_mu(mu,n,Chi,m,dop,T,0,0.001*eta)
Chi,m,nz=sum(map(x->hz(x[1],x[2],Chi,m,dop,mu,T),MK))
println(Chi,"  "," ",m)
for i in 1:5000

	
		mu=find_mu(mu,n,Chi,m,dop,T,i,0.001*eta)

        	
	jg=sum(map(x->hz(x[1],x[2],Chi,m,dop,mu,T),MK))
   

	wc=hcat((jg[1]-Chi)/(Chi+im*eta),(m-jg[2])/(m+im*eta)
		,(nz-jg[3])/(nz+im*eta))
	
	wc=map(abs,wc)
	
 	Chi,m,nz=jg
	# Chi,m=(map(x->x*0.7,transpose(jg[1:2]))+map(x->0.3*x,[Chi  m]))
	# nz=jg[3]
	#Chi,m,nz=jg
	 # if abs(n-nz) > eta
  #  		if  (n-nz) > 0 
  #  			mumin=mu
  #  			#println("dd")
  #  		else
  #  			mumax=mu
  #  			#println("mm")
  #  		end

   		#(mumin,mumax)=( n-nz > 0 ?  (mu,mumax) :  (mumin,mu)   )# (mumin,mu) : (mu,mumax)  
   		# mu=(mumin+mumax)/2
   	 # end
   	 wnc=abs(n-nz)

	

    # bc=wnc*bc1

	if (sum(wc.>eta) == 0  &&   abs(n-nz) <=eta )
		sl=1
		println("ddd",i)
		break

	end

	if (i<20&&i>5) || i%50 == 0
		println("i ",i)
		println("-----------------------------------------------------------------------------")
		println("t:",t," tt:",tt," ttt:",ttt," J:",J," T:",T," eta:",eta)
		println(wc)
		println(jg)
		println("wnc:",wnc)
		println("nn:",n , " n:" ,nz)
		println("mu:",mu)
		println("-----------------------------------------------------------------------------")

	end
end

return (Chi,m,nz,sl,mu)
end

#println(start())
#println(MF(0.01))


function caculate(T=0.)
	println("start")
fCx=open(string("Chi2z",".txt"),"w")
fM=open(string("M2z",".txt"),"w")
fn=open(string("n2z",".txt"),"w")
fsl=open(string("sl2z",".txt"),"w")
mu=0
close(fCx)
close(fsl)
close(fn)
close(fM)
for i in  0:10
	fCx=open(string("Chi2z",".txt"),"a")
fM=open(string("M2z",".txt"),"a")
fn=open(string("n2z",".txt"),"a")
fsl=open(string("sl2z",".txt"),"a")
println(i) 
dop=(i*0.02)
jg=MF(dop,mu,T)
println(jg)
write(fCx,string(jg[1]))
write(fCx,"\n")
write(fM,string(jg[2]))
write(fM,"\n")
write(fsl,string(jg[4]))
write(fsl,"\n")
write(fn,string(jg[3]))
write(fn,"\n")
println("t:",t," tt:",tt," ttt:",ttt," J:",J," T:",T)
println("------------------------------")
mu=jg[5]
close(fCx)
close(fsl)
close(fn)
close(fM)
end



end

caculate(0.001)

# Chi,m=rand(2).+0.0001
# for i in -1:0.1:1
# jg=sum(map(x->hz(x[1],x[2],Chi,m,i,0),MK))
# println(jg)

# end

# jg=sum(map(x->hz(x[1],x[2],Chi,m,0,-0.25,0),MK))
# println(jg)

# jg=sum(map(x->hz(x[1],x[2],Chi,m,0,-0.28,0),MK))
# println(jg)
# for i in 1:10
# T=0.001
# Chi,m=rand(2).+0.001
# dop=0.02
# mu=0
# n=1-dop
# nz=sum(map(x->hz(x[1],x[2],Chi,m,dop,mu,T),MK))[end]
# mu=find_mu(mu,nz,n,Chi,m,dop,T,0,0.001*eta)
# jg=sum(map(x->hz(x[1],x[2],Chi,m,dop,mu,T),MK))
# Chi,m,nz=jg
# println(nz)
# end
