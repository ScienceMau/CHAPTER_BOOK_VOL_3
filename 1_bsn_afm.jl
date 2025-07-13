using DynamicalSystems, SharedArrays, MAT, Statistics, OrdinaryDiffEq, ProgressMeter


#####################################################
##  function for mathematical modeling
######################################################
 @inline @inbounds function MathModel(u,par,t)

   X1  = 0.05;
  eta  = 0.01;
	k  = 0.5;
   f0  = 0.2;
lambda = 0.01;
   b0  = 0.5;
omega0 = 1.0;
   a0  = 0.2;
	D  = 2.0;
	n  = 3;
   k1  = 0.25;
	a  = 0.5;
	b  = 0.55;
	A  = 1.0;
gamma  = 0.45;
    p  = 0.2;
    delta  = par[1];
    phi = par[2];

 du1 = u[2];
 du2 = -2*eta*u[2]+0.5*u[1]*(1+2*delta*u[1]-u[1]^2)+X1*u[3]+f0*cos(omega0*t-a0*cos(b0*omega0*t))-a*k1*u[1]-(1-a)*D*k1*u[4]-p*sin(phi);
 du3 = -lambda*u[3]-k*u[2];
 du4 = D^(-1)*(A*u[2]-b*abs(u[2])*(abs(u[4])^(n-1)*u[4]-gamma*u[2]*abs(u[4])^n));

	return SVector{4}(du1, du2, du3, du4);
end


############################################
# Function write text in screen
###########################################
function text_screen(text)
	
	for i in text
	print(i)		
	sleep(0.05)
	end
end


 function scatter_attractors(attractors, i1, i2)
    for k in keys(attractors)
        x, y = columns(attractors[k])
        save_1(k, x, y, i1, i2)
    end
end

 function save_1(k, k2, k3, i1, i2)
   
    matwrite("$(k)_attractors_$(i1)_$(i2).mat", Dict(
        "x" => collect(k2),
        "y" => collect(k3),
        ))
end

 function save_basins(k2,i1,i2)
 
   matwrite("basins_$(i1)_$(i2).mat", Dict(
        "bsn" => collect(k2),
        ))
end


text_screen("Bem vindo ao Basins_attraction_V0.01.jl \n")




delta = [-0.15];
phi = [0.0;5.0;10.0]; 

for i=1:1:length(phi)
	for j = 1:1:length(delta)
		xg = yg = range(-2.0, 2.0; length = 800)
		zg = range(0.0, 1.0; length = 3)
		wg = range(-1.0, 1.0; length = 3)
		u0= [0.1; 0.0; 0.0; 0.0]
		ds   = ContinuousDynamicalSystem(MathModel, u0, [delta[j]; phi[i]])
		diffeq = ( Δt = 0.01, reltol = 1e-4, abstol =1e-4)
		mapper =  AttractorsViaRecurrences(ds, (xg, yg, zg, wg); diffeq)
		basins, attractors = basins_of_attraction(mapper; show_progress = true )
 
		scatter_attractors(attractors,phi[i],delta[j])
		save_basins(basins,phi[i],delta[j])
	end
end



text_screen("Developed by:Mauricio A. Ribeiro\n")
text_screen("Email:mau.ap.ribeiro@gmail.com\n")
