% in this the particles are placed in a square area
% it turns out there are some errors in this code, the values don't match
% with the papar. In the next file we'll try to set the particles 
% in a circular pattern from the beginning
%%






tic


Np=30;                        %particles
Nt=20;                           %tests
Nsims=1000;     %MC simulations
Runmax=100000;
Potcenter=[0 0];
arrangement=[3 9 15 20 30 40 50 60 70 80 90 100 ];

%units for conversions
kB=1.3806503*10^(-23);%boltzmann constant (J/K)
hbar=6.62606957*10^(-34);  % (J*s)
q=1.60217646*10^(-19);%charge of electron
epsilon=1;%
m=9.10938188*10^(-31);%mass of electron
omega=23*10^9/hbar; 
alpha=(m*omega^2)/2;


Rnul=(q^2/epsilon)^(1/3)*alpha^(-1/3);
Enul=(q^2/epsilon)^(2/3)*alpha^(1/3);
Tnul=(q^2/epsilon)^(2/3)*alpha^(1/3)/kB;


Temp=0.003; % moet je nog met die speciale eenheden  linken aan een echte temp.

Stepmax=0.25;%10^30*Rnul;
%%



randpos=zeros(Np,2);
Hpotential=zeros(1,Np);
Hinteraction=zeros(Np,Np-1);

Epotential=zeros(1,Np);
Einteraction=zeros(Np,Np-1);

H=zeros(Nsims,Np);
Etot=zeros(1,Nsims);
x=zeros(1,Nsims);

X=zeros(Np,Nsims);
Y=zeros(Np,Nsims);

    

%% In this part we'll try to arrange the particles in a certain configuration

shellsum=0;

shells=0;

while shellsum < Np
    
   shellsum=shellsum+arrangement(shells+1);    
    shells=shells+1;
end


%this loop will see that there aren't too many empty spots in the last
%shell, if so well put more particles in the lower shells
if shellsum-Np> arrangement(shells)/2 
    
    shells=shells-1;
    
    
end

    

% the shells parameter tells us how many shells we will fill for the first config. 

% distance between shells, this is all quite arbitrarry
% the first shell will be placed at half the intershell distance
% from the center of the potential

intershell=1/(shells);







%%

inshellparticle=0;
shelliter=1;
angleseparation=zeros(1,shells);
radialdistance=zeros(1,shells);
particleamount=0;
preX=zeros(1,Np);
preY=zeros(1,Np);
for iter=1:1:1 
    
    for i=1:1:shells  %in this loop there is a subloop so that all the remaining
                      % particles are placed in the last shell 
                      % (no opening of new low filled shell)
        
       if i<shells
        angleseparation(i)=2*pi/arrangement(i);
        radialdistance(i)=intershell/2+intershell*(i-1);
        particleamount=particleamount+arrangement(i);
       else
           angleseparation(i)=2*pi/(Np-particleamount);
        radialdistance(i)=intershell/2+intershell*(i-1);
           
        
       end
       
       
       
    end
    
    
    
   for j=1:1:Np  
      
       if shelliter~=shells && inshellparticle-arrangement(shelliter)>=0;
            
              shelliter=shelliter+1;
             
              inshellparticle=0;
       end
       
           
       
       
       
       preX(j)=radialdistance(shelliter)*cos(angleseparation(shelliter)*inshellparticle);
       preY(j)=radialdistance(shelliter)*sin(angleseparation(shelliter)*inshellparticle);
       
       inshellparticle=inshellparticle+1;
      
       
       X(j,iter)=preX(j)+0.08*(rand()*2-1);
       Y(j,iter)=preY(j)+0.08*(rand()*2-1);
   end
    
   
   % if E(i-1)>E(i): accept new positions and iter=iter+1
    
   for k=1:1:Np
       
       Hpotential(k)=(X(k,iter)-Potcenter(1))^2+(Y(k,iter)-Potcenter(2))^2;
      
       
       for l=1:1:Np-k
           
           Hinteraction(k,l)=1/sqrt((X(k,iter)-X(k+l,iter))^2+(Y(k,iter)-Y(k+l,iter))^2);
           
       end
       
   end
   
    H(iter)=sum(Hpotential)+sum(sum(Hinteraction));
    Etot(iter)=H(iter);
    x(iter)=iter;  % dummy variable for easy plotting
    iter=iter+1; %#ok<FXSET>
    x(iter)=iter;
end
    

%%

totalruns=0;
faulty=0;
correctly=0;
wrong=0;


 %stap moet kleiner worden naar mate je meer naar het dal gaat

while iter<=Nsims
        
    
    if Np>300
         disp('Number of particles is too much')
        break
    end
    
    
    
    
        totalruns=totalruns+1; 
       
        
        for j=1:1:Np  % for each simulation we generate random positions
       % this loop can be implemented otherwise (only once 
       % using rand in main loop and then selecting the numbers)
      
       %here you need to implement temperature
       
       
       randpos(j,1)=X(j,iter-1)+rand()*Stepmax-Stepmax/2;  %x position based on previous x position
       randpos(j,2)=Y(j,iter-1)+rand()*Stepmax-Stepmax/2;  %y position based on previous y position
       
        
    
        for k=1:1:Np
       
             Hpotential(k)=(randpos(k,1)-Potcenter(1))^2+(randpos(k,2)-Potcenter(2))^2;
      
       
                for l=1:1:Np-k
           
                    Hinteraction(k,l)=1/sqrt((randpos(k,1)-randpos(k+l,1))^2+(randpos(k,2)-randpos(k+l,2))^2);
           
                end
       
        end
  
        H(iter,j)=sum(Hpotential)+sum(sum(Hinteraction));
    
        
        
        
        if j~=1  % j not equal to 1
        
        
        if H(iter,j-1)>H(iter,j)
            correctly=correctly+1;
          
                 % for each simulation we generate random positions
                                 % this loop can be implemented otherwise (only once 
                                 % using rand in main loop and then selecting the numbers)
      
                    X(j,iter)=randpos(j,1);
                    Y(j,iter)=randpos(j,2);
                    
                
                
                
        else if rand()<=exp(-(H(iter,j)-H(iter,j-1))*Enul/(Temp*Tnul*kB))
                faulty=faulty+1;
                
                
          % zet stap groter
          % dummy variable for easy plotting
                  % for each simulation we generate random positions
                                 % using rand in main loop and then selecting the numbers)
      
                    X(j,iter)=randpos(j,1);
                    Y(j,iter)=randpos(j,2);
                    
                    
                    
            else
                X(j,iter)=X(j,iter-1);
                Y(j,iter)=Y(j,iter-1);
                wrong=wrong+1;
                
                        
  
                         H(iter,j)=H(iter-1,j);
                    
                
            end
            
            
        end 
        
        
        
        
        
        
        
        else if j==1
            
            if H(iter-1,Np)>H(iter,j)
            correctly=correctly+1;
          
                 % for each simulation we generate random positions
                                 % this loop can be implemented otherwise (only once 
                                 % using rand in main loop and then selecting the numbers)
      
                    X(j,iter)=randpos(j,1);
                    Y(j,iter)=randpos(j,2);
                    
                
                
                
               
                
            else if rand()<=exp(-(H(iter,j)-H(iter-1,Np))*Enul/(Temp*Tnul*kB))
                faulty=faulty+1;
                
                
          % zet stap groter
          % dummy variable for easy plotting
                  % for each simulation we generate random positions
                                 % using rand in main loop and then selecting the numbers)
      
                    X(j,iter)=randpos(j,1);
                    Y(j,iter)=randpos(j,2);
                    
                
                    
                    
                    
                else
                    X(j,iter)=X(j,iter-1);
                    Y(j,iter)=Y(j,iter-1);
                    wrong=wrong+1;
                
                        
  
                         H(iter,j)=H(iter-1,j);
                    
                
                end
               
                
            
            end
            
            
            end 
            
            
            
            
        end 
        end
        
        
        
        
        
        
        for k=1:1:Np
       
             Epotential(k)=(X(k,iter)-Potcenter(1))^2+(Y(k,iter)-Potcenter(2))^2;
      
       
                for l=1:1:Np-k
           
                    Einteraction(k,l)=1/sqrt((X(k,iter)-X(k+l,iter))^2+(Y(k,iter)-Y(k+l,iter))^2);
           
                end
       
        end
  
        Etot(iter)=sum(Epotential)+sum(sum(Einteraction));
        
        x(iter)=iter;
        iter=iter+1;
        
        
        
        
        throwawayratio=faulty/wrong;
        
        if throwawayratio>=0.5
            Stepmax=Stepmax*(1005/1000);
            faulty=0;
            wrong=0;
        else
        Stepmax=Stepmax*(995/1000);
            faulty=0;
            wrong=0;
        end

        
        
        if totalruns>Runmax
            disp('Number of iterations is more than Runmax, process interupted')
            break
        end
        
end



EoverN=Etot(Nsims)/Np;



toc
%%
figure(1)
plot(x,Etot)
xlabel('simulation step')
ylabel('H (energy)')


% plot the random walks

figure(2)
clf
hold on
for plotiter = 1:1:Np
    
    line(X(plotiter,:),Y(plotiter,:),'color','k');
    
end
hold off

%%
figure(3)

plot(X(:,1),Y(:,1),'r.')
title('starting positions')

%%
figure(4)
plot(X(:,Nsims),Y(:,Nsims),'k.')
title('ending positions')


figure(5)
voronoi(X(:,Nsims),Y(:,Nsims),'k.')
title('voronoi diagram')










% to do: zie dat de eerste deeltjes op een gestructureerde manier random
% worden gekozen (pseudorandom dus: verdeel ze wat over de kwadranten).


