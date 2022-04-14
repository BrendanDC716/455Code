
function c = weights(z,x,m)

%%%%%^^^^^INPUTS HERE^^^^^^%%%%%



% Calculates FD weights. The parameters are:
%  z   location where approximations are to be accurate,
%  x   vector with x-coordinates for grid points,
%  m   highest derivative that we want to find weights for  (4 grid points,
%  3 m, 5 points 4 m etc
%  c   array size m+1,lentgh(x) containing (as output) in 
%      successive rows the weights for derivatives 0,1,...,m.
%test this method with a dummy function....polynomial%

%%%%Here we notice that x is the values for U, Z is the point we are trying
%%%%to approximate and m is the highest order derivative we want to
%%%%calculate, that runs from 0 to m, so for 5 grid points we want m=4 etc



%dx=.3142;
%h=.2566;


%%xi,1 position to the left of our interfact boundary



%%%left ui%%%%%

%x=[-h,dx-h,2*(dx)-h,3*(dx)-h,4*(dx)-h];


%%%%%right ui%%%%

%%%%contact dr li here!!!!!%%%%%
%x=[4*deltax+h,3*deltax+h,2*deltax+h,deltax+h,h];


%%xi-1, 2 positions of the left of the interface point


%%%%Left u i-1%%%
%x=[-deltax-h,deltax-h,2*(deltax)-h,3*(deltax)-h,4*(deltax)-h];


%%%xL=-1
%%x(ixl)=-1.256637
%%=h=-1--1.256637=deltz
%m=4;
%z=0;
    
n=length(x); c=zeros(m+1,n); c1=1; c4=x(1)-z; c(1,1)=1;
for i=2:n
   mn=min(i,m+1); c2=1; c5=c4; c4=x(i)-z;
   for j=1:i-1
      c3=x(i)-x(j);  c2=c2*c3;
      if j==i-1 
         c(2:mn,i)=c1*((1:mn-1)'.*c(1:mn-1,i-1)-c5*c(2:mn,i-1))/c2;
         c(1,i)=-c1*c5*c(1,i-1)/c2;
      end
      c(2:mn,j)=(c4*c(2:mn,j)-(1:mn-1)'.*c(1:mn-1,j))/c3;
      c(1,j)=c4*c(1,j)/c3;
   end
   c1=c2;
   
  
end


%c first set of weights
