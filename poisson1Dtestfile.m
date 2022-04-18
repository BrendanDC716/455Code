clear all
close all
format long

nbp = 2;        % number of boundary points, xl and xr 
njp = 5;        % number of jump conditions at each boundary pt
xl = -1;     % left  boundary pt
xr =  1;     % right boundary pt
a  = -pi;        % left  extended boundary pt a < xl
b  =  pi;        % right extended boundary pt xr < b
nx =  20;       % # of intervals

dx = (b-a)/nx;  % mesh size
x  = (a+dx):dx:(b-dx); % vecotr containing all interior mesh grids, size(x) = nx - 1

n1 = nx - 1;
n2 = nbp*njp;


%%%%%%%%%%%%%%%%%%
%%%%Construct A%%%
%%%%%%%%%%%%%%%%%%%


c1A=-(5/2); % main diagonal of A
c2A=4/3; % offset by 1 diagonal of A
c3A=-(1/12); % offset by 2 diagonal of A
c4A=1/(dx); % coefficient for entire A
% populate matrix A diagonals
A=diag(c1A*ones(1,nx-1)) + diag(c2A*ones(1,nx-2),1) + diag(c3A*ones(1,nx-3),2) + diag(c2A*ones(1,nx-2),-1) + diag(c3A*ones(1,nx-3),-2);
% apply antisymmetric property
A(1,1)=c1A-c3A; 
A(nx-1,nx-1)=c1A-c3A;
% multiply coefficient of A
A=c4A*A;

%%placeholder for dimensions of matrices
matrix_B = zeros(n1,n2); %%%HARD
matrix_C = zeros(n2,n1); %%%HARD
U = zeros(n1,1);
Q = zeros(n2,1);

ixl = 0; ixr = 0;  % ixl and ixr are indices of the grids on the left of the two boundary pts

% calculate the grid to the left of the two boundary pts
for i=1:nx-1
   if x(i) <= xl
      ixl = i;
   end

   if x(i) <= xr
      ixr = i;
   end
end

% need at least 2 grids between the original and extended boundaries and 4
% grids between xl and xr
if (ixl < 3) || (ixr - ixl < 4) || (ixr > nx-3)
    error('not enough mesh grids. refinement is requixred!')
end


xlp=x(ixl);
xrp=x(ixr);

%%%%Compute coeffecients weights for left Ui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hl=abs(xl-xlp); %h for left fictitious values

xl1=[-hl,dx-hl,2*(dx)-hl,3*(dx)-hl,4*(dx)-hl]; %first fict. value on left (ui)

uilorig=weights(0,xl1,4);


%dummyuilorig=[uilorig(1,1);uilorig(2,1);uilorig(3,1);uilorig(4,1);uilorig(5,1)];

%loop to solve for wbars
for n=1:5
    for i=1:5
    uildummy(i,n)=-uilorig(i,n)/uilorig(i,1);
    
    end
end

for n=1:4
    for i=1:5
uil(i,n)=uildummy(i,n+1);

    end
end
 


%%%%Compute weights for left Ui-1%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xl2=[-dx-hl,dx-hl,2*(dx)-hl,3*(dx)-hl,4*(dx)-hl]; %2nd fict value on left (ui-1)
ui2lorig=weights(0,xl2,4);


%%%loop to solve for wbars
for n=1:5
    for i=1:5
    ui2ldummy(i,n)=-ui2lorig(i,n)/ui2lorig(i,1);
    
    end
end

for n=1:4
    for i=1:5
ui2l(i,n)=ui2ldummy(i,n+1);

    end
end

h2=abs(xr-xrp); %h for right ficititios values
%%%%%%%%


%xl1=[-hl,dx-hl,2*(dx)-hl,3*(dx)-hl,4*(dx)-hl]

%%z=dx-h

%%[-3dx,-2dx,-dx,0,dx]
xr1=[-3*dx,-2*(dx),-(dx),0,dx];
zstar=h2;
uirorig=weights(zstar,xr1,4);
%%%loop to solve for wbars
for n=1:5
    for i=1:5
    uirdummy(i,n)=-uirorig(i,n)/uirorig(i,5);
    
    end
end

for n=1:4
    for i=1:5
uir(i,n)=uirdummy(i,n);

    end
end

%%second fict value for right side
xr2=[-3*dx,-2*(dx),-(dx),0,2*dx];
uir2orig=weights(zstar,xr2,4);
%%%loop to solve for wbars
for n=1:5
    for i=1:5
    uir2dummy(i,n)=-uir2orig(i,n)/uir2orig(i,5);
    
    end
end

for n=1:4
    for i=1:5
uir2(i,n)=uir2dummy(i,n);

    end
end

UIL1hat=zeros(1,5);
%%%program UI hats and UI-1 hats from paper, equation 1,2,3 etc

for i=1:4
    
        
    UIL1hat(1,i)=uil(1,i);
end
    
    UIL1hat(1,5)=1/uilorig(1,1); %%%add 1/Wi



UIL2hat=zeros(1,5);

for i=1:4
    
        
    UIL2hat(1,i)=ui2l(1,i);
end
    
    UIL2hat(1,5)=1/ui2lorig(1,1); %%%add 1/Wi-1
    
    
%%Right interface Ui
UIR1hat=zeros(1,5);

for i=1:4
    
        
    UIR1hat(1,i)=uir(1,i);
end
    
    UIR1hat(1,5)=1/uirorig(1,5); %%%add 1/Wi
    
    %%right interface ui
UIR2hat=zeros(1,5);
%%%program UI hats and UI-1 hats from paper, equation 1,2,3 etc

for i=1:4
    
        
    UIR2hat(1,i)=uir2(1,i);
end
    
    UIR2hat(1,5)=1/uir2orig(1,5); %%%add 1/Wi-1
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%attempt to construct C%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  
    
    
    %%construct left[xi-1,xi,xi+1,xi+2,xi+3] 5by5 W  Interface at Xl
   WXL=[-dx,0,dx,2*(dx),3*(dx)];
   WXL5by5= weights(hl,WXL,4); 
   
   

   %%construct right[xi-3,xi-2,xi-1,xi,xi+1] interface at Xr
   WXR=[-2*(dx),-(dx),0,dx,2*dx];
   WXR5by5 = weights(h2,WXR,4);
    
    
   %%Construct each sub matrix for Xl and Xr in matrix C
    for k=1:5
        
       C1(:,k) = (WXL5by5(:,1)*UIL2hat(k)+WXL5by5(:,2)*UIL1hat(k)+WXL5by5(:,3));
       C2(:,k) = (WXR5by5(:,1)*UIR2hat(k)+WXR5by5(:,2)*UIR1hat(k)+WXR5by5(:,3));
    end
    
   for j=1:3
       C1(:,j) = (C1(:,j)+WXL5by5(:,j));
        C2(:,j) = (C2(:,j)+WXR5by5(:,j));
        
   end
    
    
    
    
    %%%Construct Phi
    phi=zeros(10,1);
   for i=1:5
       phi(i,1)=WXL5by5(i,5)*xl;
       phi(i+5,1)=WXR5by5(i,5)*xr;
     
   end
 
   
   

      






