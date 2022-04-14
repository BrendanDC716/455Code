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



%%%%Construct A%%%


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
%%for loop to calculate ui-1l weights
h2=abs(xr-xrp); %h for right ficititios values
xr1=[-3*dx-h2,-2*(dx)-h2,-(dx)-h2,-h2,dx-h2];
uirorig=weights(0,xr1,4);
%%%loop to solve for wbars
for n=1:5
    for i=1:5
    uirdummy(i,n)=-uirorig(i,n)/uirorig(i,1);
    
    end
end

for n=1:4
    for i=1:5
uir(i,n)=uirdummy(i,n+1);

    end
end

%%second fict value for right side
xr2=[-3*dx-h2,-2*(dx)-h2,-(dx)-h2,-h2,2*dx-h2];
uir2orig=weights(0,xr2,4);
%%%loop to solve for wbars
for n=1:5
    for i=1:5
    uir2dummy(i,n)=-uir2orig(i,n)/uir2orig(i,1);
    
    end
end

for n=1:4
    for i=1:5
uir2(i,n)=uir2dummy(i,n+1);

    end
end

UIL1hat=zeros(1,5);
%%%program UI hats and UI-1 hats from paper, equation 1,2,3 etc

for i=1:4
    
        
    UIL1hat(1,i)=uil(1,i);
end
    
    UIL1hat(1,5)=1/uilorig(1,1); %%%add 1/Wi

%%%check corrected taylor expansio 


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
    
    UIR1hat(1,5)=1/uirorig(1,1); %%%add 1/Wi
    
    %%right interface ui
UIR2hat=zeros(1,5);
%%%program UI hats and UI-1 hats from paper, equation 1,2,3 etc

for i=1:4
    
        
    UIR2hat(1,i)=uir2(1,i);
end
    
    UIR2hat(1,5)=1/uir2orig(1,1); %%%add 1/Wi-1
    
    
    %%deriv jumps ck
    
    %%%%k=0:4
    
    %%%%%attempt to construct C%%%%%
    
    for k=1:4
        
       C1(:,k) = -(uilorig(:,1)*UIL2hat(k)+uilorig(:,2)*UIL1hat(k)+uilorig(:,3)+uilorig(:,4)+uilorig(:,5));
       C2(:,k) = -(uirorig(:,1)*UIR2hat(k)+uirorig(:,2)*UIR1hat(k)+uirorig(:,3)+uirorig(:,4)+uirorig(:,5));
    end

      






