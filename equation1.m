function u = equation1(n,L,vit,vi)
kN=ones(n,1);
kS=ones(n,1);
vitA=0;
vitB=0;
for i=1:n
    if i==1
    kN(i)=(vi(i)+(vit(i)+vit(i+1))./2);
    kS(i)=(vi(i)+(vitA));
    
    elseif i==n
     kN(i)=(vi(i)+vitB);
    kS(i)=(vi(i)+(vit(i)+vit(i-1))./2);
    
    else
     kN(i)=(vi(i)+(vit(i)+vit(i+1))./2);
    kS(i)=(vi(i)+(vit(i)+vit(i-1))./2);
    end
end
%k(1:n)= vi(1:n)+vit(1:n); %coefficeint of visocus stress+ reynolds stress
del_y=(L)./(n); %stencile size
% Finite Volume
for i=1:n
    if i==1 %At Bottom Boundary
    ae(i)=kN(i)/del_y;
    aw(i)=0;
    sp(i)=-2*vi(i)./del_y;
    ap(i)=ae(i)+aw(i)-sp(i);
    su(i)=0+del_y;
     elseif i==n %At Upper Boundary
    ae(i)=0;
    aw(i)=kS(i)./del_y;
    sp(i)=-2*vi(i)./del_y;
    ap(i)=ae(i)+aw(i)-sp(i);
    su(i)=0+del_y;
     else
    ae(i)=kN(i)./del_y;
    aw(i)=kS(i)./del_y;
    ap(i)=ae(i)+aw(i);
    su(i)=1*del_y;
    end
end

%For solving TDMA
u=tdma1(aw,ap,ae,su);
end