function X= tdma1(aW,aP,aE,su);
n = length(aW);
X = zeros(n,1);
  A1=-aE;  
  C=su;              
   % Gauss elimination    
        A1(1)=aE(1)/aP(1);   
        C(1)=su(1)/aP(1);
        
for i=2:n
     A1(i)=aE(i)/(aP(i)-(aW(i)*A1(i-1)));
     C(i)=(su(i)+((aW(i)*C(i-1))))/(aP(i)-(aW(i)*A1(i-1)));
            
end
X(n)=C(n);


  %Backward substitution
for i=n-1:-1:1
    
    X(i)=C(i)+(A1(i)*X(i+1));
        
end
      


