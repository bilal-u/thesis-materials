

T=0.12;

p=4;

% time step

a=0.01; 

% solutions from ode15s 

Store_p_val=zeros(p+2,p);

Store_er_val=zeros(p+2,p);

for i=1:p
    
    if i==1 
       
      % results from idc1  
             
      [Or1,Er1]=idc1(a,T);
      
      
      Store_p_val(:,1)=Or1(1:p+2,1);
      
      Store_er_val(:,1)=Er1(1:p+2,1);
      
      
    else  
        
      % integartion matrix calling  
       
      s=matrix(i);
      
      % results from       
         
     [Or,Er]=idc(a,i,s,T);
     
     
     
     Store_p_val(:,i)=Or(1:p+2,1);
      
      Store_er_val(:,i)=Er(1:p+2,1);
    
    end
    
end

% creating latex tables

format short e
    
Er=Store_er_val;

digits(4)

format short g

Or=Store_p_val;

digits(4)

A=zeros(p+2,2*p);

A(:,1)=Er(:,1);
A(:,2)=Or(:,1);

A(:,3)=Er(:,2);
A(:,4)=Or(:,2);

A(:,5)=Er(:,3);
A(:,6)=Or(:,3);

A(:,7)=Er(:,4);
A(:,8)=Or(:,4);

Vector= A

digits(4)
 
result=latex(sym(vpa(A)))

