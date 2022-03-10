

T=1.2;

a=0.02;

p=4;

Store_p_val=zeros(p+1,p);

Store_er_val=zeros(p+2,p);

for i=1:p
    
    if i==1 
        
      [Or1,Er1]=idc1(a,T);
      
      Store_p_val(:,1)=Or1(1:p+1,1);
      
      Store_er_val(:,1)=Er1(1:p+2,1);
      
      
    else         
       
      s=matrix(i);
    
     [Or,Er]=idc(a,i,s,T);
     
     Store_p_val(:,i)=Or(1:p+1,1);
      
     Store_er_val(:,i)=Er(1:p+2,1);
    
    end
    
end

format short e
    
Er=Store_er_val;

digits(4)

format short g

Or=zeros(p+2,p);

Or(1,:)=0;

Or(2:p+2,:)=Store_p_val(1:p+1,:);

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

Vector= A;

disp(Vector)

%sprintf('%e  %f %e  %f %e  %f %e  %f\n',Vector)

 digits(4)
 
 result=latex(sym(vpa(A)))
 
 