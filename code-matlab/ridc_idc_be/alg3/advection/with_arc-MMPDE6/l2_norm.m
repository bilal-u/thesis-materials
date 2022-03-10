
function nr=l2_norm(w,n)
     
         accum = 0.0;
         for i=1:n
            accum = accum+ w(i)*w(i);
         end
         nr=sqrt(accum);
         
end