function A = LDPC(m,n)

     q = sqrt(n);
     
     P = eye(q);
     P = P([q 1:q-1],:);    %permutation with shift 1.
     l=m/q;
     A=zeros(m,n);
     % constructing as mentioned
     for i=1:l
         for j=1:q
             A([(i-1)*q+1:i*q],[(j-1)*q+1:j*q]) = P^((i-1)*(j-1));
         end
     end
end