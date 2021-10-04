function theta =  omp(A,b,eps)
      r = b;
      indices = [];
      theta = zeros(size(A,2),1);
      normA = A./sqrt(sum(A.*A, 1));
      while norm(r)^2 > eps
          [~, index] =  max(abs(r'*normA));
          indices = [indices; index];
          A1 = A(:,indices);
          theta(indices)= A1\b;
          r = b - A1*theta(indices);
      end 
end