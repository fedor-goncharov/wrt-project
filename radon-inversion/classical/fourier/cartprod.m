function P = cartprod( varargin )
     if ( nargin < 1 )
  	     error('usage: cartprod(A,B,C, ... )')
     end
      
     %% kludge allows to use input matrix arguments
     if ( rows(varargin{1}) > 1 && columns(varargin{1}) > 1 )
 	     A = varargin{1};
	    for itr = 1:columns(A)
		  varargin{itr} = A(:,itr);
 	    end
     end
      
     itr = 1; P = [];
     while ( itr <= length(varargin) )
	  B = varargin{itr};
          B = reshape(B,[length(B),1]);      
	  rP = rows(P); lB=length(B);
	  if ( rP == 0 ) 
		P = B; itr = itr + 1;
		continue;
	  end
	  Q = repmat(P,[lB,1]);	  
          R = reshape(repmat(B,[1,rP]).',[lB*rP,1]);
          P = [Q,R];
          itr = itr + 1;
     end     
     return
end