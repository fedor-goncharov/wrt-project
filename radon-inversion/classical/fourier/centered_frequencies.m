function freqs_cent = centered_frequencies(n, dfreq)
  
  freqs_cent = zeros(1, n);
  for k = 0 : n-1
    if (k < n/2)
      freqs_cent(k+1) = k;
    elseif (k > n/2)
      freqs_cent(k+1) = k - n;
    endif
  endfor
  
  if (mod(n,2) == 0)
    freqs_cent(n/2 + 1) = -n/2;
  endif
    
  freqs_cent = fftshift(freqs_cent * dfreq);
  
endfunction