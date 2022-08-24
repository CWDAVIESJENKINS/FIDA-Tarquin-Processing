function [Spectrum,sv] = op_ApplySVD(IN , FFTres, N )


if ~exist('N','var')
    N = 20;
end
keyboard

Spectrum = zeros([IN.sz,N]);
Sv = zeros(IN.sz(IN.dims.averages),N);

for NAV=1:10%IN.sz(IN.dims.averages)
    [U,S,V] = svd(hankel(IN.fids)); %U*S*V'  = H(FID)
    Sv(NAV,:) = diag(S); %singular values

    for NSing=1:N
        Ns = NAV;
        FID = U(:,Ns)*diag(Sv(Ns))*V(:,Ns)';
        Spectrum(:,NAV,NSing) = fftshift(fft(FID(:,1),FFTres));
    end

end

end