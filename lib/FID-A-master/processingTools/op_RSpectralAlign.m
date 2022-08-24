function fids = op_RSpectralAlign(fids, water_flag, MRS_struct)

ii = MRS_struct.ii;
freqRange = MRS_struct.p.sw(ii)/MRS_struct.p.LarmorFreq(ii);
freq = (size(fids,1) + 1 - (1:size(fids,1))) / size(fids,1) * freqRange + 4.68 - freqRange/2;

n = 2;

D = zeros(size(fids,2)/n);
w = cell(1,n);
data = complex(zeros(size(fids,1),n));
dataLim = [];
time = (0:(MRS_struct.p.npoints(ii)-1))'/MRS_struct.p.sw(ii);
tMax = find(time <= 0.1,1,'last');

for jj = 1:n
        ind = find(MRS_struct.fids.ON_OFF == abs(jj-2));
    for kk = 1:size(fids,2)/n
        for ll = 1:size(fids,2)/n
            tmp = sum((real(fids(1:tMax,ind(kk))) - real(fids(1:tMax,ind(ll)))).^2) / 200;
            if tmp == 0
                D(kk,ll) = NaN;
            else
                D(kk,ll) = tmp;
            end
        end
    end
    d = nanmean(D);
    w{jj} = 1./d.^2;
    w{jj} = w{jj}/sum(w{jj});
    w{jj} = repmat(w{jj}, [size(fids,1) 1]);

        data(:,jj) = sum(w{jj} .* fids(:,ind),2);
end

clear tmp
[~, data] = FlattenData(data);



    ind = 2;

fids = PhaseCorrection(data, fids, freq, ind, MRS_struct);


data = WeightedAveraging(fids, data, n, water_flag, dataLim, w, MRS_struct);

[flatdata, data] = FlattenData(data);


    
switch MRS_struct.p.target{1}
    case {'GABAGlx','Lac','EtOH'}
        % Water
        freqLim = freq <= 4.68+0.22 & freq >= 4.68-0.22;
        [~,i] = max(abs(data(freqLim,:)));
        freq2 = freq(freqLim);
        maxFreq = freq2(i);
        for jj = 1:2
            freqLim(jj,:) = freq <= maxFreq(jj)+0.22 & freq >= maxFreq(jj)-0.22;
        end
end
freqLim = or(freqLim(1,:), freqLim(2,:));
f0 = (maxFreq(1) - maxFreq(2)) * MRS_struct.p.LarmorFreq(ii);
x0 = [f0 0];
    


lsqnonlinopts = optimoptions(@lsqnonlin);
lsqnonlinopts = optimoptions(lsqnonlinopts,'Algorithm','levenberg-marquardt','Display','off');
t = 0:(1/MRS_struct.p.sw(ii)):(size(fids,1)-1)*(1/MRS_struct.p.sw(ii));

   
a = max(max(max(flatdata)));
switch MRS_struct.p.target{1}
    case {'GABAGlx','Lac','EtOH'}
        fun = @(x) objFunc(flatdata./a, freqLim, t, x);
end
param = lsqnonlin(fun, x0, [], [], lsqnonlinopts);



ind = find(MRS_struct.fids.ON_OFF == 1);
for jj = 1:length(ind)
    fids(:,ind(jj)) = fids(:,ind(jj)) .* exp(1i*param(1)*2*pi*t') * exp(1i*pi/180*param(2));
end

if ishandle(44)
    close(44);
end
if ishandle(3)
    close(3);
end

end


function [flatdata, data] = FlattenData(data)

flatdata(:,1,:) = real(data);
flatdata(:,2,:) = imag(data);
data = real(fftshift(fft(data,[],1),1));

end


function fids = PhaseCorrection(data, fids, freq, ind, MRS_struct)

OFF = data(:,ind);
ii = MRS_struct.ii;

freqLim = freq <= 3.02+0.15 & freq >= 3.02-0.15;
[~,i] = max(abs(OFF(freqLim)));
freq2 = freq(freqLim);
maxFreq = freq2(i);
freqLim = freq <= maxFreq+0.58 & freq >= maxFreq-0.42;
OFF = OFF(freqLim);
Baseline = (OFF(1) + OFF(end))/2;
Width = 0.05;
Area = (max(OFF) - min(OFF)) * Width * 4;

x0 = [Area Width maxFreq 0 Baseline 0 1] .* [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];
ModelParamChoCr = FitChoCr(freq(freqLim), OFF, x0, MRS_struct.p.LarmorFreq(ii));

fids = fids * exp(1i*pi/180*ModelParamChoCr(4));

end


function data = WeightedAveraging(fids, data, n, water_flag, dataLim, w, MRS_struct)

for ii = 1:n

    ind = find(MRS_struct.fids.ON_OFF == abs(ii-2));

    if water_flag
        data(:,ii) = sum(w{ii}(:,1:dataLim) .* fids(:,ind(1:dataLim)),2);
    else
        data(:,ii) = sum(w{ii} .* fids(:,ind),2);
    end
end

end


function out = objFunc(in, freqLim, t, x)

f   = x(1);
phi = x(2);

y = complex(in(:,1,1), in(:,2,1));
y = y .* exp(1i*pi*(t'*f*2+phi/180));

a = real(fftshift(fft(y)));
b = real(fftshift(fft(complex(in(:,1,2), in(:,2,2)))));

DIFF = a - b;
out = DIFF(freqLim);

end



