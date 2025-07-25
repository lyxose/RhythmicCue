
function stream = genStream(ITI, cSOA, cFreq, tSOA, tFreq, maxRT, stiD, sampRate, noiseAmp, cueAmp, tgAmp, ramp, noiseSeed)

    if nargin < 13 || isempty(noiseSeed)
        noiseSeed = randi([0, intmax('uint32')], 1, 'uint32'); 
    end
    originalRNG = rng;  % Preserve global random number generator state 
    rng(noiseSeed);  
    stream = noiseAmp.*(rand(1,round(sampRate*(ITI+sum(cSOA)+tSOA+maxRT))).*2-1);
    rng(originalRNG);  % Maintain global random number consistency

    cueSti = MakeBeep(cFreq,stiD,sampRate).*cueAmp;
    tgSti = MakeBeep(tFreq,stiD,sampRate).*tgAmp;
    sti_len = length(cueSti);
    % add ramp to each beep
    ramp_env = ones(1,round(sti_len));   
    ramp_len = round(sampRate*ramp);
    ramp_env(1:ramp_len)=linspace(0,1,ramp_len);
    ramp_env(length(ramp_env)-ramp_len+1:end)=linspace(1,0,sampRate*ramp);
    cueSti = cueSti.*ramp_env;
    tgSti = tgSti.*ramp_env;
    for c = round(sampRate.*(cumsum(cSOA)+ITI))
        stream(c+1:c+sti_len) = stream(c+1:c+sti_len)+cueSti;
    end
    t = c+round(sampRate*tSOA);
    stream(t+1:t+sti_len) = stream(t+1:t+sti_len)+tgSti;
end
