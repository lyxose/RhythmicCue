
function stiTexture = genStimTex(wptr, winRect, ut, tgContrast, tgCenter, GaborSF, GaborWidth, GaborOrient)
    % to generate gray matrix of stimulus (Gabor texture with soft edge)
    % then make it the PTB texture
    % tgCenter (i.e. target center coor.) is in pixel unit, in image axis (right/down as positive, [0,0] for upleft corner)
    % all last 3 parameters are in degree! (GaborSF, GaborWidth, GaborOrient)
    scWidth  = winRect(3);  
    scHeight = winRect(4);

    bgCenter = winRect(3:4)/2;
    lambda = ut.deg2pix(1/GaborSF);

    Texture = grating(winRect(3:4), tgCenter, 1/lambda, ...
                         GaborOrient, tgContrast);
    stimulus = winOverlap(zeros([scHeight,scWidth])+0.5, Texture, ...
                          ut.deg2pix(GaborWidth), bgCenter, 'cos');
    stiTexture = Screen('MakeTexture', wptr, cat(3,stimulus,stimulus,stimulus)*255);
end     