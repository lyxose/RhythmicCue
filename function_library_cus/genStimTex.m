
function stiTexture = genStimTex(wptr, winRect, ut, bgContrast, tgContrast, tgCenter, GaborSF, GaborWidth, GaborOrient, bgWidth, seed)
    % to generate gray matrix of stimulus (Gabor in a round 1/f pinkNoise texture)
    % then make it the PTB texture
    % tgCenter (i.e. target center coor.) is in pixel unit, in image axis (right/down as positive, [0,0] for upleft corner)
    % all last 4 parameters are in degree! (GaborSF, GaborWidth, GaborOrient, bgWidth)

    scWidth  = winRect(3);  
    scHeight = winRect(4);
    bgCenter = [scWidth/2, scHeight/2];
    
    background = tPinkNoise(scWidth, seed, bgContrast); % full screen background texture
    background = background(1:scHeight,1:scWidth);
    lambda = ut.deg2pix(1/GaborSF);
    if lambda ~= 0
        Texture = grating(size(background), tgCenter, ...
                            1/lambda, GaborOrient, tgContrast);
        Texture = winOverlap(background, Texture, ut.deg2pix(GaborWidth), ...
                              tgCenter, 'cos'); 
    else
        fprintf('Spatial frequency was too large (%.4f) that grating cannot be generated at current screen resolution',GaborSF)
        Texture = background;
    end
    stimulus = winOverlap(zeros([scHeight,scWidth])+0.5, Texture, ...
                          ut.deg2pix(bgWidth), bgCenter, 'hard');
    stiTexture = Screen('MakeTexture', wptr, cat(3,stimulus,stimulus,stimulus)*255);
end     