function [oper, texture] = showInstruc_Inp(wptr, instName, instFolder, Key1, Key2, backKey, minSlideT)
    % Show instruction images in sequence based on filename order.
    % 
    % showInstruc_Inp
    % This function displays a sequence of instruction images from a specified 
    % folder, allowing the user to navigate through them using keyboard inputs. 
    % The function supports navigation to the next or previous image and enforces 
    % a minimal reading duration for each slide.
    %
    % Parameters:
    % wptr        : Window pointer (returned by Screen('OpenWindow')).
    % instName    : Partial name to match instruction image files. Only files 
    %               containing this string in their names will be loaded.
    % instFolder  : Folder path containing the instruction image files.
    % Key1        : Key name for navigating to the next image (e.g., 'space').
    % backKey     : Key name for navigating to the previous image (e.g., 'backspace').
    % minSlideT   : Minimal reading duration (in seconds) for each slide. Default 
    %               is 0, except for the last slide when there are multiple slides.
    %
    % Returns:
    % oper        : Operation result indicating user action:
    %               1  - User pressed the first key.
    %               2  - User pressed the second key.
    %              -1  - User pressed the back key.
    % texture     : Texture object of the currently displayed image, which can 
    %               be used for further processing or display.

    % Read the list of image file names from the specified folder
    afiles = dir(fullfile(instFolder, 'EmptyMatchedStruct'));
    for ftype={'*.png', '*.jpg', '*.jpeg'}
        files = dir(fullfile(instFolder, ftype{1})); % Assuming JPEG images
        afiles = [afiles; files];
    end
    fileNames = {afiles.name};
    
    % Match the file path based on the given name
    matchedFiles = fileNames(contains(fileNames, instName));
    
    % Check if there is at least one match
    if isempty(matchedFiles)
        error('No file matched found for the instruction image.');
    end
    
    % Get the first matched image file path and draw it
    imgPointer = 1;
    imgPath = fullfile(instFolder, matchedFiles{imgPointer});
    texture = drawCentImg(wptr, imgPath, 'fit');
    t0 = Screen('Flip', wptr);
    
    % Ensure Key1, Key2, and backKey are cell arrays
    if ischar(Key1)
        Key1 = {Key1};
    end
    if ischar(Key2)
        Key2 = {Key2};
    end
    if ischar(backKey)
        backKey = {backKey};
    end

    % Wait for a key press
    while 1
        oper = 0;
        while nargin > 5 && GetSecs - t0 < minSlideT
            WaitSecs(0.1);  % wait until reach the minimal time of each slide
        end
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            % Check if Key1 is a cell array or a single string
            if any(cellfun(@(k) keyCode(KbName(k)), Key1))
                oper = 1; % first key pressed
            end
            % Check if Key2 is a cell array or a single string
            if any(cellfun(@(k) keyCode(KbName(k)), Key2))
                oper = 2; % second key pressed
            end
            
            if any(cellfun(@(k) keyCode(KbName(k)), backKey))
                oper = -1; % Back key pressed
            end
            if oper ~= 0
                break
            end
        end
        checkend;
        while keyIsDown % wait key release
            keyIsDown = KbCheck;
            WaitSecs(0.1);
        end
        WaitSecs(0.1);
    end

end
