function checkend

[~,~,kc]=KbCheck;
if kc(KbName('escape'))
    ListenChar(0);
    ShowCursor;
    sca;
    disp("Aborted by ESC Key!!")
    error('Interprated by user; ”√ªß÷’÷π≥Ã–Ú');
end