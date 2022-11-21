% Wait at least 15 Seconds between Blocks (only after Block 1 has finished)

waitTime = 15;
intervalTime = 1;
timePassed = 0;
printTime = 15;

waitTimeText = ['Please wait for ' num2str(printTime) ' seconds. ...' ...
                ' \n\n ' ...
                ' \n\n Block ' (num2str(BLOCK)) ' will start afterwards.'];

DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
Screen('Flip',ptbWindow);

while timePassed < waitTime
    pause(intervalTime);
    timePassed = timePassed + intervalTime;
    printTime = waitTime - timePassed;
    waitTimeText = ['Please wait for ' num2str(printTime) ' seconds. ...' ...
                    ' \n\n ' ...
                    ' \n\n Block ' (num2str(BLOCK)) ' will start afterwards.'];
    DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
    Screen('Flip',ptbWindow);
end