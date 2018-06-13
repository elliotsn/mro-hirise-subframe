%
% Function to remove the 2nd token in a string
% and convert the first token to a double.
%
function valout = str1sttok2double(strin)

valout = str2double(strtok(strin,' '));