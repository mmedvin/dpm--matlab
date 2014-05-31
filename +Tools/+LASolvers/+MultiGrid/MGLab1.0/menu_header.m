%MENU_HEADER Generate a menu header
%
%       F_M = MENU_HEADER(PARENT,LABEL,SEPARATOR,ENABLE,BGCOLOR) creates
%       a menu with the given parameters and returns its handle.

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function f_m = menu_header(parent,label,separator,enable,bgcolor)

f_m = uimenu(parent,...
     'Label',label,...
     'Separator', separator,...
     'Enable', enable,...
     'BackGroundColor',bgcolor);
