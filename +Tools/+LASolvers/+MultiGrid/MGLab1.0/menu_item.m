%MENU_ITEM Generate a menu item.
%
%       F_M = MENU_ITEM(PARENT,LABEL,SEPARATOR,ENABLE,BGCOLOR,CALLBACK) creates
%       a menu item with the given parameters and returns its handle.

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function f_m=menu_item(parent,label,separator,enable,bgcolor,callback)

f_m = uimenu(parent,...
     'Label',label,...
     'Separator', separator,...
     'Enable',enable,...
     'BackGroundColor',bgcolor,...
     'Callback', callback);
